"""
Discrete Morse Theory for Persistent Homology of Filtered Cosheaves
===================================================================
Bottom-up CoScythe verification on the C4 example.

Key design decisions (matching the paper):
  - Working poset is maintained dynamically (Curry-style in-place update)
  - CoReducePair uses dynamic adjacency for BOTH sigma^+ and tau^-
  - Matching decisions and queue use original K adjacency
  - All n+1 filtration levels updated simultaneously

Expected:  H0: [0, inf), [0, 2)    H1: [3, inf)
"""

import numpy as np
from fractions import Fraction

np.set_printoptions(precision=4, suppress=True, linewidth=120)

# ============================================================
# 1. Simplicial complex K = C4
# ============================================================

VERTICES = [(0,), (1,), (2,), (3,)]
EDGES    = [(0,1), (1,2), (2,3), (0,3)]
ALL      = VERTICES + EDGES
NM       = {(0,): 'a', (1,): 'b', (2,): 'c', (3,): 'd',
            (0,1): 'ab', (1,2): 'bc', (2,3): 'cd', (0,3): 'ad'}
TOL      = 1e-10

def sdim(s):
    return len(s) - 1

def is_face(s, t):
    return sdim(s) + 1 == sdim(t) and set(s).issubset(set(t))

def incidence(s, t):
    """Signed incidence number [s : t]."""
    if not is_face(s, t):
        return 0
    for i in range(len(t)):
        if tuple(t[:i] + t[i+1:]) == s:
            return (-1) ** i
    return 0

# Original K adjacency (static, for matching decisions and queue)
def cofaces_K(eta):
    return [t for t in ALL if is_face(eta, t)]

def faces_K(eta):
    return [s for s in ALL if is_face(s, eta)]

# ============================================================
# 2. Cosheaf data
# ============================================================

def build_cosheaf():
    """Build the 4-level filtered cosheaf on C4."""
    L = 4
    cd = [{s: 0 for s in ALL} for _ in range(L)]
    em = [{} for _ in range(L)]

    # Level 0: vertices a, d
    cd[0][(0,)] = cd[0][(3,)] = 1

    # Level 1: all vertices, edges bc, cd
    for v in VERTICES:
        cd[1][v] = 1
    cd[1][(1,2)] = cd[1][(2,3)] = 1

    # Level 2: all vertices, edges ab, bc, cd
    for v in VERTICES:
        cd[2][v] = 1
    cd[2][(0,1)] = cd[2][(1,2)] = cd[2][(2,3)] = 1

    # Level 3: all simplices; b and bc get F^2
    for s in ALL:
        cd[3][s] = 1
    cd[3][(1,)] = cd[3][(1,2)] = 2

    # Default extension maps: identity where dims match
    for lv in range(L):
        for t in EDGES:
            for s in VERTICES:
                if is_face(s, t):
                    dt, ds = cd[lv][t], cd[lv][s]
                    if dt > 0 and ds > 0 and dt == ds:
                        em[lv][(t, s)] = np.eye(ds)

    # Level 3 overrides
    em[3][((1,2), (1,))] = np.eye(2)               # bc >= b : id_{F^2}
    em[3][((1,2), (2,))] = np.array([[1., 0.]])     # bc >= c : pi_1
    em[3][((0,1), (0,))] = np.eye(1)               # ab >= a : id_F
    em[3][((0,1), (1,))] = np.array([[1.], [0.]])   # ab >= b : iota
    em[3][((2,3), (2,))] = np.eye(1)               # cd >= c
    em[3][((2,3), (3,))] = np.eye(1)               # cd >= d
    em[3][((0,3), (0,))] = np.eye(1)               # ad >= a
    em[3][((0,3), (3,))] = np.eye(1)               # ad >= d

    return cd, em

# ============================================================
# 3. Signed boundary blocks with proper zero-block handling
# ============================================================

def get_signed(lv, tau, sigma, cd, em):
    """Return the CURRENT signed boundary block C^r_{tau,sigma}.

    Returns:
      - None if either costalk dimension is zero (block does not exist)
      - A numpy matrix otherwise (possibly the zero matrix)

    Blocks updated by Schur complement are stored under ('S', tau, sigma).
    Original blocks are stored under (tau, sigma) and multiplied by
    the incidence sign.  If a face relation exists but no block is
    stored, we return the zero matrix (defensive default).
    """
    dt = cd[lv].get(tau, 0)
    ds = cd[lv].get(sigma, 0)
    if dt == 0 or ds == 0:
        return None

    # Schur-complement updated block takes priority
    S = em[lv].get(('S', tau, sigma))
    if S is not None:
        return S

    # Original block
    E = em[lv].get((tau, sigma))
    if E is not None:
        inc = incidence(sigma, tau)
        return float(inc) * E if inc != 0 else np.zeros((ds, dt))

    # No stored block: return zero matrix
    # (covers both non-adjacent pairs and missing-but-adjacent pairs)
    return np.zeros((ds, dt))

def set_signed(lv, tau, sigma, val, em):
    """Store a Schur-complement updated block."""
    em[lv][('S', tau, sigma)] = val

def block_is_nonzero(M, tol=TOL):
    """Check whether a block matrix is present and has nonzero entries."""
    return M is not None and np.any(np.abs(M) > tol)

# ============================================================
# 4. Dynamic adjacency (working poset)
# ============================================================

class DynAdj:
    """Track current nonzero boundary blocks in the working poset.

    This is the data structure corresponding to the covering relations
    of the dynamic working poset in the paper.  A relation sigma < tau
    exists iff at least one filtration level carries a nonzero block.
    """

    def __init__(self, simplices, cd, em, L):
        self.up = {s: set() for s in simplices}
        self.dn = {s: set() for s in simplices}
        # Initialise from nonzero blocks
        for tau in simplices:
            for sigma in simplices:
                if sigma == tau:
                    continue
                for lv in range(L):
                    if block_is_nonzero(get_signed(lv, tau, sigma, cd, em)):
                        self._link(sigma, tau)
                        break

    def _link(self, sigma, tau):
        self.up.setdefault(sigma, set()).add(tau)
        self.dn.setdefault(tau, set()).add(sigma)

    def _unlink(self, sigma, tau):
        self.up.get(sigma, set()).discard(tau)
        self.dn.get(tau, set()).discard(sigma)

    def sync(self, sigma, tau, cd, em, L):
        """Re-check whether the relation sigma < tau is nonzero."""
        alive = False
        for lv in range(L):
            if block_is_nonzero(get_signed(lv, tau, sigma, cd, em)):
                alive = True
                break
        if alive:
            self._link(sigma, tau)
        else:
            self._unlink(sigma, tau)

    def cofaces(self, sigma, working):
        """Current cofaces of sigma in the working poset."""
        return [tau for tau in self.up.get(sigma, set()) if tau in working]

    def faces(self, tau, working):
        """Current faces of tau in the working poset."""
        return [sigma for sigma in self.dn.get(tau, set()) if sigma in working]

    def remove(self, s):
        """Remove a simplex from all adjacency lists."""
        for tau in list(self.up.get(s, set())):
            self.dn.get(tau, set()).discard(s)
        for sigma in list(self.dn.get(s, set())):
            self.up.get(sigma, set()).discard(s)
        self.up.pop(s, None)
        self.dn.pop(s, None)

# ============================================================
# 5. Invertibility check
# ============================================================

def invertible_all_levels(tau, sigma, cd, em, L):
    """Check C^r_{tau,sigma} is square and invertible at every level."""
    for lv in range(L):
        dt = cd[lv].get(tau, 0)
        ds = cd[lv].get(sigma, 0)
        if dt == 0 and ds == 0:
            continue
        if dt == 0 or ds == 0 or dt != ds:
            return False
        C = get_signed(lv, tau, sigma, cd, em)
        if C is None:
            return False
        if abs(np.linalg.det(C)) < TOL:
            return False
    return True

# ============================================================
# 6. CoReducePair (dynamic adjacency for BOTH loops)
# ============================================================

def co_reduce_pair(sigma, tau, cd, em, working, dadj, L):
    """Schur-complement cancellation of the pair (sigma < tau).

    Both sigma^+ and tau^- are taken from the DYNAMIC working poset,
    following CGN's ReducePair which operates entirely on the
    current poset.  Missing blocks are treated as zero.
    """
    s_plus  = [a for a in dadj.cofaces(sigma, working) if a != tau]
    t_minus = [o for o in dadj.faces(tau, working)     if o != sigma]

    for alpha in s_plus:
        for omega in t_minus:
            for lv in range(L):
                da = cd[lv].get(alpha, 0)
                do = cd[lv].get(omega, 0)
                dt = cd[lv].get(tau, 0)
                ds = cd[lv].get(sigma, 0)
                if da == 0 or do == 0 or dt == 0 or ds == 0:
                    continue

                Ctw = get_signed(lv, tau, omega, cd, em)
                Cts = get_signed(lv, tau, sigma, cd, em)
                Cas = get_signed(lv, alpha, sigma, cd, em)

                # Skip if pivot is zero (should not happen for
                # compatible pairs, but defensive)
                if Cts is None or abs(np.linalg.det(Cts)) < TOL:
                    continue
                if Ctw is None:
                    Ctw = np.zeros((do, dt))
                if Cas is None:
                    Cas = np.zeros((ds, da))

                # Schur complement: C_{alpha,omega} -= C_{tau,omega} * C_{tau,sigma}^{-1} * C_{alpha,sigma}
                upd = Ctw @ np.linalg.solve(Cts, Cas)
                Caw = get_signed(lv, alpha, omega, cd, em)
                if Caw is None:
                    Caw = np.zeros((do, da))

                new_val = Caw - upd
                new_val[np.abs(new_val) < TOL] = 0.0
                set_signed(lv, alpha, omega, new_val, em)

            # Sync dynamic adjacency after updating all levels
            dadj.sync(omega, alpha, cd, em, L)

    # Remove sigma and tau from working poset
    working.discard(sigma)
    working.discard(tau)
    dadj.remove(sigma)
    dadj.remove(tau)

# ============================================================
# 7. Bottom-up CoScythe
# ============================================================

def co_scythe(cd, em, L):
    """Bottom-up CoScythe paralleling Scythe from CGN Section 4.

    - Select minimal unclassified simplex as critical root
    - Matching decisions use original K faces
    - Queue propagation uses original K cofaces
    - CoReducePair uses dynamic working poset adjacency
    """
    working      = set(ALL)
    unclassified = set(ALL)
    dadj         = DynAdj(ALL, cd, em, L)
    crit         = set()
    match        = []

    while unclassified:
        # Line 03: minimal unclassified, reverse-lex tie-break
        mn = min(sdim(s) for s in unclassified)
        c = max((s for s in unclassified if sdim(s) == mn),
                key=lambda s: s)
        unclassified.discard(c)
        crit.add(c)
        print(f"  Critical: {NM[c]}")

        # Lines 05-06: BFS queue with flags
        que = [c]
        flagged = {c}

        while que:
            y = que.pop(0)

            # Skip if already removed from working poset
            if y not in working:
                continue

            # Line 09: check original K faces for unique unclassified
            all_uncl = [s for s in faces_K(y) if s in unclassified]
            if len(all_uncl) == 1:
                sigma = all_uncl[0]
                if invertible_all_levels(y, sigma, cd, em, L):
                    # Line 10: match
                    match.append((sigma, y))
                    print(f"  Matched: ({NM[sigma]} < {NM[y]})")

                    # Line 11: enqueue sigma^+ \ {y} from original K
                    for a in cofaces_K(sigma):
                        if a != y and a in working and a not in flagged:
                            que.append(a)
                            flagged.add(a)

                    # Line 12: Schur-complement reduction
                    co_reduce_pair(sigma, y, cd, em, working, dadj, L)
                    unclassified.discard(sigma)
                    unclassified.discard(y)

            # Line 14: enqueue y^+ from original K
            for a in cofaces_K(y):
                if a in working and a not in flagged:
                    que.append(a)
                    flagged.add(a)

    print(f"  Critical simplices: {sorted(NM[s] for s in crit)}")
    return match, crit

# ============================================================
# 8. Persistence computation
# ============================================================

def build_persistence_matrix(cd, em, simplices, L):
    """Build persistence boundary matrix at the top filtration level.

    Each basis vector corresponds to a coordinate of a costalk that
    first appears at some filtration level (its birth time).
    """
    top = L - 1
    basis = []  # (simplex, vector_index, birth_level, homological_dim)

    for s in simplices:
        d_top = cd[top].get(s, 0)
        if d_top == 0:
            continue
        prev = 0
        for lv in range(L):
            cur = cd[lv].get(s, 0)
            for vi in range(prev, cur):
                basis.append((s, vi, lv, sdim(s)))
            prev = cur

    # Sort by (birth, dimension, simplex, vector index)
    basis.sort(key=lambda x: (x[2], x[3], x[0], x[1]))
    n = len(basis)
    births = [b[2] for b in basis]
    dims   = [b[3] for b in basis]

    print(f"\n  Basis ({n} vectors):")
    for i, (s, vi, br, hd) in enumerate(basis):
        print(f"    {i:2d}: {NM[s]:>3s}[{vi}] dim={hd} birth={br}")

    # Build boundary matrix at top level
    D = np.zeros((n, n))
    for j, (tau, vj, _, dj) in enumerate(basis):
        if dj == 0:
            continue
        for i, (sigma, vi, _, di) in enumerate(basis):
            if di != dj - 1:
                continue
            blk = get_signed(top, tau, sigma, cd, em)
            if blk is None:
                continue
            if vi < blk.shape[0] and vj < blk.shape[1]:
                D[i, j] = blk[vi, vj]

    return D, births, dims

def column_reduce(D, births, dims):
    """Standard persistence column reduction."""
    n = D.shape[1]
    R = D.copy()

    def low(j):
        nz = np.where(np.abs(R[:, j]) > TOL)[0]
        return int(nz[-1]) if len(nz) else -1

    piv = {}
    for j in range(n):
        while True:
            l = low(j)
            if l < 0 or l not in piv:
                break
            r = piv[l]
            R[:, j] -= (R[l, j] / R[l, r]) * R[:, r]
        l = low(j)
        if l >= 0:
            piv[l] = j

    # Extract barcodes
    paired = set()
    bars = {}
    for j in range(n):
        l = low(j)
        if l >= 0:
            b, d = births[l], births[j]
            k = dims[l]
            if b != d:
                bars.setdefault(k, []).append((b, d))
            paired.update([l, j])
    for j in range(n):
        if j not in paired and low(j) < 0:
            bars.setdefault(dims[j], []).append((births[j], float('inf')))
    return bars

# ============================================================
# 9. Display and verification
# ============================================================

def fmt(bars):
    out = []
    for k in sorted(bars):
        if bars[k]:
            pieces = [f"[{b},{'inf' if d == float('inf') else int(d)})"
                      for b, d in sorted(bars[k])]
            out.append(f"  H{k}: {', '.join(pieces)}")
    return '\n'.join(out)

def main():
    L = 4
    print("=" * 64)
    print(" Cosheaf Morse Theory - C4 Verification (Bottom-Up)")
    print("=" * 64)

    # --- Original barcodes ---
    print("\n[1] Original barcodes")
    cd1, em1 = build_cosheaf()
    D1, b1, d1 = build_persistence_matrix(cd1, em1, VERTICES + EDGES, L)
    bars_orig = column_reduce(D1, b1, d1)
    print(f"\n{fmt(bars_orig)}")

    # --- CoScythe ---
    print("\n[2] Bottom-Up CoScythe")
    cd2, em2 = build_cosheaf()
    match, crit = co_scythe(cd2, em2, L)

    # --- Morse barcodes (all dimensions, not just 0 and 1) ---
    print("\n[3] Morse barcodes")
    max_dim = max(sdim(s) for s in crit) if crit else 0
    morse_simps = []
    for d in range(max_dim + 1):
        morse_simps += sorted([s for s in crit if sdim(s) == d])

    D2, b2, d2 = build_persistence_matrix(cd2, em2, morse_simps, L)
    print(f"\n  Morse boundary:\n{D2}")
    bars_morse = column_reduce(D2, b2, d2)
    print(f"\n{fmt(bars_morse)}")

    # --- Verification ---
    print("\n[4] Verification")
    ok = True
    for k in sorted(set(list(bars_orig) + list(bars_morse))):
        a = sorted(bars_orig.get(k, []))
        b = sorted(bars_morse.get(k, []))
        same = (a == b)
        ok &= same
        print(f"  H{k}: orig={a} morse={b} {'OK' if same else 'MISMATCH'}")

    print(f"\n  {'ALL BARCODES MATCH!' if ok else 'MISMATCH DETECTED'}")

    # --- Summary ---
    print(f"\n  Matching: {{{', '.join(f'({NM[s]}<{NM[t]})' for s, t in match)}}}")
    print(f"  Critical: {sorted(NM[s] for s in crit)}")
    print(f"  Reduction: {len(ALL)} -> {len(crit)} simplices")

    # --- Assert expected values ---
    expected = {0: [(0, 2), (0, float('inf'))], 1: [(3, float('inf'))]}
    for k in sorted(set(list(expected) + list(bars_morse))):
        got = sorted(bars_morse.get(k, []))
        exp = sorted(expected.get(k, []))
        assert got == exp, f"H{k}: expected {exp}, got {got}"
    print("  Expected barcode assertion passed.")

if __name__ == "__main__":
    main()