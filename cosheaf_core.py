from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Any
import json
import math
import random
import time
from collections import deque, defaultdict

import sympy as sp

Simplex = Tuple[Any, ...]
Matrix = sp.Matrix


def sorted_simplex(vertices: Sequence[Any]) -> Simplex:
    return tuple(sorted(vertices))


def simplex_dim(s: Simplex) -> int:
    return len(s) - 1


def codim1_faces(s: Simplex) -> List[Simplex]:
    if len(s) <= 1:
        return []
    return [s[:i] + s[i+1:] for i in range(len(s))]


def incidence_sign(face: Simplex, simplex: Simplex) -> int:
    if len(face) + 1 != len(simplex):
        return 0
    for i in range(len(simplex)):
        if simplex[:i] + simplex[i+1:] == face:
            return 1 if i % 2 == 0 else -1
    return 0


def to_matrix(data) -> Matrix:
    if isinstance(data, Matrix):
        return data
    return Matrix(data)


def zero_matrix(rows: int, cols: int) -> Matrix:
    return sp.zeros(rows, cols)


def identity_matrix(n: int) -> Matrix:
    return sp.eye(n)


def basis_matrix(vectors: List[Matrix], rows: int) -> Matrix:
    if not vectors:
        return sp.zeros(rows, 0)
    cols = []
    for v in vectors:
        mv = to_matrix(v)
        if mv.cols != 1:
            raise ValueError("Basis vector must be a column vector")
        cols.append(mv)
    return Matrix.hstack(*cols)


def matrix_rank(M: Matrix) -> int:
    return int(to_matrix(M).rank())


def invertible(M: Matrix) -> bool:
    M = to_matrix(M)
    return M.rows == M.cols and M.rank() == M.rows


def ensure_invertible(M: Matrix) -> Matrix:
    M = to_matrix(M)
    if not invertible(M):
        raise ValueError(f"Matrix not invertible: {M}")
    return M


@dataclass
class FilteredCosheaf:
    simplices: List[Simplex]
    levels: int  # number of filtration levels = n+1
    costalk_dims: List[Dict[Simplex, int]]
    extension_maps: List[Dict[Tuple[Simplex, Simplex], Matrix]]  # only codim-1 coface relations tau>=sigma
    level_maps: List[Dict[Simplex, Matrix]]  # maps from level i to i+1 on each simplex, length levels-1
    name: str = "example"

    simplices_by_dim: Dict[int, List[Simplex]] = field(init=False)
    max_dim: int = field(init=False)

    def __post_init__(self):
        self.simplices = [tuple(s) for s in self.simplices]
        self.simplices_by_dim = defaultdict(list)
        for s in self.simplices:
            self.simplices_by_dim[simplex_dim(s)].append(s)
        for k in self.simplices_by_dim:
            self.simplices_by_dim[k] = sorted(self.simplices_by_dim[k])
        self.max_dim = max(self.simplices_by_dim) if self.simplices_by_dim else -1
        assert len(self.costalk_dims) == self.levels
        assert len(self.extension_maps) == self.levels
        assert len(self.level_maps) == max(0, self.levels - 1)
        self._validate_shapes()
        self._validate_naturality()

    def _validate_shapes(self):
        for r in range(self.levels):
            dims = self.costalk_dims[r]
            for s in self.simplices:
                dims.setdefault(s, 0)
            for tau in self.simplices:
                for sigma in codim1_faces(tau):
                    key = (tau, sigma)
                    if key not in self.extension_maps[r]:
                        rows = dims[sigma]
                        cols = dims[tau]
                        self.extension_maps[r][key] = sp.zeros(rows, cols)
                    else:
                        M = to_matrix(self.extension_maps[r][key])
                        if M.rows != dims[sigma] or M.cols != dims[tau]:
                            raise ValueError(f"Bad shape at level {r} for {tau}>={sigma}: {M.shape}, expected {(dims[sigma], dims[tau])}")
                        self.extension_maps[r][key] = M
        for i in range(self.levels - 1):
            for s in self.simplices:
                if s not in self.level_maps[i]:
                    self.level_maps[i][s] = sp.zeros(self.costalk_dims[i+1][s], self.costalk_dims[i][s])
                else:
                    M = to_matrix(self.level_maps[i][s])
                    if M.rows != self.costalk_dims[i+1][s] or M.cols != self.costalk_dims[i][s]:
                        raise ValueError(f"Bad level-map shape at {i} simplex {s}: {M.shape}")
                    self.level_maps[i][s] = M

    def _validate_naturality(self):
        # only codim-1 relations needed since all others omitted in paper/code
        for i in range(self.levels - 1):
            for tau in self.simplices:
                for sigma in codim1_faces(tau):
                    left = self.extension_maps[i+1][(tau, sigma)] * self.level_maps[i][tau]
                    right = self.level_maps[i][sigma] * self.extension_maps[i][(tau, sigma)]
                    if left != right:
                        raise ValueError(f"Naturality failed at level {i}, relation {tau}>={sigma}\nleft={left}\nright={right}")

    def signed_block(self, level: int, alpha: Simplex, beta: Simplex) -> Matrix:
        sign = incidence_sign(beta, alpha)
        if sign == 0:
            return sp.zeros(self.costalk_dims[level][beta], self.costalk_dims[level][alpha])
        return sign * self.extension_maps[level][(alpha, beta)]

    def chain_group_dim(self, level: int, k: int) -> int:
        return sum(self.costalk_dims[level][s] for s in self.simplices_by_dim.get(k, []))

    def simplex_offsets(self, level: int, k: int) -> Dict[Simplex, Tuple[int, int]]:
        offs = {}
        cur = 0
        for s in self.simplices_by_dim.get(k, []):
            d = self.costalk_dims[level][s]
            offs[s] = (cur, cur + d)
            cur += d
        return offs

    def boundary_matrix(self, level: int, k: int) -> Matrix:
        # D_k: C_k -> C_{k-1}
        domain = self.simplices_by_dim.get(k, [])
        target = self.simplices_by_dim.get(k-1, [])
        rows = self.chain_group_dim(level, k-1)
        cols = self.chain_group_dim(level, k)
        M = sp.zeros(rows, cols)
        col_off = self.simplex_offsets(level, k)
        row_off = self.simplex_offsets(level, k-1)
        for alpha in domain:
            c0, c1 = col_off[alpha]
            for beta in target:
                block = self.signed_block(level, alpha, beta)
                if block.rows == 0 or block.cols == 0 or block == sp.zeros(block.rows, block.cols):
                    continue
                r0, r1 = row_off[beta]
                M[r0:r1, c0:c1] = block
        return M

    def chain_map_matrix(self, i: int, k: int) -> Matrix:
        # from level i to i+1 on C_k
        src = self.simplices_by_dim.get(k, [])
        dst = self.simplices_by_dim.get(k, [])
        rows = self.chain_group_dim(i+1, k)
        cols = self.chain_group_dim(i, k)
        M = sp.zeros(rows, cols)
        src_off = self.simplex_offsets(i, k)
        dst_off = self.simplex_offsets(i+1, k)
        for s in src:
            c0, c1 = src_off[s]
            r0, r1 = dst_off[s]
            block = self.level_maps[i][s]
            M[r0:r1, c0:c1] = block
        return M

    def composite_chain_map(self, i: int, j: int, k: int) -> Matrix:
        if i == j:
            return sp.eye(self.chain_group_dim(i, k))
        M = self.chain_map_matrix(i, k)
        for t in range(i + 1, j):
            M = self.chain_map_matrix(t, k) * M
        return M

    def births(self) -> Dict[Simplex, Optional[int]]:
        out = {}
        for s in self.simplices:
            b = None
            for r in range(self.levels):
                if self.costalk_dims[r][s] > 0:
                    b = r
                    break
            out[s] = b
        return out


def boundary_data(fc: FilteredCosheaf):
    D = []
    for r in range(fc.levels):
        levelD = {}
        for k in range(fc.max_dim + 1):
            levelD[k] = fc.boundary_matrix(r, k)
        D.append(levelD)
    return D


def homology_persistent_ranks(fc: FilteredCosheaf, use_morse: bool = False, morse_data=None):
    if use_morse:
        if morse_data is None:
            raise ValueError("Need morse_data")
        dims = morse_data['dims']
        boundaries = morse_data['boundaries']
        chain_maps = morse_data['chain_maps']
        max_dim = morse_data['max_dim']
        levels = morse_data['levels']
    else:
        max_dim = fc.max_dim
        levels = fc.levels
        boundaries = [{k: fc.boundary_matrix(r, k) for k in range(max_dim + 1)} for r in range(levels)]
        chain_maps = [{k: fc.chain_map_matrix(r, k) for k in range(max_dim + 1)} for r in range(levels - 1)]
        dims = [{k: fc.chain_group_dim(r, k) for k in range(max_dim + 1)} for r in range(levels)]

    Z_basis = [[None]*(max_dim+1) for _ in range(levels)]
    B_basis = [[None]*(max_dim+1) for _ in range(levels)]
    for r in range(levels):
        for k in range(max_dim + 1):
            Dk = boundaries[r][k]
            Z_basis[r][k] = basis_matrix(Dk.nullspace(), dims[r][k])
            Dkp1 = boundaries[r].get(k+1, sp.zeros(dims[r][k], 0))
            B_basis[r][k] = basis_matrix(Dkp1.columnspace(), dims[r][k])

    def composite(i: int, j: int, k: int) -> Matrix:
        if i == j:
            return sp.eye(dims[i][k])
        M = chain_maps[i][k]
        for t in range(i + 1, j):
            M = chain_maps[t][k] * M
        return M

    rho = {(k, i, j): 0 for k in range(max_dim + 1) for i in range(levels) for j in range(i, levels)}
    for k in range(max_dim + 1):
        for i in range(levels):
            Zi = Z_basis[i][k]
            for j in range(i, levels):
                Fij = composite(i, j, k)
                img_cycles = Fij * Zi if Zi.cols else sp.zeros(dims[j][k], 0)
                Bj = B_basis[j][k]
                if Bj.cols == 0:
                    rho[(k, i, j)] = matrix_rank(img_cycles)
                else:
                    rho[(k, i, j)] = matrix_rank(Matrix.hstack(Bj, img_cycles)) - matrix_rank(Bj)
    return rho


def barcodes_from_ranks(levels: int, rho: Dict[Tuple[int, int, int], int], max_dim: int):
    bars: Dict[int, List[Tuple[int, Optional[int]]]] = {k: [] for k in range(max_dim + 1)}
    def r(k, i, j):
        if i < 0 or j < 0 or i >= levels or j >= levels or i > j:
            return 0
        return int(rho.get((k, i, j), 0))

    for k in range(max_dim + 1):
        # finite intervals [b,d)
        for b in range(levels):
            for d in range(b + 1, levels):
                mult = r(k, b, d-1) - r(k, b-1, d-1) - r(k, b, d) + r(k, b-1, d)
                for _ in range(max(0, mult)):
                    bars[k].append((b, d))
            inf_mult = r(k, b, levels-1) - r(k, b-1, levels-1)
            # subtract finite intervals already counted with same birth b extending to finite d
            finite_with_birth = sum(1 for bb, dd in bars[k] if bb == b)
            # This subtraction is wrong if other births counted. Use formula directly for infinite part below after separate pass.
        # redo properly using direct formula without side-effects
        bars[k] = []
        for b in range(levels):
            for d in range(b + 1, levels):
                mult = r(k, b, d-1) - r(k, b-1, d-1) - r(k, b, d) + r(k, b-1, d)
                for _ in range(max(0, mult)):
                    bars[k].append((b, d))
            inf_mult = r(k, b, levels-1) - r(k, b-1, levels-1)
            for _ in range(max(0, inf_mult)):
                bars[k].append((b, None))
        bars[k].sort(key=lambda x: (x[0], float('inf') if x[1] is None else x[1]))
    return bars


def format_barcodes(bars: Dict[int, List[Tuple[int, Optional[int]]]]) -> Dict[str, List[str]]:
    out = {}
    for k, lst in bars.items():
        if not lst:
            continue
        out[f"H{k}"] = [f"[{b},∞)" if d is None else f"[{b},{d})" for b, d in lst]
    return out


@dataclass
class CoScytheResult:
    matching: List[Tuple[Simplex, Simplex]]
    critical: List[Simplex]
    morse_boundaries: List[Dict[int, Matrix]]
    morse_chain_maps: List[Dict[int, Matrix]]
    dims: List[Dict[int, int]]
    critical_by_dim: Dict[int, List[Simplex]]
    trace: List[dict]
    timings: Dict[str, float]


def run_coscythe(fc: FilteredCosheaf, tie_breaker: str = "minlex", record_trace: bool = False) -> CoScytheResult:
    t0 = time.perf_counter()
    # current simplices in working poset
    present = set(fc.simplices)
    critical = set()
    matching: List[Tuple[Simplex, Simplex]] = []
    matched = set()
    births = fc.births()

    # dynamic adjacency from initial simplicial complex
    faces = {s: set(codim1_faces(s)) for s in fc.simplices}
    cofaces = {s: set() for s in fc.simplices}
    for tau in fc.simplices:
        for sigma in codim1_faces(tau):
            cofaces[sigma].add(tau)

    # boundary blocks at each level, for all current codim1 adjacency pairs
    blocks: List[Dict[Tuple[Simplex, Simplex], Matrix]] = []
    for r in range(fc.levels):
        bd = {}
        for alpha in fc.simplices:
            for beta in codim1_faces(alpha):
                bd[(alpha, beta)] = fc.signed_block(r, alpha, beta)
        blocks.append(bd)

    trace = []

    def current_unclassified() -> List[Simplex]:
        return [s for s in present if s not in critical and s not in matched]

    def choose_minimal(cands: List[Simplex]) -> Simplex:
        mins = []
        min_dim = min(simplex_dim(s) for s in cands)
        mins = [s for s in cands if simplex_dim(s) == min_dim]
        if tie_breaker == "random":
            return random.choice(mins)
        if tie_breaker == "reverselex":
            return sorted(mins, reverse=True)[0]
        return sorted(mins)[0]

    def is_compatible_pair(sigma: Simplex, tau: Simplex) -> bool:
        # quick birth-time rejection
        if births[sigma] is None or births[tau] is None:
            return False
        if births[sigma] != births[tau]:
            return False
        for r in range(fc.levels):
            M = blocks[r].get((tau, sigma), sp.zeros(fc.costalk_dims[r][sigma], fc.costalk_dims[r][tau]))
            if not invertible(M):
                return False
        return True

    def coreduce_pair(sigma: Simplex, tau: Simplex):
        # insert dynamic relations and update blocks
        sigma_cofaces = [a for a in list(cofaces.get(sigma, set())) if a in present and a != tau]
        tau_faces = [w for w in list(faces.get(tau, set())) if w in present and w != sigma]
        for alpha in sigma_cofaces:
            for omega in tau_faces:
                # add dynamic adjacency relation omega < alpha
                faces.setdefault(alpha, set()).add(omega)
                cofaces.setdefault(omega, set()).add(alpha)
                for r in range(fc.levels):
                    A = blocks[r].get((tau, sigma), sp.zeros(fc.costalk_dims[r][sigma], fc.costalk_dims[r][tau]))
                    T_omega = blocks[r].get((tau, omega), sp.zeros(fc.costalk_dims[r][omega], fc.costalk_dims[r][tau]))
                    A_sigma = blocks[r].get((alpha, sigma), sp.zeros(fc.costalk_dims[r][sigma], fc.costalk_dims[r][alpha]))
                    old = blocks[r].get((alpha, omega), sp.zeros(fc.costalk_dims[r][omega], fc.costalk_dims[r][alpha]))
                    new = old - T_omega * ensure_invertible(A).inv() * A_sigma
                    blocks[r][(alpha, omega)] = new
        # remove sigma and tau from working poset and adjacency
        present.discard(sigma)
        present.discard(tau)
        matched.add(sigma)
        matched.add(tau)
        for x in list(faces.get(sigma, set())):
            cofaces[x].discard(sigma)
        for x in list(cofaces.get(sigma, set())):
            faces[x].discard(sigma)
        for x in list(faces.get(tau, set())):
            cofaces[x].discard(tau)
        for x in list(cofaces.get(tau, set())):
            faces[x].discard(tau)
        faces[sigma] = set(); cofaces[sigma] = set()
        faces[tau] = set(); cofaces[tau] = set()

    while current_unclassified():
        c = choose_minimal(current_unclassified())
        critical.add(c)
        que = deque([c])
        flagged = {c}
        if record_trace:
            trace.append({"action": "mark_critical", "simplex": str(c), "critical": [str(x) for x in sorted(critical)], "matching": [[str(a), str(b)] for a,b in matching], "queue": [str(x) for x in que]})
        while que:
            y = que.popleft()
            if y not in present:
                continue
            unclassified_faces = [s for s in faces.get(y, set()) if s in present and s not in critical and s not in matched]
            paired = False
            if len(unclassified_faces) == 1:
                sigma = unclassified_faces[0]
                if is_compatible_pair(sigma, y):
                    matching.append((sigma, y))
                    # enqueue cofaces of sigma other than y before reduction
                    for a in sorted(cofaces.get(sigma, set())):
                        if a != y and a in present and a not in flagged:
                            que.append(a); flagged.add(a)
                    if record_trace:
                        trace.append({"action": "pair", "sigma": str(sigma), "tau": str(y), "critical": [str(x) for x in sorted(critical)], "matching": [[str(a), str(b)] for a,b in matching], "queue_before_reduce": [str(x) for x in que]})
                    coreduce_pair(sigma, y)
                    paired = True
            # regardless, enqueue cofaces of y (still present)
            for a in sorted(cofaces.get(y, set())):
                if a in present and a not in flagged:
                    que.append(a); flagged.add(a)
            if record_trace:
                trace.append({"action": "dequeue", "simplex": str(y), "paired": paired, "queue": [str(x) for x in que], "critical": [str(x) for x in sorted(critical)], "matching": [[str(a), str(b)] for a,b in matching]})

    critical_list = sorted(critical)
    critical_by_dim = defaultdict(list)
    for s in critical_list:
        critical_by_dim[simplex_dim(s)].append(s)

    morse_boundaries: List[Dict[int, Matrix]] = []
    dims: List[Dict[int, int]] = []
    offsets: List[Dict[int, Dict[Simplex, Tuple[int,int]]]] = []
    for r in range(fc.levels):
        level_dims = {}
        level_offsets = {}
        for k in range(fc.max_dim + 1):
            cur = 0
            level_offsets[k] = {}
            for s in critical_by_dim.get(k, []):
                d = fc.costalk_dims[r][s]
                level_offsets[k][s] = (cur, cur + d)
                cur += d
            level_dims[k] = cur
        dims.append(level_dims)
        offsets.append(level_offsets)
        level_bd = {}
        for k in range(fc.max_dim + 1):
            M = sp.zeros(level_dims.get(k-1, 0), level_dims.get(k, 0))
            for alpha in critical_by_dim.get(k, []):
                c0,c1 = level_offsets[k][alpha]
                for omega in critical_by_dim.get(k-1, []):
                    block = blocks[r].get((alpha, omega), sp.zeros(fc.costalk_dims[r][omega], fc.costalk_dims[r][alpha]))
                    if block.rows == 0 or block.cols == 0:
                        continue
                    r0,r1 = level_offsets[k-1][omega]
                    M[r0:r1, c0:c1] = block
            level_bd[k] = M
        morse_boundaries.append(level_bd)

    # Morse chain maps between levels: restrict original level maps to critical simplices
    morse_chain_maps: List[Dict[int, Matrix]] = []
    for i in range(fc.levels - 1):
        level_maps_by_k = {}
        for k in range(fc.max_dim + 1):
            rows = dims[i+1][k]
            cols = dims[i][k]
            M = sp.zeros(rows, cols)
            for s in critical_by_dim.get(k, []):
                c0,c1 = offsets[i][k][s]
                r0,r1 = offsets[i+1][k][s]
                M[r0:r1, c0:c1] = fc.level_maps[i][s]
            level_maps_by_k[k] = M
        morse_chain_maps.append(level_maps_by_k)

    timings = {"coscythe_seconds": time.perf_counter() - t0}
    return CoScytheResult(matching=matching, critical=critical_list, morse_boundaries=morse_boundaries, morse_chain_maps=morse_chain_maps, dims=dims, critical_by_dim=dict(critical_by_dim), trace=trace, timings=timings)


def morse_persistent_barcodes(fc: FilteredCosheaf, result: CoScytheResult):
    rho = homology_persistent_ranks(
        fc,
        use_morse=True,
        morse_data={
            'levels': fc.levels,
            'max_dim': fc.max_dim,
            'boundaries': result.morse_boundaries,
            'chain_maps': result.morse_chain_maps,
            'dims': result.dims,
        },
    )
    return barcodes_from_ranks(fc.levels, rho, fc.max_dim), rho


def original_persistent_barcodes(fc: FilteredCosheaf):
    rho = homology_persistent_ranks(fc)
    return barcodes_from_ranks(fc.levels, rho, fc.max_dim), rho


# ---- Example builders ----

def _zero_extended_identity_filtration(simplices: List[Simplex], births: Dict[Simplex, int], levels: int, name: str) -> FilteredCosheaf:
    costalk_dims = []
    ext_maps = []
    for r in range(levels):
        dims = {}
        ext = {}
        for s in simplices:
            dims[s] = 1 if births[s] <= r else 0
        for tau in simplices:
            for sigma in codim1_faces(tau):
                if dims[tau] == dims[sigma] == 1:
                    ext[(tau, sigma)] = sp.eye(1)
                else:
                    ext[(tau, sigma)] = sp.zeros(dims[sigma], dims[tau])
        costalk_dims.append(dims)
        ext_maps.append(ext)
    level_maps = []
    for r in range(levels - 1):
        mp = {}
        for s in simplices:
            d0 = costalk_dims[r][s]
            d1 = costalk_dims[r+1][s]
            if d0 == d1 == 1:
                mp[s] = sp.eye(1)
            else:
                mp[s] = sp.zeros(d1, d0)
        level_maps.append(mp)
    return FilteredCosheaf(simplices=simplices, levels=levels, costalk_dims=costalk_dims, extension_maps=ext_maps, level_maps=level_maps, name=name)


def build_c4_nonconstant_example() -> FilteredCosheaf:
    simplices = [
        ('a',), ('b',), ('c',), ('d',),
        ('a','b'), ('a','d'), ('b','c'), ('c','d')
    ]
    levels = 4
    costalk_dims = []
    ext_maps = []
    nonzero_levels = [
        {('a',), ('d',)},
        {('a',), ('b',), ('c',), ('d',), ('b','c'), ('c','d')},
        {('a',), ('b',), ('c',), ('d',), ('a','b'), ('b','c'), ('c','d')},
        set(simplices),
    ]
    for r in range(levels):
        dims = {s: 0 for s in simplices}
        for s in nonzero_levels[r]:
            dims[s] = 1
        if r == 3:
            dims[('b',)] = 2
            dims[('b','c')] = 2
        ext = {}
        for tau in simplices:
            for sigma in codim1_faces(tau):
                rows, cols = dims[sigma], dims[tau]
                M = sp.zeros(rows, cols)
                if rows == 0 or cols == 0:
                    ext[(tau, sigma)] = M
                    continue
                if r < 3:
                    M = sp.eye(1)
                else:
                    if tau == ('b','c') and sigma == ('b',):
                        M = sp.eye(2)
                    elif tau == ('b','c') and sigma == ('c',):
                        M = Matrix([[1, 0]])
                    elif tau == ('a','b') and sigma == ('b',):
                        M = Matrix([[1], [0]])
                    else:
                        if rows == cols:
                            M = sp.eye(rows)
                ext[(tau, sigma)] = M
        costalk_dims.append(dims)
        ext_maps.append(ext)
    level_maps = []
    for r in range(levels - 1):
        mp = {}
        for s in simplices:
            d0 = costalk_dims[r][s]
            d1 = costalk_dims[r+1][s]
            M = sp.zeros(d1, d0)
            if d0 == d1 and d0 > 0:
                M = sp.eye(d0)
            elif d0 == 1 and d1 == 2 and s in {('b',), ('b','c')}:
                M = Matrix([[1], [0]])
            mp[s] = M
        level_maps.append(mp)
    return FilteredCosheaf(simplices=simplices, levels=levels, costalk_dims=costalk_dims, extension_maps=ext_maps, level_maps=level_maps, name="c4_nonconstant")


def build_path_large_reduction_example(length: int = 4) -> FilteredCosheaf:
    vertices = [(i,) for i in range(length)]
    edges = [(i, i+1) for i in range(length-1)]
    simplices = vertices + edges
    births = {s: 0 for s in simplices}
    return _zero_extended_identity_filtration(simplices, births, levels=2, name=f"path_{length}_large_reduction")


def build_filled_triangle_no_reduction_example() -> FilteredCosheaf:
    simplices = [(0,), (1,), (2,), (0,1), (0,2), (1,2), (0,1,2)]
    births = {(0,):0, (1,):0, (2,):0, (0,1):1, (0,2):1, (1,2):1, (0,1,2):2}
    return _zero_extended_identity_filtration(simplices, births, levels=3, name="filled_triangle_no_reduction")


def build_c4_zero_extended_gallery_example() -> FilteredCosheaf:
    simplices = [('a',), ('b',), ('c',), ('d',), ('a','b'), ('a','d'), ('b','c'), ('c','d')]
    births = {('a',):0, ('d',):0, ('b',):1, ('c',):1, ('b','c'):1, ('c','d'):1, ('a','b'):2, ('a','d'):3}
    return _zero_extended_identity_filtration(simplices, births, levels=4, name="c4_gallery_moderate")


def build_random_path_example(num_vertices: int, levels: int = 3, costalk_dim: int = 1, seed: int = 0) -> FilteredCosheaf:
    # For complexity experiments: all simplices present at all levels, constant costalk dimension d,
    # and level maps equal to the identity so that naturality is automatic.
    vertices = [(i,) for i in range(num_vertices)]
    edges = [(i, i+1) for i in range(num_vertices-1)]
    simplices = vertices + edges
    base_ext = {}
    for tau in simplices:
        for sigma in codim1_faces(tau):
            A = sp.eye(costalk_dim)
            if costalk_dim > 1:
                for i in range(costalk_dim-1):
                    A[i, i+1] = (seed + len(tau) + len(sigma) + i + 1) % 3
            base_ext[(tau, sigma)] = A
    costalk_dims = [{s: costalk_dim for s in simplices} for _ in range(levels)]
    ext_maps = [{k: v.copy() for k, v in base_ext.items()} for _ in range(levels)]
    level_maps = [{s: sp.eye(costalk_dim) for s in simplices} for _ in range(levels-1)]
    return FilteredCosheaf(simplices=simplices, levels=levels, costalk_dims=costalk_dims, extension_maps=ext_maps, level_maps=level_maps, name=f"random_path_v{num_vertices}_d{costalk_dim}")


def summary_dict(fc: FilteredCosheaf, result: CoScytheResult, bars_orig, bars_morse):
    return {
        "name": fc.name,
        "levels": fc.levels,
        "simplices": [list(s) for s in fc.simplices],
        "matching": [[list(a), list(b)] for a, b in result.matching],
        "critical": [list(s) for s in result.critical],
        "critical_count": len(result.critical),
        "reduction_ratio": len(result.critical) / len(fc.simplices) if fc.simplices else 0.0,
        "barcodes_original": format_barcodes(bars_orig),
        "barcodes_morse": format_barcodes(bars_morse),
        "timings": result.timings,
    }


def write_json(path: str, obj):
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)


def bars_equal(b1, b2) -> bool:
    return format_barcodes(b1) == format_barcodes(b2)
