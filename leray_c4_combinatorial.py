from __future__ import annotations

import json
from pathlib import Path
import numpy as np

VERTICES = ['a', 'b', 'c', 'd']
EDGES = ['ab', 'bc', 'cd', 'ad']
SIMPLICES = VERTICES + EDGES
DIM = {s: 0 if len(s) == 1 else 1 for s in SIMPLICES}
STAR = {
    'a':  {0, 1, 2, 3, 13, 14, 15},
    'b':  {4, 1, 2, 3, 5, 6, 7},
    'c':  {8, 5, 6, 7, 9, 10, 11},
    'd':  {12, 9, 10, 11, 13, 14, 15},
    'ab': {1, 2, 3},
    'bc': {5, 6, 7},
    'cd': {9, 10, 11},
    'ad': {13, 14, 15},
}
X = {
    0: {0, 12},
    1: {0, 4, 5, 6, 7, 8, 9, 10, 11, 12},
    2: {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
    3: set(range(16)),
}
ORDER = {'a': 0, 'b': 1, 'c': 2, 'd': 3}


def incidence(face: str, edge: str) -> int:
    if DIM[edge] != DIM[face] + 1:
        return 0
    if face not in edge:
        return 0
    x, y = edge[0], edge[1]
    return +1 if face == y else -1


def num_components(S: set[int]) -> int:
    if not S:
        return 0
    points = sorted(S)
    n = len(points)
    if n == 16:
        return 1
    gaps = 0
    for i in range(n):
        nxt = points[(i + 1) % n]
        cur = points[i]
        gap = (points[0] - cur) % 16 if i == n - 1 else nxt - cur
        if gap > 1:
            gaps += 1
    return gaps if gaps > 0 else 1


def components_of(S: set[int]) -> list[set[int]]:
    if not S:
        return []
    points = sorted(S)
    n = len(points)
    if n == 16:
        return [set(points)]
    start = 0
    for i in range(n):
        nxt = points[(i + 1) % n]
        cur = points[i]
        gap = (points[0] - cur) % 16 if i == n - 1 else nxt - cur
        if gap > 1:
            start = (i + 1) % n
            break
    else:
        return [set(points)]
    comps: list[set[int]] = []
    cur_comp = [points[start]]
    for k in range(1, n):
        idx = (start + k) % n
        prev_idx = (start + k - 1) % n
        if (points[idx] - points[prev_idx]) % 16 == 1:
            cur_comp.append(points[idx])
        else:
            comps.append(set(cur_comp))
            cur_comp = [points[idx]]
    comps.append(set(cur_comp))
    return comps


def costalk_dim(sigma: str, level: int) -> int:
    return num_components(STAR[sigma] & X[level])


def extension_matrix(sigma: str, tau: str, level: int) -> np.ndarray:
    src = components_of(STAR[tau] & X[level])
    tgt = components_of(STAR[sigma] & X[level])
    M = np.zeros((len(tgt), len(src)), dtype=int)
    for k, comp_src in enumerate(src):
        for j, comp_tgt in enumerate(tgt):
            if comp_src.issubset(comp_tgt):
                M[j, k] = 1
                break
        else:
            raise ValueError(f"component of st({tau}) not contained in any component of st({sigma}) at level {level}")
    return M


def boundary_matrix(level: int) -> np.ndarray:
    v_dims = [costalk_dim(v, level) for v in VERTICES]
    e_dims = [costalk_dim(e, level) for e in EDGES]
    rows, cols = sum(v_dims), sum(e_dims)
    M = np.zeros((rows, cols), dtype=int)
    if rows == 0 or cols == 0:
        return M
    v_off, e_off = {}, {}
    r = 0
    for v, d in zip(VERTICES, v_dims):
        v_off[v] = r
        r += d
    c = 0
    for e, d in zip(EDGES, e_dims):
        e_off[e] = c
        c += d
    for e, de in zip(EDGES, e_dims):
        if de == 0:
            continue
        for v, dv in zip(VERTICES, v_dims):
            if dv == 0:
                continue
            sign = incidence(v, e)
            if sign == 0:
                continue
            ext = extension_matrix(v, e, level)
            M[v_off[v]:v_off[v]+dv, e_off[e]:e_off[e]+de] = sign * ext
    return M


def filtration_matrix(sigma: str, i: int) -> np.ndarray:
    src = components_of(STAR[sigma] & X[i])
    tgt = components_of(STAR[sigma] & X[i + 1])
    M = np.zeros((len(tgt), len(src)), dtype=int)
    for k, comp_src in enumerate(src):
        for j, comp_tgt in enumerate(tgt):
            if comp_src.issubset(comp_tgt):
                M[j, k] = 1
                break
        else:
            raise ValueError("filtration map ill-defined")
    return M


def verify_monomorphism() -> bool:
    for i in range(3):
        for sigma in SIMPLICES:
            M = filtration_matrix(sigma, i)
            if M.size == 0:
                continue
            if np.linalg.matrix_rank(M) < M.shape[1]:
                return False
    return True


def verify_naturality() -> bool:
    for e in EDGES:
        for v in VERTICES:
            if incidence(v, e) == 0:
                continue
            for i in range(3):
                ext_i = extension_matrix(v, e, i)
                ext_i1 = extension_matrix(v, e, i + 1)
                psi_e = filtration_matrix(e, i)
                psi_v = filtration_matrix(v, i)
                if not np.array_equal(psi_v @ ext_i, ext_i1 @ psi_e):
                    return False
    return True


def homology_dims(level: int) -> tuple[int, int]:
    M = boundary_matrix(level)
    rows, cols = M.shape
    rk = int(np.linalg.matrix_rank(M)) if cols else 0
    ker = cols - rk
    return rows - rk, ker


def chain_map(i: int, q: int) -> np.ndarray:
    items = VERTICES if q == 0 else EDGES
    blocks = [filtration_matrix(s, i) for s in items]
    rows = sum(b.shape[0] for b in blocks)
    cols = sum(b.shape[1] for b in blocks)
    M = np.zeros((rows, cols), dtype=int)
    r = c = 0
    for b in blocks:
        if b.size > 0:
            M[r:r+b.shape[0], c:c+b.shape[1]] = b
        r += b.shape[0]
        c += b.shape[1]
    return M


def persistent_rank(i: int, j: int, q: int) -> int:
    # all examples here are 1-dimensional, so we can compute ranks using linear algebra on cycles
    Mi = boundary_matrix(i)
    if q == 0:
        Zi = np.eye(Mi.shape[0], dtype=float)
    else:
        # cycles in degree 1 are kernel of boundary_matrix(i)
        if Mi.shape[1] == 0:
            Zi = np.zeros((Mi.shape[1], 0), dtype=float)
        else:
            # use SVD-based nullspace over R; ranks here are all exact for these tiny integer matrices
            u, s, vh = np.linalg.svd(Mi.astype(float), full_matrices=True)
            rank = int((s > 1e-10).sum())
            Zi = vh[rank:].T.copy()
            if Zi.size == 0:
                Zi = np.zeros((Mi.shape[1], 0), dtype=float)
    cur = Zi
    for k in range(i, j):
        cm = chain_map(k, q).astype(float)
        if cur.shape[1] == 0:
            cur = np.zeros((cm.shape[0], 0), dtype=float)
        else:
            cur = cm @ cur
    if q == 0:
        Bj = boundary_matrix(j).astype(float)
    else:
        Bj = np.zeros((boundary_matrix(j).shape[1], 0), dtype=float)
    if cur.shape[1] == 0:
        return 0
    if Bj.size == 0:
        return int(np.linalg.matrix_rank(cur))
    aug = np.hstack([cur, Bj])
    return int(np.linalg.matrix_rank(aug)) - int(np.linalg.matrix_rank(Bj))


def barcode(q: int, n_levels: int = 4) -> list[tuple[int, float]]:
    def rk(i: int, j: int) -> int:
        if i < 0 or j >= n_levels or i > j:
            return 0
        return persistent_rank(i, j, q)
    bars: list[tuple[int, float]] = []
    for b in range(n_levels):
        for d in range(b + 1, n_levels):
            count = rk(b, d - 1) - rk(b, d) - rk(b - 1, d - 1) + rk(b - 1, d)
            for _ in range(count):
                bars.append((b, d))
        infn = rk(b, n_levels - 1) - rk(b - 1, n_levels - 1)
        for _ in range(infn):
            bars.append((b, float('inf')))
    return bars


def format_barcode(bars: list[tuple[int, float]]) -> list[str]:
    out = []
    for b, d in sorted(bars):
        d_str = '∞' if d == float('inf') else str(int(d))
        out.append(f'[{b},{d_str})')
    return out


def report_dict() -> dict:
    return {
        'costalk_dimensions': {s: [costalk_dim(s, i) for i in range(4)] for s in SIMPLICES},
        'monomorphism_check': verify_monomorphism(),
        'naturality_check': verify_naturality(),
        'homology_dimensions': {f'level_{i}': {'H0': homology_dims(i)[0], 'H1': homology_dims(i)[1]} for i in range(4)},
        'barcode_H0': format_barcode(barcode(0, 4)),
        'barcode_H1': format_barcode(barcode(1, 4)),
    }


def write_outputs(out_dir: Path | str) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rep = report_dict()
    (out_dir / 'leray_c4_report.json').write_text(json.dumps(rep, indent=2, ensure_ascii=False), encoding='utf-8')
    md = []
    md.append('# Leray C4 combinatorial example\n')
    md.append('This script realises the discussion in Section 2.4 by constructing a zeroth Leray cosheaf filtration on the cycle graph $C_4$ from a refined triangulation of $S^1$.')
    md.append('')
    md.append('## Costalk dimensions')
    for s in SIMPLICES:
        dims = ', '.join(f'i={i}: {costalk_dim(s, i)}' for i in range(4))
        md.append(f'- `{s}`: {dims}')
    md.append('')
    md.append(f"- Monomorphism check: {'PASS' if rep['monomorphism_check'] else 'FAIL'}")
    md.append(f"- Naturality check: {'PASS' if rep['naturality_check'] else 'FAIL'}")
    md.append('')
    md.append('## Homology dimensions by level')
    for lvl, vals in rep['homology_dimensions'].items():
        md.append(f"- {lvl}: H0 = F^{vals['H0']}, H1 = F^{vals['H1']}")
    md.append('')
    md.append('## Barcodes')
    md.append(f"- H0: {rep['barcode_H0']}")
    md.append(f"- H1: {rep['barcode_H1']}")
    md.append('')
    md.append('These agree with the zero-extended constant filtration pattern of Example 3.5: H0 = [0,∞), [0,2) and H1 = [3,∞).')
    (out_dir / 'leray_c4_report.md').write_text('\n'.join(md), encoding='utf-8')
    return rep


def main() -> None:
    rep = write_outputs(Path(__file__).resolve().parent / 'out')
    print(json.dumps(rep, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
