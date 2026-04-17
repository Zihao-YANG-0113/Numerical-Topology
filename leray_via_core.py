from __future__ import annotations

from typing import Tuple

import numpy as np
import sympy as sp

from cosheaf_core import FilteredCosheaf
from leray_c4_combinatorial import (
    EDGES,
    SIMPLICES,
    costalk_dim,
    extension_matrix,
    filtration_matrix,
)


def _simplex_tuple(name: str) -> Tuple[str, ...]:
    if len(name) == 1:
        return (name,)
    return tuple(sorted(name))


def _np_to_sympy(M: np.ndarray) -> sp.Matrix:
    rows, cols = M.shape
    if rows == 0 or cols == 0:
        return sp.zeros(rows, cols)
    return sp.Matrix(M.tolist())


def build_leray_c4_cosheaf() -> FilteredCosheaf:
    """Assemble the Leray C4 data into a FilteredCosheaf object that can be
    consumed by cosheaf_core.run_coscythe."""
    levels = 4
    simplex_tuples = [_simplex_tuple(s) for s in SIMPLICES]
    name_of = {_simplex_tuple(s): s for s in SIMPLICES}

    costalk_dims = []
    extension_maps = []
    for r in range(levels):
        dims = {st: costalk_dim(name_of[st], r) for st in simplex_tuples}
        ext = {}
        for edge_name in EDGES:
            tau = _simplex_tuple(edge_name)
            for vertex_char in edge_name:
                sigma = (vertex_char,)
                ext[(tau, sigma)] = _np_to_sympy(
                    extension_matrix(vertex_char, edge_name, r)
                )
        costalk_dims.append(dims)
        extension_maps.append(ext)

    level_maps = []
    for i in range(levels - 1):
        mp = {}
        for st in simplex_tuples:
            mp[st] = _np_to_sympy(filtration_matrix(name_of[st], i))
        level_maps.append(mp)

    return FilteredCosheaf(
        simplices=simplex_tuples,
        levels=levels,
        costalk_dims=costalk_dims,
        extension_maps=extension_maps,
        level_maps=level_maps,
        name="leray_c4_via_core",
    )
