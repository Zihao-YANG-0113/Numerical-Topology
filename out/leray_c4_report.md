# Leray C4 combinatorial example

This script realises the discussion in Section 2.4 by constructing a zeroth Leray cosheaf filtration on the cycle graph $C_4$ from a refined triangulation of $S^1$.

## Costalk dimensions
- `a`: i=0: 1, i=1: 1, i=2: 1, i=3: 1
- `b`: i=0: 0, i=1: 1, i=2: 1, i=3: 1
- `c`: i=0: 0, i=1: 1, i=2: 1, i=3: 1
- `d`: i=0: 1, i=1: 1, i=2: 1, i=3: 1
- `ab`: i=0: 0, i=1: 0, i=2: 1, i=3: 1
- `bc`: i=0: 0, i=1: 1, i=2: 1, i=3: 1
- `cd`: i=0: 0, i=1: 1, i=2: 1, i=3: 1
- `ad`: i=0: 0, i=1: 0, i=2: 0, i=3: 1

- Monomorphism check: PASS
- Naturality check: PASS

## Homology dimensions by level
- level_0: H0 = F^2, H1 = F^0
- level_1: H0 = F^2, H1 = F^0
- level_2: H0 = F^1, H1 = F^0
- level_3: H0 = F^1, H1 = F^1

## Barcodes
- H0: ['[0,2)', '[0,∞)']
- H1: ['[3,∞)']

These agree with the zero-extended constant filtration pattern of Example 3.5: H0 = [0,∞), [0,2) and H1 = [3,∞).