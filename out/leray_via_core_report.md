# Leray C4 example through the FilteredCosheaf / CoScythe pipeline

This report routes the Leray-constructed C4 cosheaf (see `leray_c4_combinatorial.py` and Section 2.4 of the essay) through the same pipeline as the main Example 5.3, verifying Theorem 4.10 on Leray data.

## Costalk dimensions

| simplex | L^0 | L^1 | L^2 | L^3 |
|---------|-----|-----|-----|-----|
| `a` | 1 | 1 | 1 | 1 |
| `b` | 0 | 1 | 1 | 1 |
| `c` | 0 | 1 | 1 | 1 |
| `d` | 1 | 1 | 1 | 1 |
| `ab` | 0 | 0 | 1 | 1 |
| `bc` | 0 | 1 | 1 | 1 |
| `cd` | 0 | 1 | 1 | 1 |
| `ad` | 0 | 0 | 0 | 1 |

## CoScythe output

- Input simplices: 8
- Critical simplices: 6
- Reduction ratio: 0.750
- Matching: [['c', 'bc']]
- Critical: ['a', 'ab', 'ad', 'b', 'cd', 'd']
- Matching is acyclic (Prop 5.1): PASS

## Barcodes

- Original: {'H0': ['[0,2)', '[0,∞)'], 'H1': ['[3,∞)']}
- Morse:    {'H0': ['[0,2)', '[0,∞)'], 'H1': ['[3,∞)']}

## Theorem 4.10 verified: original barcode = Morse barcode -> PASS
