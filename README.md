# Cosheaf experiments: supplementary code for

### "Discrete Morse Theory for Persistent Homology of Filtered Cosheaves"

This repository contains a runnable, self-contained implementation of
the constructions and algorithms in the essay, organised around three
separable components:

1. **The main worked example** (Example 5.3, essay § 5.3):
   a filtered cosheaf on the cycle graph C4 with non-trivial costalks
   and extension maps at level 3.

2. **A Leray cosheaf example** (essay § 2.4):
   an alternative filtered cosheaf on C4 arising from the standard
   Leray construction applied to a refined triangulation of S^1.

3. **Complexity experiment** (essay § 5.2, Theorem 5.3):
   empirical scaling of CoScythe against input size N and maximum
   costalk dimension d.

All three components are built on the same core library
`cosheaf_core.py`, which provides `FilteredCosheaf`, the `CoScythe`
and `CoReducePair` routines, and the computation of persistent
barcodes by direct linear algebra. The two worked examples are
independently checked to satisfy Theorem 4.10 (Morse reduction
preserves the persistent barcode).

## Quick start

```bash
cd cosheaf_complete_suite
pip install -r requirements.txt
PYTHONPATH=. python scripts/run_all.py
PYTHONPATH=. pytest -q
```

All outputs are written under `out/`. Delete `out/` before rerunning
to fully reproduce the included results from scratch.

## Repository layout

```
cosheaf_complete_suite/
├── cosheaf_core.py                  core library
├── leray_c4_combinatorial.py        Leray geometric data + sanity checks
├── leray_via_core.py                bridge: Leray data → FilteredCosheaf
│
├── scripts/
│   ├── run_c4_and_gallery.py        main C4 example + compatibility gallery
│   ├── run_leray_c4_combinatorial.py   standalone Leray sanity checks
│   ├── run_leray_via_core.py        Leray example through the full pipeline
│   ├── run_zero_extended_check.py   constant-coefficient sanity check
│   ├── run_complexity_experiment.py CoScythe runtime vs N, d
│   ├── plot_figure1.py              reproduces essay Figure 1(a) barcode plot
│   └── run_all.py                   runs everything in order
│
├── tests/                           pytest regression tests
├── out/                             reproducible outputs (JSON / MD / CSV / PNG)
└── requirements.txt
```

## Component 1 — Main C4 example (essay § 5.3)

The main example is specifically designed to exercise the algorithm
**beyond the zero-extended constant case** of Proposition 3.4. At level 3,
the costalks at `b` and `bc` are F^2 (rather than 0 or F), and the
extension maps `C^3(bc ≥ c) = π_1` and `C^3(ab ≥ b) = ι` are non-identity
maps between costalks of different dimension. This takes the input
outside the reach of the classical Mischaikow–Nanda algorithm [3],
which operates on filtrations of subcomplexes with constant coefficients,
and makes the example a genuine test of the cosheaf-theoretic
generalisation developed in the essay.

The filtered cosheaf of Example 5.3 is assembled directly from the
cosheaf data given in the essay.

`scripts/run_c4_and_gallery.py` performs the following steps:

1. Build the filtered cosheaf object
   (`build_c4_nonconstant_example` in `cosheaf_core.py`).
2. Compute the original persistent barcode by direct computation of
   boundary matrices, ranks, and quotient maps
   (`original_persistent_barcodes`).
3. Run CoScythe to obtain a compatible acyclic matching Σ and the
   Morse boundary operators (`run_coscythe`).
4. Compute the Morse persistent barcode from the reduced complex
   (`morse_persistent_barcodes`).
5. Assert that the two barcodes coincide (`bars_equal`), verifying
   Theorem 4.10 on this example.

The script additionally runs three smaller examples spanning the full
range of reduction strength — reduction ratios ≈ 14%, 75%, 100% — each
independently verified to satisfy `original barcode = Morse barcode`.
Together they illustrate the limitation flagged in the essay's
Discussion: Definition 3.3's compatibility requirement can restrict
the algorithm's efficiency, with `filled_triangle_no_reduction` being
the extreme case (100%) where no pair is matchable.

Expected output for the main example: H0 = {[0, ∞), [0, 2)}   H1 = {[3, ∞)} with four critical simplices `a, d, ab, ad`.

## Component 2 — Leray C4 example (essay § 2.4)

We realise a filtered cosheaf on C4 from the zeroth Leray
construction applied to a continuous map `f : S^1 → |K|`.

### Construction

We model S^1 as the cyclic quotient of 16 points labelled
`v_0, v_1, …, v_15`, with the K-vertices placed at

  `a = v_0,   b = v_4,   c = v_8,   d = v_12`,

and the K-edges as the three consecutive arcs between adjacent
vertices:

  `ab = {v_1, v_2, v_3}`,     `bc = {v_5, v_6, v_7}`,
  `cd = {v_9, v_10, v_11}`,   `ad = {v_13, v_14, v_15}`.

The open star of each K-simplex is the union of the simplex itself
and its incident higher simplices. The filtration of subspaces is

  `X_0 = {v_0, v_12}`                            (just `a` and `d`)
  `X_1 = X_0 ∪ {v_4, …, v_11}`                   (adds `b`, `c`, arcs `bc`, `cd`)
  `X_2 = X_1 ∪ {v_1, v_2, v_3}`                  (adds arc `ab`)
  `X_3 = X_2 ∪ {v_13, v_14, v_15} = S^1`         (adds arc `ad`)

The zeroth Leray cosheaf assigns to each simplex σ the vector space

  `L^i(σ) = H_0( X_i ∩ st(σ) ; F )`,

whose dimension equals the number of connected components of
`X_i ∩ st(σ)` under the cyclic adjacency on the 16 points. Costalk
dimensions at each filtration level are:

| level | `a` | `b` | `c` | `d` | `ab` | `bc` | `cd` | `ad` |
|-------|-----|-----|-----|-----|------|------|------|------|
| `L^0` |  1  |  0  |  0  |  1  |  0   |  0   |  0   |  0   |
| `L^1` |  1  |  1  |  1  |  1  |  0   |  1   |  1   |  0   |
| `L^2` |  1  |  1  |  1  |  1  |  1   |  1   |  1   |  0   |
| `L^3` |  1  |  1  |  1  |  1  |  1   |  1   |  1   |  1   |

### Worked derivations of table entries

Unrolling the formula `L^i(σ) = H_0(X_i ∩ st(σ); F)` on four representative entries:

- **`L^0(a)`**: `st(a) = {v_0, v_1, v_2, v_3, v_13, v_14, v_15}`, and
  `X_0 = {v_0, v_12}`. The intersection is `{v_0}` — a single point,
  hence a single connected component. So `dim L^0(a) = 1`.

- **`L^1(b)`**: `st(b) = {v_1, v_2, v_3, v_4, v_5, v_6, v_7}`, and
  `X_1 = {v_0, v_4, v_5, …, v_11, v_12}`. The intersection is
  `{v_4, v_5, v_6, v_7}`, four *cyclically consecutive* points forming
  a single arc. So `dim L^1(b) = 1`.

- **`L^1(bc)`**: `st(bc) = {v_5, v_6, v_7}` (an edge's open star is
  the edge itself), so `st(bc) ∩ X_1 = {v_5, v_6, v_7}`, again a
  single arc, giving `dim L^1(bc) = 1`.

The cyclic adjacency is the reason `{v_0, v_1, v_2, v_3, v_13, v_14, v_15}`
counts as *one* component and not two. In the code, this is handled by
`num_components` in [leray_c4_combinatorial.py](leray_c4_combinatorial.py),
which measures gaps modulo 16 and treats the special case `n = 16`
(the whole circle) as a single component. Every entry of the table
above is the exact output of `costalk_dim(σ, i) = num_components(STAR[σ] & X[i])`
on the `F_MAP`-derived `STAR` and `X` sets; the test
`test_leray_c4_barcode_and_checks` re-runs the full recomputation on
every CI pass.

Extension maps `L^i(τ ≥ σ)` are induced by the inclusion
`X_i ∩ st(τ) ⊆ X_i ∩ st(σ)` at the H_0 level: each component of
`X_i ∩ st(τ)` is sent to the unique component of `X_i ∩ st(σ)`
containing it. Filtration maps `Ψ^i_σ : L^i(σ) → L^{i+1}(σ)` are
defined analogously from `X_i ⊆ X_{i+1}`.

### Scripts

Two scripts operate on this data, at different levels of the framework:

- **`scripts/run_leray_c4_combinatorial.py`** — a lightweight,
  standalone verification:
    - Computes costalk dimensions.
    - Verifies that every filtration map is injective
      (the monomorphism condition of Definition 2.5).
    - Verifies that each commuting square `Ψ^i_σ ∘ L^i(τ ≥ σ) =
      L^{i+1}(τ ≥ σ) ∘ Ψ^i_τ` commutes (naturality of Ψ^i).
    - Computes the barcode directly from boundary matrices and
      persistent ranks, without Morse reduction.

- **`scripts/run_leray_via_core.py`** — routes the Leray data through
  the full pipeline used for the main example:
    - Packages the Leray data as a `FilteredCosheaf` object
      (whose `__post_init__` re-checks monomorphism and naturality).
    - Runs `run_coscythe` to obtain a compatible acyclic matching Σ.
    - Computes the Morse boundary and the Morse persistent barcode.
    - Asserts that the Morse barcode equals the original barcode,
      verifying Theorem 4.10 on Leray data.

Both scripts agree on the barcode:

  H0 = {[0, ∞), [0, 2)}   H1 = {[3, ∞)}

## Component 3 — Complexity experiment (essay § 5.2)

`scripts/run_complexity_experiment.py` measures the end-to-end runtime
of CoScythe (search, reduction, and barcode extraction) on families of
randomly generated filtered cosheaves supported on path graphs.

**Test complex.** A path graph on `V` vertices (so `V−1` edges and
`N = 2V − 1` simplices in total), 3 filtration levels, every simplex
carrying an `F^d` costalk at every level, identity level maps (so
naturality is automatic), and near-identity invertible extension maps.
Path graphs are chosen so that the coface-count bound `p` in
Theorem 5.3 stays small and constant across `V`, isolating the effect
of `V` (equivalently `N`) and `d` on the runtime.

**Setup.** `V ∈ {10, 20, 30, 40, 50}` (i.e. `N ∈ {19, 39, 59, 79, 99}`)
with `d = 1`; and `d ∈ {1, 2, 3, 4}` with `V = 30`. 2 trials per point,
averaged. In the CSV, `V` appears under the self-explanatory column
name `num_vertices`; the sweep plot filename `complexity_vs_N.png` uses
`N` as a shorthand for the swept variable. Outputs:

- `out/complexity_experiment.csv`: one row per `(num_vertices, d)` pair.
- `out/complexity_vs_N.png`: runtime against `V` with `d` fixed.
- `out/complexity_vs_d.png`: runtime against `d` with `V` fixed.

Over this range, runtime grows monotonically in both `V` and `d` and
stays well below the `O(n N p m̃ d^θ + M_Σ^3)` bound of Theorem 5.3.
The experiment verifies the **direction** of scaling rather than
fitting exponents, since `V ≤ 50` and `d ≤ 4` are too small for tight
asymptotic estimates.

## Tests

`tests/` contains pytest regression tests covering:

- `test_examples.py`: main C4 example critical simplices and barcode
  agreement; zero-extended gallery agreement; CoScythe-produced
  matchings are acyclic (Proposition 5.1).
- `test_leray_example.py`: Leray standalone monomorphism, naturality,
  `STAR = f^{-1}(st(σ))` consistency, and barcode.
- `test_leray_via_core.py`: agreement between original and Morse
  barcodes on Leray data routed through the full pipeline.

Run:

```bash
PYTHONPATH=. pytest -q
```

Expected: 6 tests passed.

## Outputs

Contents of `out/`:

- `gallery_report.md`, `gallery_report.json` — compatibility gallery.
- `c4_nonconstant_trace.json` — step-by-step trace of the main
  CoScythe run on Example 5.3.
- `zero_extended_check.json` — agreement check on zero-extended
  constant filtrations (Proposition 3.4 regime).
- `leray_c4_report.md`, `leray_c4_report.json` — standalone Leray
  sanity report.
- `leray_via_core_report.md`, `leray_via_core_report.json` — Leray
  example through the full pipeline, with Theorem 4.10 verification
  and matching-acyclicity check (Proposition 5.1).
- `figure1_barcodes.png` — reproduction of essay Figure 1(a)
  (persistent barcodes of the filtered Morse complex for both the
  main C4 non-constant example and the Leray C4 example).
- `complexity_experiment.csv`, `complexity_vs_N.png`,
  `complexity_vs_d.png` — complexity experiment outputs.

## Dependencies

See `requirements.txt`: `numpy`, `sympy`, `scipy`, `matplotlib`, `pytest`.
