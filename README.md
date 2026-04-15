# Cosheaf experiments: complete package

This folder contains a runnable standalone implementation for the examples, checks, and algorithmic supplements discussed in the essay.

## Run everything

```bash
cd cosheaf_complete_suite
PYTHONPATH=. python scripts/run_all.py
PYTHONPATH=. pytest -q
```

You can also run the parts individually:

```bash
PYTHONPATH=. python scripts/run_c4_and_gallery.py
PYTHONPATH=. python scripts/run_zero_extended_check.py
PYTHONPATH=. python scripts/run_complexity_experiment.py
PYTHONPATH=. python scripts/run_leray_c4_combinatorial.py
```

## Included parts

### Core implementation
- `cosheaf_core.py`: filtered cosheaves, chain complexes, persistent-rank computation, and the `CoScythe` / `CoReducePair` routines.

### Main worked examples and verification
- `scripts/run_c4_and_gallery.py`: runs the paper's non-constant `C4` example together with three additional compatibility-gallery examples.
- `scripts/run_zero_extended_check.py`: verifies that original and Morse barcodes agree on zero-extended constant filtrations.
- `scripts/run_leray_c4_combinatorial.py`: runs a combinatorial Leray-cosheaf example on `C4`, built from a refined triangulation of `S^1`, and writes a report.
- `leray_c4_combinatorial.py`: standalone module implementing the Leray example, including costalk dimensions, naturality checks, homology dimensions, and barcodes.

### Complexity experiment
- `scripts/run_complexity_experiment.py`: simple runtime experiments with output CSV and plots.

### Tests
- `tests/test_examples.py`: regression tests for the original `C4` example and the zero-extended checks.
- `tests/test_leray_example.py`: regression test for the Leray `C4` example.

### Outputs generated / included
- `out/gallery_report.md`, `out/gallery_report.json`: compatibility gallery summary.
- `out/c4_nonconstant_trace.json`: step-by-step trace of the main `C4` run.
- `out/zero_extended_check.json`: agreement check for zero-extended constant filtrations.
- `out/complexity_experiment.csv`, `out/complexity_vs_N.png`, `out/complexity_vs_d.png`: runtime experiment outputs.
- `out/leray_c4_report.json`, `out/leray_c4_report.md`: Leray `C4` report.

## Dependencies

See `requirements.txt`.
