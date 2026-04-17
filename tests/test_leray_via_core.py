from cosheaf_core import (
    bars_equal,
    format_barcodes,
    morse_persistent_barcodes,
    original_persistent_barcodes,
    run_coscythe,
    verify_matching_acyclic,
)
from leray_via_core import build_leray_c4_cosheaf


def test_leray_via_core_barcodes_agree_and_match_paper():
    fc = build_leray_c4_cosheaf()
    res = run_coscythe(fc, tie_breaker="minlex")
    bo, _ = original_persistent_barcodes(fc)
    bm, _ = morse_persistent_barcodes(fc, res)
    assert bars_equal(bo, bm)
    assert verify_matching_acyclic(res.matching)
    fmt = format_barcodes(bo)
    assert fmt["H0"] == ["[0,2)", "[0,∞)"] or fmt["H0"] == ["[0,∞)", "[0,2)"]
    assert fmt["H1"] == ["[3,∞)"]
