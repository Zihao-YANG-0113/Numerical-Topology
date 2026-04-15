from cosheaf_core import (
    build_c4_nonconstant_example,
    build_path_large_reduction_example,
    build_filled_triangle_no_reduction_example,
    build_c4_zero_extended_gallery_example,
    original_persistent_barcodes,
    morse_persistent_barcodes,
    run_coscythe,
    bars_equal,
    format_barcodes,
)


def test_c4_expected_critical_simplices():
    fc = build_c4_nonconstant_example()
    res = run_coscythe(fc, tie_breaker='reverselex')
    assert set(res.critical) == {('a',), ('d',), ('a','b'), ('a','d')}


def test_c4_barcodes_agree_and_match_paper():
    fc = build_c4_nonconstant_example()
    res = run_coscythe(fc, tie_breaker='reverselex')
    bo, _ = original_persistent_barcodes(fc)
    bm, _ = morse_persistent_barcodes(fc, res)
    assert bars_equal(bo, bm)
    fmt = format_barcodes(bo)
    assert fmt['H0'] == ['[0,2)', '[0,∞)'] or fmt['H0'] == ['[0,∞)', '[0,2)']
    assert fmt['H1'] == ['[3,∞)']


def test_zero_extended_gallery_agrees():
    for builder in [build_path_large_reduction_example, build_c4_zero_extended_gallery_example, build_filled_triangle_no_reduction_example]:
        fc = builder()
        res = run_coscythe(fc)
        bo, _ = original_persistent_barcodes(fc)
        bm, _ = morse_persistent_barcodes(fc, res)
        assert bars_equal(bo, bm)
