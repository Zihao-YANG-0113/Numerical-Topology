from leray_c4_combinatorial import report_dict


def test_leray_c4_barcode_and_checks():
    rep = report_dict()
    assert rep['monomorphism_check'] is True
    assert rep['naturality_check'] is True
    assert rep['barcode_H0'] == ['[0,2)', '[0,∞)'] or rep['barcode_H0'] == ['[0,∞)', '[0,2)']
    assert rep['barcode_H1'] == ['[3,∞)']
