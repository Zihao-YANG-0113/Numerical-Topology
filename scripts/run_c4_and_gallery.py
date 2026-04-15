from pathlib import Path
import json

from cosheaf_core import (
    build_c4_nonconstant_example,
    build_path_large_reduction_example,
    build_filled_triangle_no_reduction_example,
    build_c4_zero_extended_gallery_example,
    original_persistent_barcodes,
    morse_persistent_barcodes,
    run_coscythe,
    summary_dict,
    write_json,
)

OUT = Path(__file__).resolve().parents[1] / "out"
OUT.mkdir(parents=True, exist_ok=True)


def run_example(builder, trace=False):
    fc = builder()
    result = run_coscythe(fc, record_trace=trace)
    bars_orig, _ = original_persistent_barcodes(fc)
    bars_morse, _ = morse_persistent_barcodes(fc, result)
    return fc, result, bars_orig, bars_morse


def main():
    report = {"examples": []}
    for builder, name, do_trace, tie_breaker in [
        (build_c4_nonconstant_example, "c4_nonconstant", True, "reverselex"),
        (build_path_large_reduction_example, "path_large_reduction", False, "minlex"),
        (build_c4_zero_extended_gallery_example, "c4_gallery_moderate", False, "minlex"),
        (build_filled_triangle_no_reduction_example, "filled_triangle_no_reduction", False, "minlex"),
    ]:
        fc = builder()
        result = run_coscythe(fc, tie_breaker=tie_breaker, record_trace=do_trace)
        bo, _ = original_persistent_barcodes(fc)
        bm, _ = morse_persistent_barcodes(fc, result)
        entry = summary_dict(fc, result, bo, bm)
        report["examples"].append(entry)
        if do_trace:
            write_json(str(OUT / f"{name}_trace.json"), result.trace)
    write_json(str(OUT / "gallery_report.json"), report)
    md = []
    md.append("# Compatibility gallery and worked examples\n")
    for entry in report["examples"]:
        md.append(f"## {entry['name']}")
        md.append(f"- Critical simplices: {entry['critical_count']} / {len(entry['simplices'])}")
        md.append(f"- Reduction ratio: {entry['reduction_ratio']:.3f}")
        md.append(f"- Matching: {entry['matching']}")
        md.append(f"- Critical: {entry['critical']}")
        md.append(f"- Original barcodes: {entry['barcodes_original']}")
        md.append(f"- Morse barcodes: {entry['barcodes_morse']}\n")
    (OUT / "gallery_report.md").write_text("\n".join(md), encoding="utf-8")
    print((OUT / "gallery_report.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
