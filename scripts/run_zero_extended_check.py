from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cosheaf_core import (
    build_path_large_reduction_example,
    build_filled_triangle_no_reduction_example,
    build_c4_zero_extended_gallery_example,
    original_persistent_barcodes,
    morse_persistent_barcodes,
    run_coscythe,
    bars_equal,
    write_json,
    format_barcodes,
)

OUT = Path(__file__).resolve().parents[1] / "out"
OUT.mkdir(parents=True, exist_ok=True)


def check(builder):
    fc = builder()
    result = run_coscythe(fc)
    bo, _ = original_persistent_barcodes(fc)
    bm, _ = morse_persistent_barcodes(fc, result)
    return {
        "name": fc.name,
        "barcodes_original": format_barcodes(bo),
        "barcodes_morse": format_barcodes(bm),
        "agree": bars_equal(bo, bm),
        "critical_count": len(result.critical),
        "total_simplices": len(fc.simplices),
    }


def main():
    data = {
        "checks": [
            check(build_path_large_reduction_example),
            check(build_c4_zero_extended_gallery_example),
            check(build_filled_triangle_no_reduction_example),
        ]
    }
    write_json(str(OUT / "zero_extended_check.json"), data)
    for c in data['checks']:
        print(c)


if __name__ == "__main__":
    main()
