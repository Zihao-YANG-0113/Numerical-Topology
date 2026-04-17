"""Leray C4 example routed through the full FilteredCosheaf / CoScythe
pipeline, verifying Theorem 4.10 on Leray-constructed data.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cosheaf_core import (
    FilteredCosheaf,
    bars_equal,
    format_barcodes,
    morse_persistent_barcodes,
    original_persistent_barcodes,
    run_coscythe,
    verify_matching_acyclic,
)
from leray_via_core import build_leray_c4_cosheaf


def _format_simplex(s) -> str:
    return "".join(s)


def _build_report(fc: FilteredCosheaf, result, bars_orig, bars_morse, agree: bool) -> dict:
    return {
        "name": fc.name,
        "levels": fc.levels,
        "simplices": [_format_simplex(s) for s in fc.simplices],
        "costalk_dimensions": {
            _format_simplex(s): [fc.costalk_dims[r][s] for r in range(fc.levels)]
            for s in fc.simplices
        },
        "matching": [[_format_simplex(a), _format_simplex(b)] for a, b in result.matching],
        "critical": [_format_simplex(s) for s in result.critical],
        "critical_count": len(result.critical),
        "reduction_ratio": len(result.critical) / len(fc.simplices),
        "matching_is_acyclic": verify_matching_acyclic(result.matching),
        "barcodes_original": format_barcodes(bars_orig),
        "barcodes_morse": format_barcodes(bars_morse),
        "theorem_4_10_verified": agree,
        "timings": result.timings,
    }


def _write_markdown(report: dict, path: Path) -> None:
    md = []
    md.append("# Leray C4 example through the FilteredCosheaf / CoScythe pipeline\n")
    md.append(
        "This report routes the Leray-constructed C4 cosheaf "
        "(see `leray_c4_combinatorial.py` and Section 2.4 of the essay) "
        "through the same pipeline as the main Example 5.3, verifying "
        "Theorem 4.10 on Leray data.\n"
    )

    md.append("## Costalk dimensions\n")
    md.append("| simplex | L^0 | L^1 | L^2 | L^3 |")
    md.append("|---------|-----|-----|-----|-----|")
    for s, dims in report["costalk_dimensions"].items():
        md.append(f"| `{s}` | " + " | ".join(str(d) for d in dims) + " |")
    md.append("")

    md.append("## CoScythe output\n")
    md.append(f"- Input simplices: {len(report['simplices'])}")
    md.append(f"- Critical simplices: {report['critical_count']}")
    md.append(f"- Reduction ratio: {report['reduction_ratio']:.3f}")
    md.append(f"- Matching: {report['matching']}")
    md.append(f"- Critical: {report['critical']}")
    md.append(f"- Matching is acyclic (Prop 5.1): {'PASS' if report['matching_is_acyclic'] else 'FAIL'}\n")

    md.append("## Barcodes\n")
    md.append(f"- Original: {report['barcodes_original']}")
    md.append(f"- Morse:    {report['barcodes_morse']}\n")

    status = "PASS" if report["theorem_4_10_verified"] else "FAIL"
    md.append(
        f"## Theorem 4.10 verified: original barcode = Morse barcode -> {status}\n"
    )

    path.write_text("\n".join(md), encoding="utf-8")


def main() -> dict:
    fc = build_leray_c4_cosheaf()
    result = run_coscythe(fc, tie_breaker="minlex")
    bars_orig, _ = original_persistent_barcodes(fc)
    bars_morse, _ = morse_persistent_barcodes(fc, result)
    agree = bars_equal(bars_orig, bars_morse)
    assert agree, (
        f"Theorem 4.10 violated on Leray data: "
        f"{format_barcodes(bars_orig)} != {format_barcodes(bars_morse)}"
    )

    report = _build_report(fc, result, bars_orig, bars_morse, agree)
    out_dir = ROOT / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "leray_via_core_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    _write_markdown(report, out_dir / "leray_via_core_report.md")
    return report


if __name__ == "__main__":
    rep = main()
    out = ROOT / "out"
    print("Wrote:", out / "leray_via_core_report.json")
    print("Wrote:", out / "leray_via_core_report.md")
    print("Theorem 4.10 verified:", rep["theorem_4_10_verified"])
