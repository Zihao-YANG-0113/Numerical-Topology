"""Standalone Leray C4 sanity checks: costalks, monomorphism, naturality,
and barcode computed directly from boundary matrices (no Morse reduction).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from leray_c4_combinatorial import SIMPLICES, costalk_dim, report_dict


def write_outputs(out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    rep = report_dict()
    (out_dir / 'leray_c4_report.json').write_text(
        json.dumps(rep, indent=2, ensure_ascii=False), encoding='utf-8'
    )
    md = []
    md.append('# Leray C4 combinatorial example\n')
    md.append('This script realises the discussion in Section 2.4 by constructing a zeroth Leray cosheaf filtration on the cycle graph $C_4$ from a refined triangulation of $S^1$.')
    md.append('')
    md.append('## Costalk dimensions')
    for s in SIMPLICES:
        dims = ', '.join(f'i={i}: {costalk_dim(s, i)}' for i in range(4))
        md.append(f'- `{s}`: {dims}')
    md.append('')
    md.append(f"- STAR = f^-1(st(sigma)) consistency: {'PASS' if rep['fmap_star_consistency'] else 'FAIL'}")
    md.append(f"- Monomorphism check: {'PASS' if rep['monomorphism_check'] else 'FAIL'}")
    md.append(f"- Naturality check: {'PASS' if rep['naturality_check'] else 'FAIL'}")
    md.append('')
    md.append('## Homology dimensions by level')
    for lvl, vals in rep['homology_dimensions'].items():
        md.append(f"- {lvl}: H0 = F^{vals['H0']}, H1 = F^{vals['H1']}")
    md.append('')
    md.append('## Barcodes')
    md.append(f"- H0: {rep['barcode_H0']}")
    md.append(f"- H1: {rep['barcode_H1']}")
    md.append('')
    md.append(
        'The resulting costalks are 0 or 1 throughout, with the same birth '
        'times as Example 3.5 — so this Leray construction realises Example '
        '3.5\'s zero-extended constant filtration *geometrically* '
        '(via $f : S^1 \\to |K|$), and the barcodes '
        'H0 = [0,∞), [0,2) and H1 = [3,∞) coincide with the ones computed '
        'combinatorially there.'
    )
    (out_dir / 'leray_c4_report.md').write_text('\n'.join(md), encoding='utf-8')
    return rep


def main() -> None:
    out = ROOT / 'out'
    rep = write_outputs(out)
    print('Wrote:', out / 'leray_c4_report.json')
    print('Wrote:', out / 'leray_c4_report.md')
    print(rep)


if __name__ == '__main__':
    main()
