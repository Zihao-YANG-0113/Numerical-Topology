from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from leray_c4_combinatorial import write_outputs

if __name__ == '__main__':
    out = ROOT / 'out'
    rep = write_outputs(out)
    print('Wrote:', out / 'leray_c4_report.json')
    print('Wrote:', out / 'leray_c4_report.md')
    print(rep)
