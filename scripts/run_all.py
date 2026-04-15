from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [
    'run_c4_and_gallery.py',
    'run_zero_extended_check.py',
    'run_complexity_experiment.py',
    'run_leray_c4_combinatorial.py',
]

for script in SCRIPTS:
    print(f'=== running {script} ===')
    subprocess.run([sys.executable, str(ROOT / 'scripts' / script)], check=True)
