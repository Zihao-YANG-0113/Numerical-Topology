from pathlib import Path
import csv
import time
import matplotlib.pyplot as plt

from cosheaf_core import build_random_path_example, run_coscythe

OUT = Path(__file__).resolve().parents[1] / "out"
OUT.mkdir(parents=True, exist_ok=True)


def time_run(num_vertices: int, costalk_dim: int, levels: int = 3, trials: int = 3):
    times = []
    for t in range(trials):
        fc = build_random_path_example(num_vertices=num_vertices, levels=levels, costalk_dim=costalk_dim, seed=t)
        start = time.perf_counter()
        run_coscythe(fc)
        times.append(time.perf_counter() - start)
    return sum(times) / len(times)


def main():
    rows = []
    for nv in [10, 20, 30, 40, 50]:
        avg = time_run(nv, costalk_dim=1, trials=2)
        rows.append({"mode": "vary_N", "num_vertices": nv, "costalk_dim": 1, "seconds": avg})
    for d in [1, 2, 3, 4]:
        avg = time_run(30, costalk_dim=d, trials=2)
        rows.append({"mode": "vary_d", "num_vertices": 30, "costalk_dim": d, "seconds": avg})

    with open(OUT / "complexity_experiment.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["mode", "num_vertices", "costalk_dim", "seconds"])
        writer.writeheader()
        writer.writerows(rows)

    xs = [r['num_vertices'] for r in rows if r['mode'] == 'vary_N']
    ys = [r['seconds'] for r in rows if r['mode'] == 'vary_N']
    plt.figure()
    plt.plot(xs, ys, marker='o')
    plt.xlabel('number of vertices in path complex')
    plt.ylabel('runtime (seconds)')
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(OUT / 'complexity_vs_N.png', dpi=160)
    plt.close()

    xs = [r['costalk_dim'] for r in rows if r['mode'] == 'vary_d']
    ys = [r['seconds'] for r in rows if r['mode'] == 'vary_d']
    plt.figure()
    plt.plot(xs, ys, marker='o')
    plt.xlabel('costalk dimension d')
    plt.ylabel('runtime (seconds)')
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(OUT / 'complexity_vs_d.png', dpi=160)
    plt.close()

    print('Wrote complexity_experiment.csv, complexity_vs_N.png, complexity_vs_d.png')


if __name__ == '__main__':
    main()
