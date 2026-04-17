"""Reproduces Figure 1(a) of the essay: the persistent barcode of the main
non-constant C4 example, plus a companion plot for the Leray C4 example.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt

from cosheaf_core import (
    build_c4_nonconstant_example,
    morse_persistent_barcodes,
    original_persistent_barcodes,
    run_coscythe,
)
from leray_via_core import build_leray_c4_cosheaf

OUT = ROOT / "out"
OUT.mkdir(parents=True, exist_ok=True)


def _plot_one_example(ax_h0, ax_h1, bars, levels, title):
    inf_x = levels - 0.5
    for ax, k, colour in [(ax_h0, 0, "#1f77b4"), (ax_h1, 1, "#d62728")]:
        blist = sorted(bars.get(k, []), key=lambda bd: (bd[0], float("inf") if bd[1] is None else bd[1]))
        for i, (b, d) in enumerate(blist):
            end = inf_x if d is None else d
            ax.hlines(i, b, end, linewidth=5, color=colour)
            if d is None:
                ax.plot([inf_x], [i], marker=">", color=colour, markersize=8, clip_on=False)
                label = f"[{b}, inf)"
            else:
                label = f"[{b}, {d})"
            ax.text(inf_x + 0.15, i, label, va="center", fontsize=9)
        ax.set_ylabel(f"H_{k}", rotation=0, labelpad=18, fontsize=11)
        ax.set_xlim(-0.3, inf_x + 1.0)
        ax.set_ylim(-0.6, max(1, len(blist)) - 0.4)
        ax.set_yticks([])
        ax.set_xticks(range(levels))
        ax.grid(axis="x", linestyle=":", alpha=0.6)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
    ax_h0.set_title(title, fontsize=11)
    ax_h1.set_xlabel("filtration level i")


def _run(fc):
    result = run_coscythe(fc, tie_breaker="minlex")
    bars_orig, _ = original_persistent_barcodes(fc)
    bars_morse, _ = morse_persistent_barcodes(fc, result)
    return bars_orig, bars_morse, result


def main() -> None:
    fc_main = build_c4_nonconstant_example()
    bars_orig_main, bars_morse_main, _ = _run(fc_main)

    fc_leray = build_leray_c4_cosheaf()
    bars_orig_leray, bars_morse_leray, _ = _run(fc_leray)

    fig, axes = plt.subplots(2, 2, figsize=(10, 5), sharex="col")
    _plot_one_example(axes[0, 0], axes[1, 0], bars_morse_main, fc_main.levels,
                      "Main C4 non-constant example (essay Figure 1(a))")
    _plot_one_example(axes[0, 1], axes[1, 1], bars_morse_leray, fc_leray.levels,
                      "Leray C4 example (essay Section 2.4)")
    fig.suptitle("Persistent barcodes of the filtered Morse complex", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "figure1_barcodes.png", dpi=160)
    plt.close(fig)

    # Sanity: original and Morse barcodes agree, so either set reproduces Fig 1.
    assert bars_orig_main == bars_morse_main
    assert bars_orig_leray == bars_morse_leray
    print("Wrote:", OUT / "figure1_barcodes.png")


if __name__ == "__main__":
    main()
