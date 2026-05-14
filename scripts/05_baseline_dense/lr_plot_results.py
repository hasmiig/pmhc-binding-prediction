"""
lr_plot_results.py
==================
Plots ablation study results for the logistic regression baseline.

Mirrors the box-plot + diamond style of plddt_cv_performance.py.

Rows  = thresholds (50, 60, 70, 80)
Cols  = metrics    (AUC, AP, F1, Precision, Recall, TPR@FPR10)
X-axis per panel = ablation configurations:
    baseline | layer_norm | l2 | l1 | layer_norm_l2 | layer_norm_l1

Box + jitter = 5 CV folds (val rows)
Diamond      = test row for that threshold / config

Usage:
    python lr_plot_results.py \\
        --metrics_csv  /path/to/all_metrics.csv \\
        --output_png   /path/to/lr_ablation.png
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--metrics_csv", required=True,
                    help="Path to all_metrics.csv produced by lr_model.py")
parser.add_argument("--output_png",  required=True,
                    help="Where to save the output PNG")
args = parser.parse_args()

# ── load ───────────────────────────────────────────────────────────────────────
df = pd.read_csv(args.metrics_csv)

# ── build config label ─────────────────────────────────────────────────────────
# Each unique (layer_norm, regularization) combo gets a short human-readable name.

def config_label(row):
    ln  = row["layer_norm"]
    reg = str(row["regularization"]) if pd.notna(row["regularization"]) else "none"
    if not ln and reg == "none":
        return "baseline"
    parts = []
    if ln:
        parts.append("LN")
    if reg != "none":
        parts.append(reg.upper())
    return "+".join(parts)

df["config"] = df.apply(config_label, axis=1)

# Canonical order for the x-axis
CONFIG_ORDER = ["baseline", "LN", "L2", "L1", "LN+L2", "LN+L1"]

# Keep only configs present in the data (in order)
present_configs = [c for c in CONFIG_ORDER if c in df["config"].unique()]

THRESHOLDS  = sorted(df["threshold"].unique())
METRICS     = ["auc", "ap", "f1", "precision", "recall", "tpr_at_fpr10"]
METRIC_LABELS = {
    "auc":          "AUC",
    "ap":           "AP",
    "f1":           "F1",
    "precision":    "Precision",
    "recall":       "Recall",
    "tpr_at_fpr10": "TPR@FPR10",
}

n_configs = len(present_configs)
n_thresholds = len(THRESHOLDS)
n_metrics    = len(METRICS)

# ── colour palette (one colour per config) ─────────────────────────────────────
palette   = plt.rcParams["axes.prop_cycle"].by_key()["color"]
color_map = {cfg: palette[i % len(palette)] for i, cfg in enumerate(present_configs)}

# ── layout constants ───────────────────────────────────────────────────────────
BREAK_GAP = 0.15          # gap between CV and test panels to trigger y-axis break
_TEXT_BOX = dict(boxstyle="round,pad=0.15", facecolor="white",
                 alpha=0.75, edgecolor="none")
rng = np.random.default_rng(42)

col_width = max(4.0, 2.0 * n_configs)
fig = plt.figure(figsize=(col_width * n_metrics, 4.5 * n_thresholds))

outer_gs = GridSpec(
    n_thresholds, n_metrics, figure=fig,
    hspace=0.65, wspace=0.38,
    bottom=0.16, top=0.92,
)

# ── helper functions ───────────────────────────────────────────────────────────

def _style(ax, show_xticks):
    ax.set_xticks(np.arange(n_configs))
    ax.set_xlim(-0.6, n_configs - 0.4)
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(axis="y", lw=0.5, alpha=0.4, linestyle="--")
    ax.spines["right"].set_visible(False)
    if show_xticks:
        ax.set_xticklabels(present_configs, fontsize=9,
                           rotation=30, ha="right")
    else:
        ax.set_xticklabels([])
        ax.tick_params(bottom=False)


def _draw_break_marks(ax_top, ax_bot):
    d = 0.03
    for ax, y0 in [(ax_top, 0), (ax_bot, 1)]:
        kw = dict(transform=ax.transAxes, color="k", clip_on=False, lw=1.2)
        ax.plot((-d, +d), (y0 - 2*d, y0 + 2*d), **kw)
        ax.plot((1-d, 1+d), (y0 - 2*d, y0 + 2*d), **kw)


def _draw_boxes(ax, cv_sub, metric):
    for j, cfg in enumerate(present_configs):
        vals = pd.to_numeric(
            cv_sub[cv_sub["config"] == cfg][metric]
        ).dropna().values
        if len(vals) == 0:
            continue
        color = color_map[cfg]
        ax.boxplot(
            vals, positions=[j], widths=0.55,
            patch_artist=True,
            boxprops=dict(facecolor=color, alpha=0.35),
            medianprops=dict(color="black", lw=2),
            whiskerprops=dict(color=color, lw=1.2),
            capprops=dict(color=color, lw=1.2),
            showfliers=False,
        )
        jitter = rng.uniform(-0.13, 0.13, size=len(vals))
        ax.scatter(j + jitter, vals, color=color, s=30, zorder=5,
                   alpha=0.85, edgecolors="white", linewidths=0.4)
        median  = np.median(vals)
        y_label = np.max(vals) + 0.012
        ax.text(j, y_label, f"{median:.3f}",
                ha="center", va="bottom", fontsize=7.5,
                fontweight="bold", color="black", zorder=7, bbox=_TEXT_BOX)


def _draw_diamonds(ax, test_sub, metric, label_above=False):
    for j, cfg in enumerate(present_configs):
        row = test_sub[test_sub["config"] == cfg]
        if not len(row):
            continue
        val   = pd.to_numeric(row[metric]).values[0]
        color = color_map[cfg]
        ax.scatter([j], [val], marker="D", color=color, s=55, zorder=6,
                   edgecolors="black", linewidths=0.7)
        offset, va = (0.012, "bottom") if label_above else (-0.018, "top")
        ax.text(j, val + offset, f"{val:.3f}",
                ha="center", va=va, fontsize=7, color="black",
                style="italic", zorder=7, bbox=_TEXT_BOX)


# ── main plotting loop ─────────────────────────────────────────────────────────
for row_i, thr in enumerate(THRESHOLDS):
    thr_df   = df[df["threshold"] == thr]
    cv_sub   = thr_df[thr_df["eval_set"] == "val"]
    test_sub = thr_df[thr_df["eval_set"] == "test"]

    for col_i, metric in enumerate(METRICS):
        all_cv  = pd.to_numeric(cv_sub[metric]).dropna().values
        gt_vals = pd.to_numeric(test_sub[metric]).dropna().values \
                  if len(test_sub) else np.array([])

        needs_break = (
            len(gt_vals) > 0 and len(all_cv) > 0 and
            all_cv.min() - gt_vals.max() > BREAK_GAP
        )

        cell = outer_gs[row_i, col_i]

        if needs_break:
            inner  = GridSpecFromSubplotSpec(2, 1, subplot_spec=cell,
                                             hspace=0.1, height_ratios=[3, 1])
            ax_top = fig.add_subplot(inner[0])
            ax_bot = fig.add_subplot(inner[1])

            _draw_boxes(ax_top, cv_sub, metric)
            ax_top.set_ylim(max(0.0, all_cv.min() - 0.04),
                            min(1.0, all_cv.max() + 0.12))
            ax_top.spines["bottom"].set_visible(False)
            ax_top.spines["top"].set_visible(False)
            _style(ax_top, show_xticks=False)

            _draw_diamonds(ax_bot, test_sub, metric, label_above=True)
            ax_bot.set_ylim(max(0.0, gt_vals.min() - 0.06),
                            gt_vals.max() + 0.06)
            ax_bot.spines["top"].set_visible(False)
            ax_bot.spines["right"].set_visible(False)
            _style(ax_bot, show_xticks=True)

            _draw_break_marks(ax_top, ax_bot)

            if row_i == 0:
                ax_top.set_title(METRIC_LABELS[metric],
                                 fontsize=12, fontweight="bold", pad=8)
            if col_i == 0:
                ax_top.set_ylabel(f"thr={thr}", fontsize=11, labelpad=6)

        else:
            ax = fig.add_subplot(cell)
            _draw_boxes(ax, cv_sub, metric)
            if len(gt_vals) > 0:
                _draw_diamonds(ax, test_sub, metric, label_above=False)

            all_shown = (np.concatenate([all_cv, gt_vals])
                         if len(gt_vals) > 0 else all_cv)
            if len(all_shown) > 0:
                ax.set_ylim(max(0.0, all_shown.min() - 0.04),
                            min(1.0, all_shown.max() + 0.13))
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            _style(ax, show_xticks=True)

            if row_i == 0:
                ax.set_title(METRIC_LABELS[metric],
                             fontsize=12, fontweight="bold", pad=8)
            if col_i == 0:
                ax.set_ylabel(f"thr={thr}", fontsize=11, labelpad=6)

# ── legend ─────────────────────────────────────────────────────────────────────
patches = [
    mpatches.Patch(color=color_map[c], alpha=0.7, label=c)
    for c in present_configs
]
diamond_handle = mlines.Line2D(
    [], [], marker="D", color="grey", linestyle="None",
    markersize=7, markeredgecolor="black", label="test (diamond)"
)
fig.legend(
    handles=patches + [diamond_handle],
    loc="lower center",
    ncol=n_configs + 1,
    fontsize=10,
    frameon=False,
    bbox_to_anchor=(0.5, 0.04),
)

fig.suptitle(
    "LR Baseline — 5-fold CV ablation study\n"
    "(box+jitter = val folds, diamond = test)",
    fontsize=13, fontweight="bold", y=0.97,
)

plt.savefig(args.output_png, dpi=150, bbox_inches="tight")
print(f"Saved → {args.output_png}")