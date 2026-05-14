"""
plddt_cv_performance.py
=======================
Threshold baseline — 5-fold cross-validation performance.

Evaluation strategy:
    Dendrogram splits:
    - Per fold: find optimal τ on train, report metrics on val  (CV estimate)
    - Global:   median of per-fold τs, report metrics on test.parquet

Threshold criterion: argmin(sqrt(fpr^2 + (1-tpr)^2))

Features:
    pLDDT (pep)    — pep_mean_plddt    (higher = better)
    pLDDT (anchor) — anchor_mean_plddt (higher = better)
    PAE (pep)      — pae_pep           (lower = better → negated before use)

Results saved to cv_results/results_<split_tag>.csv.
Re-running is idempotent — existing rows for these model names are replaced.

Usage:
    python plddt_cv_performance.py \\
        --splits-root /path/to/output_<tag>/splits \\
        --output-dir  /path/to/output_dir \\
        --feature     both         # plddt | pae | both (default: both)
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from sklearn.metrics import (roc_curve, roc_auc_score, average_precision_score,
                             f1_score, precision_score, recall_score)


# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description='Threshold baseline CV evaluation')
parser.add_argument('--splits-root', required=True,
                    help='Path to splits root directory')
parser.add_argument('--output-dir', required=True,
                    help='Directory where cv_results/ and PNG are written')
parser.add_argument('--feature', choices=['plddt', 'pae', 'both'], default='both',
                    help='Which features to use as threshold baseline')
args = parser.parse_args()

SPLITS_ROOT = args.splits_root
OUTPUT_DIR  = os.path.abspath(args.output_dir)

_split_parent = os.path.basename(
    os.path.normpath(os.path.dirname(os.path.normpath(SPLITS_ROOT)))
)
_csv_tag = _split_parent[len('output_'):] if _split_parent.startswith('output_') else _split_parent
_feat_tag = args.feature

RESULTS_DIR = os.path.join(OUTPUT_DIR, 'cv_results')
RESULTS_CSV = os.path.join(RESULTS_DIR, f'results_{_csv_tag}_{_feat_tag}.csv')
OUT_PNG     = os.path.join(OUTPUT_DIR, f'cv_performance_{_split_parent}_{_feat_tag}.png')
os.makedirs(RESULTS_DIR, exist_ok=True)

print(f'splits root : {SPLITS_ROOT}')
print(f'feature     : {args.feature}')
print(f'results CSV : {RESULTS_CSV}')
print(f'output PNG  : {OUT_PNG}')

N_FOLDS = 5

PLDDT_FEATURES = {
    'pLDDT (pep)':    'pep_mean_plddt',
    'pLDDT (anchor)': 'anchor_mean_plddt',
}
PAE_FEATURES = {
    'PAE (pep)': 'pae_pep',
}

if args.feature == 'plddt':
    FEATURES = PLDDT_FEATURES
elif args.feature == 'pae':
    FEATURES = PAE_FEATURES
else:
    FEATURES = {**PLDDT_FEATURES, **PAE_FEATURES}

# Columns where lower = better — negate so roc_curve sees higher = more likely binder
PAE_COLS = {'pae_pep'}

METRICS = ['F1', 'Precision', 'Recall', 'AUC', 'AP', 'TPR@FPR10']


# ── helpers ────────────────────────────────────────────────────────────────────

def negate_pae(df: pd.DataFrame) -> pd.DataFrame:
    """Negate PAE columns so higher score = more likely binder."""
    df = df.copy()
    for col in PAE_COLS:
        if col in df.columns:
            df[col] = -df[col]
    return df


def load_split(fold, part):
    df = pd.read_parquet(f'{SPLITS_ROOT}/fold_{fold}/{part}.parquet')
    return negate_pae(df)


def load_test():
    df = pd.read_parquet(f'{SPLITS_ROOT}/test/test.parquet')
    return negate_pae(df)


def find_threshold(df, col):
    fpr, tpr, thresholds = roc_curve(df['label'], df[col])
    idx = np.argmin(np.sqrt(fpr**2 + (1 - tpr)**2))
    return float(thresholds[idx])


def compute_metrics(df, col, tau):
    scores = df[col].values
    labels = df['label'].values
    preds  = (scores >= tau).astype(int)
    fpr, tpr, _ = roc_curve(labels, scores)
    mask = fpr <= 0.10
    tpr_at_fpr10 = float(tpr[mask][-1]) if mask.any() else 0.0
    return {
        'threshold':  tau,
        'F1':         f1_score(labels, preds, zero_division=0),
        'Precision':  precision_score(labels, preds, zero_division=0),
        'Recall':     recall_score(labels, preds, zero_division=0),
        'AUC':        roc_auc_score(labels, scores),
        'AP':         average_precision_score(labels, scores),
        'TPR@FPR10':  tpr_at_fpr10,
    }


# ── evaluation ─────────────────────────────────────────────────────────────────
rows = []

print('\n── CV (train→τ, val→metrics) ──')
for fold in range(N_FOLDS):
    train = load_split(fold, 'train')
    val   = load_split(fold, 'val')
    for model_name, col in FEATURES.items():
        if col not in train.columns:
            print(f'  WARNING: {col} not in fold {fold} train — skipping')
            continue
        tau  = find_threshold(train, col)
        mets = compute_metrics(val, col, tau)
        rows.append({'model': model_name, 'eval_set': 'val', 'fold': fold, **mets})
        print(f'  fold {fold} | {model_name}: F1={mets["F1"]:.3f}  '
              f'P={mets["Precision"]:.3f}  R={mets["Recall"]:.3f}  '
              f'AUC={mets["AUC"]:.3f}  τ={tau:.3f}')

print('\n── Test (median fold-τ → test metrics) ──')
test = load_test()

# collect per-fold thresholds — never touches test data
fold_taus = {model_name: [] for model_name in FEATURES}
for fold in range(N_FOLDS):
    train = load_split(fold, 'train')
    for model_name, col in FEATURES.items():
        if col in train.columns:
            fold_taus[model_name].append(find_threshold(train, col))

for model_name, col in FEATURES.items():
    taus = fold_taus[model_name]
    if not taus:
        print(f'  WARNING: no thresholds found for {model_name} — skipping test eval')
        continue
    tau  = float(np.median(taus))
    mets = compute_metrics(test, col, tau)
    rows.append({'model': model_name, 'eval_set': 'test', 'fold': 'global', **mets})
    print(f'  {model_name}: F1={mets["F1"]:.3f}  P={mets["Precision"]:.3f}  '
          f'R={mets["Recall"]:.3f}  AUC={mets["AUC"]:.3f}  τ={tau:.3f}  '
          f'(fold τs: {[round(t, 3) for t in taus]})')

# ── save CSV ───────────────────────────────────────────────────────────────────
new_df = pd.DataFrame(rows)
these_models = new_df['model'].unique().tolist()

if os.path.exists(RESULTS_CSV):
    existing = pd.read_csv(RESULTS_CSV)
    existing = existing[~existing['model'].isin(these_models)]
    combined_df = pd.concat([existing, new_df], ignore_index=True)
else:
    combined_df = new_df

combined_df.to_csv(RESULTS_CSV, index=False)
print(f'\nsaved {RESULTS_CSV}')

# ── plotting ───────────────────────────────────────────────────────────────────
df_plot    = pd.read_csv(RESULTS_CSV)
all_models = df_plot['model'].unique().tolist()
palette    = plt.rcParams['axes.prop_cycle'].by_key()['color']
color_map  = {m: palette[i % len(palette)] for i, m in enumerate(all_models)}

PLOT_ROWS = [('val', 'CV (val)')]
n_models  = len(all_models)
col_width = max(4.0, 2.2 * n_models)
rng       = np.random.default_rng(42)
BREAK_GAP = 0.15
_TEXT_BOX = dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.75, edgecolor='none')


def _draw_break_marks(ax_top, ax_bot):
    d = 0.03
    for ax, y0 in [(ax_top, 0), (ax_bot, 1)]:
        kw = dict(transform=ax.transAxes, color='k', clip_on=False, lw=1.2)
        ax.plot((-d, +d), (y0 - 2*d, y0 + 2*d), **kw)
        ax.plot((1-d, 1+d), (y0 - 2*d, y0 + 2*d), **kw)


def _style(ax, show_xticks):
    ax.set_xticks(np.arange(n_models))
    ax.set_xlim(-0.6, n_models - 0.4)
    ax.tick_params(axis='y', labelsize=9)
    ax.grid(axis='y', lw=0.5, alpha=0.4, linestyle='--')
    ax.spines['right'].set_visible(False)
    if show_xticks:
        ax.set_xticklabels(all_models, fontsize=10, rotation=30, ha='right')
    else:
        ax.set_xticklabels([])
        ax.tick_params(bottom=False)


def _draw_boxes(ax, cv_sub, metric):
    for j, model in enumerate(all_models):
        vals = pd.to_numeric(cv_sub[cv_sub['model'] == model][metric]).dropna().values
        if len(vals) == 0:
            continue
        color = color_map[model]
        ax.boxplot(
            vals, positions=[j], widths=0.55,
            patch_artist=True,
            boxprops=dict(facecolor=color, alpha=0.35),
            medianprops=dict(color='black', lw=2),
            whiskerprops=dict(color=color, lw=1.2),
            capprops=dict(color=color, lw=1.2),
            showfliers=False,
        )
        jitter = rng.uniform(-0.13, 0.13, size=len(vals))
        ax.scatter(j + jitter, vals, color=color, s=35, zorder=5,
                   alpha=0.85, edgecolors='white', linewidths=0.4)
        median  = np.median(vals)
        y_label = np.max(vals) + 0.014
        ax.text(j, y_label, f'{median:.3f}', ha='center', va='bottom',
                fontsize=8, fontweight='bold', color='black', zorder=7,
                bbox=_TEXT_BOX)


def _draw_diamonds(ax, gt_sub, metric, label_above=False):
    for j, model in enumerate(all_models):
        row = gt_sub[gt_sub['model'] == model]
        if not len(row):
            continue
        val   = pd.to_numeric(row[metric]).values[0]
        color = color_map[model]
        ax.scatter([j], [val], marker='D', color=color, s=60, zorder=6,
                   edgecolors='black', linewidths=0.8)
        offset, va = (0.014, 'bottom') if label_above else (-0.020, 'top')
        ax.text(j, val + offset, f'{val:.3f}', ha='center', va=va,
                fontsize=7.5, color='black', style='italic', zorder=7,
                bbox=_TEXT_BOX)


fig = plt.figure(figsize=(col_width * len(METRICS), 6.0 * len(PLOT_ROWS)))
outer_gs = GridSpec(len(PLOT_ROWS), len(METRICS), figure=fig,
                    hspace=0.6, wspace=0.35, bottom=0.32, top=0.86)

for row_i, (eval_set, row_label) in enumerate(PLOT_ROWS):
    cv_sub = df_plot[df_plot['eval_set'] == eval_set]
    gt_sub = df_plot[df_plot['eval_set'] == 'test']

    for col_i, metric in enumerate(METRICS):
        all_cv  = pd.to_numeric(cv_sub[metric]).dropna().values
        gt_vals = pd.to_numeric(gt_sub[metric]).dropna().values \
                  if gt_sub is not None and len(gt_sub) else None

        needs_break = (gt_vals is not None and len(all_cv) > 0 and
                       all_cv.min() - gt_vals.max() > BREAK_GAP)

        cell = outer_gs[row_i, col_i]

        if needs_break:
            inner  = GridSpecFromSubplotSpec(2, 1, subplot_spec=cell,
                                             hspace=0.1, height_ratios=[3, 1])
            ax_top = fig.add_subplot(inner[0])
            ax_bot = fig.add_subplot(inner[1])

            _draw_boxes(ax_top, cv_sub, metric)
            ax_top.set_ylim(max(0.0, all_cv.min() - 0.04),
                            min(1.0, all_cv.max() + 0.12))
            ax_top.spines['bottom'].set_visible(False)
            ax_top.spines['top'].set_visible(False)
            _style(ax_top, show_xticks=False)

            _draw_diamonds(ax_bot, gt_sub, metric, label_above=True)
            ax_bot.set_ylim(max(0.0, gt_vals.min() - 0.08),
                            gt_vals.max() + 0.08)
            ax_bot.spines['top'].set_visible(False)
            ax_bot.spines['right'].set_visible(False)
            _style(ax_bot, show_xticks=True)

            _draw_break_marks(ax_top, ax_bot)

            if row_i == 0:
                ax_top.set_title(metric, fontsize=12, fontweight='bold', pad=8)
            if col_i == 0:
                ax_top.set_ylabel(row_label, fontsize=11, labelpad=6)

        else:
            ax = fig.add_subplot(cell)
            _draw_boxes(ax, cv_sub, metric)
            if gt_vals is not None and len(gt_vals) > 0:
                _draw_diamonds(ax, gt_sub, metric, label_above=False)

            all_shown = np.concatenate([all_cv, gt_vals]) \
                        if gt_vals is not None and len(gt_vals) > 0 else all_cv
            if len(all_shown) > 0:
                ax.set_ylim(max(0.0, all_shown.min() - 0.05),
                            min(1.0, all_shown.max() + 0.14))
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            _style(ax, show_xticks=True)

            if row_i == 0:
                ax.set_title(metric, fontsize=12, fontweight='bold', pad=8)
            if col_i == 0:
                ax.set_ylabel(row_label, fontsize=11, labelpad=6)

# legend
patches = [mpatches.Patch(color=color_map[m], alpha=0.7, label=m) for m in all_models]
diamond_handle = mlines.Line2D([], [], marker='D', color='grey', linestyle='None',
                               markersize=7, markeredgecolor='black', label='test (diamond)')
fig.legend(handles=patches + [diamond_handle], loc='lower center',
           ncol=n_models + 1, fontsize=10, frameon=False, bbox_to_anchor=(0.5, 0.06))

fig.suptitle(f'Threshold baseline — 5-fold CV  [{_csv_tag}]  [{_feat_tag}]',
             fontsize=13, fontweight='bold', y=0.96)

plt.savefig(OUT_PNG, dpi=150, bbox_inches='tight')
print(f'saved {OUT_PNG}')