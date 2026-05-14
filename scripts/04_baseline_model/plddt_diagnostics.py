"""
plddt_diagnostics.py — Q1 and Q2 diagnostic plots.

Supports two split types via --split-type flag:
  --split-type combined   (default) HLA + anchor splits, val/test eval sets
  --split-type dendrogram Single dendrogram split, val eval set only

Q1: pLDDT distributions per CV fold vs. full dataset, with optimal threshold marked.
    Checks whether fold splitting causes distribution inflation.

Q2: Per-allele binder/nonbinder ratio (class imbalance visualisation).

Usage examples:
  python plddt_diagnostics.py
  python plddt_diagnostics.py --split-type combined
  python plddt_diagnostics.py --split-type dendrogram

Outputs (saved next to this script):
  plddt_diagnostics_q1.png
  plddt_diagnostics_q2.png
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats import gaussian_kde
from sklearn.metrics import roc_curve

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description='pLDDT diagnostics Q1 + Q2')
parser.add_argument(
    '--split-type',
    choices=['combined', 'dendrogram'],
    default='combined',
    help=(
        '"combined"  → output_combined_splits, HLA-split-val + anchor-split-test rows (default)\n'
        '"dendrogram" → output_dendrogram_splits, single fold structure with val eval set'
    ),
)
args = parser.parse_args()

# ── Paths (edit BASE_DIR once; everything else derives from it) ────────────────
BASE_DIR = '/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/data'
COMBINED_PATH = f'{BASE_DIR}/output_combined_splits_from_resampled/combined.parquet'

if args.split_type == 'combined':
    COMBINED_PATH = f'{BASE_DIR}/output_combined_splits_from_resampled/combined.parquet'
    SPLITS_ROOT   = f'{BASE_DIR}/output_combined_splits_from_resampled/splits'
    # Each entry: (split_mode subfolder, eval parquet name, row label for plot)
    SPLIT_CONFIGS = [
        ('hla',    'val',  'HLA split — val sets'),
        ('anchor', 'test', 'Anchor split — test sets'),
    ]
    FOLD_RANGE = range(1, 6)   # folds 1–5
else:  # dendrogram
    COMBINED_PATH = f'{BASE_DIR}/output_dendro_splits_trainmask_val/combined.parquet'
    SPLITS_ROOT   = f'{BASE_DIR}/output_dendro_splits_trainmask_val/splits'
    # Dendrogram splits live directly under splits/fold_N/; no split_mode subfolder
    SPLIT_CONFIGS = [
        ('', 'val', 'Dendrogram split — val sets'),
    ]
    FOLD_RANGE = range(0,5)   # folds 0–4

SCRIPT_DIR    = os.path.dirname(os.path.abspath(__file__))
N_FOLDS       = 5
FEATURES      = {'pLDDT (pep)': 'pep_mean_plddt', 'pLDDT (anchor)': 'anchor_mean_plddt'}

BINDER_COLOR    = 'steelblue'
NONBINDER_COLOR = 'tomato'


# ── Helpers ────────────────────────────────────────────────────────────────────
def find_threshold(df, col):
    fpr, tpr, thresholds = roc_curve(df['label'], df[col])
    idx = np.argmin(np.sqrt(fpr**2 + (1 - tpr)**2))
    return thresholds[idx]
 
 
def load_split(split_mode, fold, part):
    """
    combined  → splits/<split_mode>/fold_<fold>/<part>.parquet
    dendrogram → splits/fold_<fold>/<part>.parquet  (split_mode is '')
    """
    if split_mode:
        path = f'{SPLITS_ROOT}/{split_mode}/fold_{fold}/{part}.parquet'
    else:
        path = f'{SPLITS_ROOT}/fold_{fold}/{part}.parquet'
    return pd.read_parquet(path)
 
 
def kde(x, grid):
    if len(x) < 2:
        return np.zeros_like(grid, dtype=float)
    k = gaussian_kde(x)
    return k(grid)
 
 
print(f'Split type : {args.split_type}')
print(f'Combined   : {COMBINED_PATH}')
print(f'Splits root: {SPLITS_ROOT}')
 
# ── Q1: pLDDT distributions ────────────────────────────────────────────────────
print('── Q1: pLDDT distributions ──')
combined = pd.read_parquet(COMBINED_PATH)
 
fold_palette = plt.cm.tab10(np.linspace(0, 0.45, len(FOLD_RANGE)))
 
fig, axes = plt.subplots(
    len(SPLIT_CONFIGS), len(FEATURES),
    figsize=(6 * len(FEATURES), 4.5 * len(SPLIT_CONFIGS)),
    squeeze=False,   # always 2-D axes array even for 1 row
)
fig.subplots_adjust(hspace=0.5, wspace=0.3, bottom=0.13, top=0.9)
 
for row_i, (split_mode, eval_part, row_label) in enumerate(SPLIT_CONFIGS):
    for col_i, (feat_name, col) in enumerate(FEATURES.items()):
        ax = axes[row_i, col_i]
 
        all_vals = combined[col].dropna().values
        x_grid   = np.linspace(np.percentile(all_vals, 0.5),
                               np.percentile(all_vals, 99.5), 400)
 
        # Full dataset reference (dashed)
        for label_val, color in [(1, BINDER_COLOR), (0, NONBINDER_COLOR)]:
            vals = combined.loc[combined['label'] == label_val, col].dropna().values
            k    = kde(vals, x_grid)
            if k.max() > 0:
                ax.plot(x_grid, k / k.max(), color=color, lw=2.2, ls='--', alpha=0.55)
 
        # Per-fold distributions + threshold
        for i, fold in enumerate(FOLD_RANGE):
            train  = load_split(split_mode, fold, 'train')
            eval_  = load_split(split_mode, fold, eval_part)
            tau    = find_threshold(train, col)
            color  = fold_palette[i]
 
            for label_val, ls in [(1, '-'), (0, ':')]:
                vals = eval_.loc[eval_['label'] == label_val, col].dropna().values
                k    = kde(vals, x_grid)
                if k.max() > 0:
                    ax.plot(x_grid, k / k.max(), color=color, lw=1.3, ls=ls, alpha=0.8)
 
            ax.axvline(tau, color=color, lw=1.0, alpha=0.6, zorder=3)
 
        ax.set_xlabel(feat_name, fontsize=10)
        ax.set_ylabel('Normalised density', fontsize=9)
        ax.set_title(f'{row_label}\n{feat_name}', fontsize=10, fontweight='bold')
        ax.spines[['top', 'right']].set_visible(False)
 
# Shared legend
handles = [
    mlines.Line2D([], [], color=BINDER_COLOR,    ls='--', lw=2.2, label='full dataset — binder'),
    mlines.Line2D([], [], color=NONBINDER_COLOR, ls='--', lw=2.2, label='full dataset — nonbinder'),
    mlines.Line2D([], [], color='grey', ls='-',  lw=1.5,  label='fold eval — binder'),
    mlines.Line2D([], [], color='grey', ls=':',  lw=1.5,  label='fold eval — nonbinder'),
    mlines.Line2D([], [], color='grey', ls='-',  lw=1.0,  alpha=0.6, label='fold threshold τ'),
]
fig.legend(handles=handles, loc='lower center', ncol=5,
           fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.01))
 
split_label = args.split_type.capitalize()
fig.suptitle(
    f'Q1 — pLDDT distributions: full dataset vs. CV folds  [{split_label} splits]',
    fontsize=12, fontweight='bold',
)
 
out_q1 = os.path.join(SCRIPT_DIR, 'plddt_diagnostics_q1.png')
plt.savefig(out_q1, dpi=150, bbox_inches='tight')
print(f'saved {out_q1}')
plt.close()
 
# ── Q2: per-allele binder/nonbinder ratio ──────────────────────────────────────
# Q2 uses only combined.parquet — identical for both split types
print('── Q2: per-allele imbalance ──')
 
allele_stats = (
    combined.groupby('allele')['label']
    .agg(binders=lambda x: (x == 1).sum(),
         nonbinders=lambda x: (x == 0).sum())
    .reset_index()
)
allele_stats['total'] = allele_stats['binders'] + allele_stats['nonbinders']
allele_stats['ratio'] = allele_stats['binders'] / allele_stats['total']
allele_stats = allele_stats.sort_values('total', ascending=False).reset_index(drop=True)
 
n = len(allele_stats)
x = np.arange(n)
overall_ratio = allele_stats['binders'].sum() / allele_stats['total'].sum()
 
fig2, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(max(14, n * 0.1), 8))
fig2.subplots_adjust(hspace=0.4, bottom=0.1, top=0.92)
 
# Stacked bar: binder + nonbinder counts
ax_top.bar(x, allele_stats['binders'],    color=BINDER_COLOR,    width=1.0, label='binders')
ax_top.bar(x, allele_stats['nonbinders'], color=NONBINDER_COLOR, width=1.0,
           bottom=allele_stats['binders'], label='nonbinders')
ax_top.set_ylabel('Sample count', fontsize=10)
ax_top.set_title('Per-allele sample counts (sorted by total, high → low)', fontsize=11, fontweight='bold')
ax_top.legend(fontsize=9, frameon=False)
ax_top.set_xticks([])
ax_top.spines[['top', 'right']].set_visible(False)
 
# Binder fraction per allele
ax_bot.bar(x, allele_stats['ratio'], color='mediumpurple', width=1.0, alpha=0.85)
ax_bot.axhline(overall_ratio, color='black', lw=1.5, ls='--',
               label=f'overall binder fraction ({overall_ratio:.2f})')
ax_bot.set_ylim(0, 1.05)
ax_bot.set_ylabel('Binder fraction', fontsize=10)
ax_bot.set_xlabel(f'Alleles  (n = {n}, sorted by total count)', fontsize=10)
ax_bot.set_xticks([])
ax_bot.legend(fontsize=9, frameon=False)
ax_bot.spines[['top', 'right']].set_visible(False)
 
fig2.suptitle(
    f'Q2 — Per-allele binder/nonbinder ratio  [{split_label} splits]',
    fontsize=12, fontweight='bold',
)
 
out_q2 = os.path.join(SCRIPT_DIR, 'plddt_diagnostics_q2.png')
plt.savefig(out_q2, dpi=150, bbox_inches='tight')
print(f'saved {out_q2}')
plt.close()