"""
plddt_per_allele_plots.py
=========================
Four analytical plots from per_allele_metrics.csv
(output of plddt_per_allele_metrics.py).

  1. n_binders vs AUC scatter      — performance vs data quantity per allele
  2. Cross-threshold stability      — allele AUC heatmap across dendrogram thresholds
  3. Feature agreement scatter      — pairwise AUC comparisons between features
  4. Allele metric heatmap          — alleles × (feature × metric) dense summary
  5. Evaluability heatmap           — per-allele evaluable/single-class/not-in-val across thresholds

Usage:
  python plddt_per_allele_plots.py \\
      --csv  .../per_allele_results/per_allele_metrics.csv \\
      --out-dir .../per_allele_results      # optional; defaults to CSV directory
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from scipy import stats

DENDRO_THRESHOLDS = ['50', '60', '70', '80']
FEATURES = ['pLDDT (pep)', 'pLDDT (anchor)', 'PAE (pep)']
METRICS  = ['AUC', 'AP', 'TPR@FPR10', 'F1', 'Precision', 'Recall']

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument('--csv', required=True, help='Path to per_allele_metrics.csv')
parser.add_argument('--out-dir', default=None,
                    help='Output directory (default: same directory as CSV)')
args = parser.parse_args()

OUT_DIR = args.out_dir or os.path.dirname(os.path.abspath(args.csv))
os.makedirs(OUT_DIR, exist_ok=True)

# ── load ───────────────────────────────────────────────────────────────────────
raw = pd.read_csv(args.csv)
for m in METRICS:
    raw[m] = pd.to_numeric(raw[m], errors='coerce')
raw['n_binders'] = pd.to_numeric(raw['n_binders'], errors='coerce')

# val rows with both classes present
val = raw[(raw['eval_set'] == 'val') & ~raw['single_class']].copy()
val['thr'] = val['split'].str.extract(r'dendrogram_(\d+)')[0]
dendro_val = val[val['split'].str.startswith('dendrogram_')].copy()

if dendro_val.empty:
    print('WARNING: no dendrogram val rows found in CSV — exiting')
    raise SystemExit(1)

# ── colours ────────────────────────────────────────────────────────────────────
palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
thr_colors = {t: palette[i % len(palette)] for i, t in enumerate(DENDRO_THRESHOLDS)}
_TEXT_BOX  = dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.75, edgecolor='none')

# ── per-allele summary: median metric across folds, sum n_binders ─────────────
_agg = {m: 'median' for m in METRICS}
_agg['n_binders'] = 'sum'
allele_sum = (
    dendro_val
    .groupby(['split', 'thr', 'feature', 'allele'])[METRICS + ['n_binders']]
    .agg(_agg)
    .reset_index()
)


# ══════════════════════════════════════════════════════════════════════════════
# Plot 1 — n_binders vs AUC scatter
# ══════════════════════════════════════════════════════════════════════════════
print('Plot 1: n_binders vs AUC...')

_n_max = allele_sum['n_binders'].max()

fig, axes = plt.subplots(1, len(FEATURES),
                         figsize=(5.5 * len(FEATURES), 5), sharey=False)

for col_i, feat in enumerate(FEATURES):
    ax   = axes[col_i]
    fsub = allele_sum[allele_sum['feature'] == feat]

    for thr in DENDRO_THRESHOLDS:
        pts = fsub[fsub['thr'] == thr].dropna(subset=['n_binders', 'AUC'])
        if pts.empty:
            continue
        s_vals = 8 + 150 * np.sqrt(pts['n_binders'] / _n_max)
        ax.scatter(pts['n_binders'], pts['AUC'],
                   color=thr_colors[thr], alpha=0.65, s=s_vals,
                   edgecolors='none', zorder=3)

    # ρ computed only on alleles with enough binders to give reliable AUC
    reliable = fsub[fsub['n_binders'] >= 10].dropna(subset=['n_binders', 'AUC'])
    n_small  = int((fsub['n_binders'] < 10).sum())
    if len(reliable) > 4:
        r, _ = stats.spearmanr(np.log1p(reliable['n_binders']), reliable['AUC'])
        ax.text(0.97, 0.04,
                f'ρ = {r:.2f}  (n≥10 binders)\n{n_small} pts with n<10 (unreliable)',
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=7.5, bbox=_TEXT_BOX)

    # size legend on last panel
    if col_i == len(FEATURES) - 1:
        for n_ref in [10, 100, 500]:
            s_ref = 8 + 150 * np.sqrt(n_ref / _n_max)
            ax.scatter([], [], color='grey', alpha=0.6, s=s_ref,
                       edgecolors='none', label=f'n = {n_ref}')
        ax.legend(fontsize=7, frameon=False,
                  title='n binders\n(point size)', title_fontsize=7,
                  loc='upper left')

    ax.set_xscale('log')
    ax.set_xlabel('Val binders — sum across folds (log scale)', fontsize=9)
    ax.set_ylabel('Median AUC (across folds)', fontsize=9)
    ax.set_title(feat, fontsize=10, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', lw=0.4, alpha=0.35, linestyle='--')

handles = [
    plt.Line2D([0], [0], marker='o', color='w',
               markerfacecolor=thr_colors[t], markersize=7, label=f'thr {t}')
    for t in DENDRO_THRESHOLDS
]
fig.legend(handles=handles, title='pep-clust threshold', title_fontsize=7,
           fontsize=8, frameon=False,
           loc='upper center', bbox_to_anchor=(0.5, 1.0),
           ncol=len(DENDRO_THRESHOLDS))
fig.suptitle('Data quantity vs performance per allele\n'
             '(dendrogram val — median AUC across folds)',
             fontsize=11, fontweight='bold', y=1.07)
plt.tight_layout(rect=[0, 0, 1, 0.90])
out1 = os.path.join(OUT_DIR, 'plot1_nbinders_vs_auc.png')
plt.savefig(out1, dpi=150, bbox_inches='tight')
plt.close()
print(f'  Saved: {out1}')


# ══════════════════════════════════════════════════════════════════════════════
# Plot 2 — Cross-threshold stability heatmap
# ══════════════════════════════════════════════════════════════════════════════
print('Plot 2: cross-threshold stability heatmap...')

# unified allele order: mean AUC across all features and thresholds
unified_allele_order = (
    allele_sum.groupby('allele')['AUC']
    .mean()
    .sort_values(ascending=False)
    .index.tolist()
)

pivots = {}
for feat in FEATURES:
    fsub  = allele_sum[allele_sum['feature'] == feat][['allele', 'thr', 'AUC']]
    pivot = fsub.pivot_table(index='allele', columns='thr', values='AUC')
    pivot = pivot.reindex(columns=[t for t in DENDRO_THRESHOLDS if t in pivot.columns])
    pivot = pivot[pivot.notna().sum(axis=1) >= 1]
    # apply unified order, keeping only alleles present in this feature
    ordered = [a for a in unified_allele_order if a in pivot.index]
    pivot = pivot.loc[ordered]
    pivots[feat] = pivot

n_alleles_max = max(len(p) for p in pivots.values()) if pivots else 0
fig_h = max(6, n_alleles_max * 0.22)

fig, axes = plt.subplots(1, len(FEATURES),
                         figsize=(4 * len(FEATURES), fig_h),
                         sharey=False)

for col_i, feat in enumerate(FEATURES):
    ax    = axes[col_i]
    pivot = pivots[feat]
    im    = ax.imshow(pivot.values, aspect='auto', cmap='RdYlGn',
                      vmin=0.4, vmax=1.0, interpolation='nearest')

    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([f'thr {t}' for t in pivot.columns],
                       fontsize=8, rotation=30, ha='right')
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=6)
    ax.set_title(feat, fontsize=10, fontweight='bold')

    # annotate cells
    for ri in range(len(pivot)):
        for ci in range(len(pivot.columns)):
            v = pivot.values[ri, ci]
            if not np.isnan(v):
                ax.text(ci, ri, f'{v:.2f}', ha='center', va='center',
                        fontsize=5, color='black' if 0.45 < v < 0.85 else 'white')

    plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02,
                 label='Median AUC' if col_i == len(FEATURES) - 1 else '')

fig.suptitle('Per-allele AUC stability across peptide-cluster thresholds\n'
             '(dendrogram val — median across folds; sorted by mean AUC)',
             fontsize=11, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.97])
out2 = os.path.join(OUT_DIR, 'plot2_threshold_stability.png')
plt.savefig(out2, dpi=150, bbox_inches='tight')
plt.close()
print(f'  Saved: {out2}')


# ══════════════════════════════════════════════════════════════════════════════
# Plot 3 — Feature agreement scatter (pairwise AUC)
# ══════════════════════════════════════════════════════════════════════════════
print('Plot 3: feature agreement scatter...')

feat_pivot = (
    allele_sum
    .pivot_table(index=['allele', 'split', 'thr'], columns='feature', values='AUC')
    .reset_index()
)

pairs = [
    ('pLDDT (pep)',    'pLDDT (anchor)'),
    ('pLDDT (pep)',    'PAE (pep)'),
    ('pLDDT (anchor)', 'PAE (pep)'),
]

fig, axes = plt.subplots(1, len(pairs),
                         figsize=(5.5 * len(pairs), 5), sharey=False)

for col_i, (fx, fy) in enumerate(pairs):
    ax  = axes[col_i]
    sub = feat_pivot[[fx, fy, 'thr']].dropna()

    for thr in DENDRO_THRESHOLDS:
        pts = sub[sub['thr'] == thr]
        if pts.empty:
            continue
        ax.scatter(pts[fx], pts[fy],
                   color=thr_colors[thr], alpha=0.55, s=30,
                   label=f'thr {thr}', edgecolors='none', zorder=3)

    # identity line
    ax.plot([0, 1], [0, 1], 'k--', lw=0.8, alpha=0.35, zorder=1)

    all_pts = sub.dropna()
    if len(all_pts) > 4:
        r_sp, _ = stats.spearmanr(all_pts[fx], all_pts[fy])
        r_pe, _ = stats.pearsonr(all_pts[fx], all_pts[fy])
        ax.text(0.03, 0.97,
                f'Spearman ρ = {r_sp:.2f}\nPearson r  = {r_pe:.2f}',
                transform=ax.transAxes, ha='left', va='top',
                fontsize=8, bbox=_TEXT_BOX)

    ax.set_xlabel(fx, fontsize=9)
    ax.set_ylabel(fy, fontsize=9)
    ax.set_title(f'{fx}\nvs {fy}', fontsize=9, fontweight='bold')
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(lw=0.4, alpha=0.3, linestyle='--')

handles = [
    plt.Line2D([0], [0], marker='o', color='w',
               markerfacecolor=thr_colors[t], markersize=7, label=f'thr {t}')
    for t in DENDRO_THRESHOLDS
]
fig.legend(handles=handles, title='pep-clust threshold', title_fontsize=7,
           fontsize=8, frameon=False,
           loc='upper center', bbox_to_anchor=(0.5, 1.0),
           ncol=len(DENDRO_THRESHOLDS))
fig.suptitle('Feature agreement — per-allele median AUC (dendrogram val)',
             fontsize=11, fontweight='bold', y=1.07)
plt.tight_layout(rect=[0, 0, 1, 0.90])
out3 = os.path.join(OUT_DIR, 'plot3_feature_agreement.png')
plt.savefig(out3, dpi=150, bbox_inches='tight')
plt.close()
print(f'  Saved: {out3}')


# ══════════════════════════════════════════════════════════════════════════════
# Plot 4 — Allele metric heatmap (alleles × feature × metric)
# ══════════════════════════════════════════════════════════════════════════════
print('Plot 4: allele metric heatmap...')

# median across all folds AND thresholds
heat_base = (
    dendro_val
    .groupby(['feature', 'allele'])[METRICS]
    .median()
    .reset_index()
)

# build wide frame: one column per (feature, metric)
feat_frames = []
for feat in FEATURES:
    fsub = heat_base[heat_base['feature'] == feat].set_index('allele')[METRICS]
    fsub.columns = pd.MultiIndex.from_tuples([(feat, m) for m in METRICS])
    feat_frames.append(fsub)

heat_df = pd.concat(feat_frames, axis=1).dropna(how='all')
# sort by pLDDT (pep) AUC descending
sort_key = ('pLDDT (pep)', 'AUC')
if sort_key in heat_df.columns:
    heat_df = heat_df.sort_values(sort_key, ascending=False)

n_alleles = len(heat_df)
n_cols    = len(heat_df.columns)            # 18 = 3 features × 6 metrics

fig, ax = plt.subplots(figsize=(n_cols * 0.85, max(8, n_alleles * 0.22)))
im = ax.imshow(heat_df.values.astype(float), aspect='auto',
               cmap='RdYlGn', vmin=0, vmax=1, interpolation='nearest')

# x-axis: metric names (bottom row of multi-index)
ax.set_xticks(range(n_cols))
ax.set_xticklabels(heat_df.columns.get_level_values(1),
                   fontsize=7, rotation=45, ha='right')

# y-axis: allele names
ax.set_yticks(range(n_alleles))
ax.set_yticklabels(heat_df.index, fontsize=6)

# vertical dividers between feature groups
for g in range(1, len(FEATURES)):
    ax.axvline(g * len(METRICS) - 0.5, color='white', lw=2)

# feature group labels above the columns
for g, feat in enumerate(FEATURES):
    x_center = g * len(METRICS) + len(METRICS) / 2 - 0.5
    ax.text(x_center, 1.02, feat, ha='center', va='bottom',
            fontsize=8, fontweight='bold',
            transform=ax.get_xaxis_transform())

plt.colorbar(im, ax=ax, fraction=0.01, pad=0.01, label='Metric value')
fig.suptitle(
    'Per-allele metric heatmap — all dendrogram thresholds, val set\n'
    '(median across folds & thresholds; sorted by pLDDT(pep) AUC)',
    fontsize=10, fontweight='bold',
)
fig.subplots_adjust(top=0.93, bottom=0.08, left=0.12, right=0.96)
out4 = os.path.join(OUT_DIR, 'plot4_allele_metric_heatmap.png')
plt.savefig(out4, dpi=150, bbox_inches='tight')
plt.close()
print(f'  Saved: {out4}')

# ══════════════════════════════════════════════════════════════════════════════
# Plot 5 — Evaluability per allele across thresholds
# ══════════════════════════════════════════════════════════════════════════════
print('Plot 5: evaluability heatmap...')

# single_class is data-dependent (not feature-dependent) — use one feature as ref
ref = raw[(raw['eval_set'] == 'val') & (raw['feature'] == 'pLDDT (pep)')].copy()
ref['thr'] = ref['split'].str.extract(r'dendrogram_(\d+)')[0]
ref = ref[ref['thr'].notna()]

thr_list = [t for t in DENDRO_THRESHOLDS if t in ref['thr'].unique()]

# status per (allele, threshold): 2=evaluable, 1=single-class, 0=not-in-val
ref['status'] = np.where(~ref['single_class'], 2, 1)
status = (
    ref.pivot_table(index='allele', columns='thr', values='status', aggfunc='max')
    .reindex(columns=thr_list)
    .fillna(0)
    .astype(int)
)

# sort: evaluable count desc, then mean AUC desc
n_eval = (status == 2).sum(axis=1)
mean_auc_all = allele_sum.groupby('allele')['AUC'].mean()
status['_n_eval']   = n_eval
status['_mean_auc'] = status.index.map(mean_auc_all).fillna(0)
status = status.sort_values(['_n_eval', '_mean_auc'], ascending=[False, False])
status = status.drop(columns=['_n_eval', '_mean_auc'])

n_alleles_5 = len(status)
fig_h5 = max(8, n_alleles_5 * 0.18)

fig5, (ax_heat, ax_hist) = plt.subplots(
    1, 2,
    figsize=(7, fig_h5),
    gridspec_kw={'width_ratios': [4, 1]},
)

# ── heatmap ──
ev_cmap = ListedColormap(['#cccccc', '#ff7f0e', '#2ca02c'])   # grey / orange / green
im5 = ax_heat.imshow(status.values, aspect='auto', cmap=ev_cmap,
                     vmin=0, vmax=2, interpolation='nearest')

ax_heat.set_xticks(range(len(thr_list)))
ax_heat.set_xticklabels([f'thr {t}' for t in thr_list], fontsize=8)
ax_heat.set_yticks(range(n_alleles_5))
ax_heat.set_yticklabels(status.index, fontsize=5.5)
ax_heat.set_title('Per-allele evaluability\n(1 val fold per threshold)',
                  fontsize=10, fontweight='bold')

legend_patches = [
    mpatches.Patch(color='#2ca02c', label='evaluable (both classes)'),
    mpatches.Patch(color='#ff7f0e', label='single-class in val'),
    mpatches.Patch(color='#cccccc', label='not in val (test / train)'),
]
ax_heat.legend(handles=legend_patches, fontsize=7, frameon=True,
               loc='lower left', bbox_to_anchor=(0, -0.06),
               ncol=1, borderpad=0.5)

# ── histogram: how many thresholds each allele is evaluable in ──
eval_counts = (status == 2).sum(axis=1).value_counts().sort_index()
ax_hist.barh(eval_counts.index.astype(str), eval_counts.values,
             color='#2ca02c', alpha=0.75, edgecolor='white')
for y, v in zip(eval_counts.index.astype(str), eval_counts.values):
    ax_hist.text(v + 0.5, y, str(v), va='center', fontsize=8)
ax_hist.set_xlabel('# alleles', fontsize=9)
ax_hist.set_ylabel('# thresholds\nevaluable', fontsize=9)
ax_hist.spines[['top', 'right']].set_visible(False)
ax_hist.set_title('Distribution', fontsize=9, fontweight='bold')

fig5.suptitle('Allele evaluability across dendrogram splits',
              fontsize=11, fontweight='bold')
fig5.subplots_adjust(top=0.96, bottom=0.10, left=0.18, right=0.97, wspace=0.35)
out5 = os.path.join(OUT_DIR, 'plot5_evaluability.png')
plt.savefig(out5, dpi=150, bbox_inches='tight')
plt.close()
print(f'  Saved: {out5}')

print('\nDone.')
