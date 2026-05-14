"""
plddt_per_allele_metrics.py
===========================
Computes per-allele classification metrics for the pLDDT baseline.

Default: dendrogram split only.  Use --split-type combined or both to include
the combined (HLA/anchor) splits.

Evaluation strategy:
  Val  — per fold: τ fit on train, per-allele metrics on that fold's val set
          → box plot per allele showing distribution across folds
  Test — median of per-fold τs, per-allele metrics on global test set
          → diamond markers overlaid on the box plots

Metrics: F1, Precision, Recall, AUC, AP, TPR@FPR10
Features: pLDDT (pep), pLDDT (anchor), PAE (pep)  [PAE negated: lower→higher]

Outputs saved next to this script:
  per_allele_results/per_allele_metrics.csv
  per_allele_results/per_allele_boxplots_<feature>_<split>_<eval>.png

Usage:
  python plddt_per_allele_metrics.py                              # all 4 dendrogram thresholds
  python plddt_per_allele_metrics.py --thresholds 70 80           # specific thresholds only
  python plddt_per_allele_metrics.py --split-type combined        # combined splits only
  python plddt_per_allele_metrics.py --split-type both            # combined + all dendrogram
"""

import argparse
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from sklearn.metrics import (
    roc_auc_score, average_precision_score,
    f1_score, precision_score, recall_score, roc_curve,
)

warnings.filterwarnings('ignore')

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument(
    '--split-type',
    choices=['combined', 'dendrogram', 'both'],
    default='dendrogram',
)
parser.add_argument(
    '--thresholds',
    nargs='+',
    choices=['50', '60', '70', '80'],
    default=['50', '60', '70', '80'],
    help='peptide-cluster thresholds to evaluate (dendrogram splits only)',
)
args = parser.parse_args()

# ── paths ──────────────────────────────────────────────────────────────────────
BASE_DIR   = '/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/data'
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR    = os.path.join('/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/pmhc-binding-prediction/results/01_baseline', 'per_allele_results')
os.makedirs(OUT_DIR, exist_ok=True)

COMBINED_SPLITS_ROOT    = f'{BASE_DIR}/output_combined_splits_from_resampled/splits'
DENDROGRAM_SPLITS_BASE  = f'{BASE_DIR}/combined_dendrogram_splits'

FEATURES = {
    'pLDDT (pep)':    'pep_mean_plddt',
    'pLDDT (anchor)': 'anchor_mean_plddt',
    'PAE (pep)':      'pae_pep',
}

PAE_COLS = {'pae_pep'}   # lower = better; negate so roc_curve sees higher = binder

N_FOLDS_COMBINED   = 5   # folds 1–5
N_FOLDS_DENDROGRAM = 5   # folds 0–4

METRICS = ['AUC', 'AP', 'TPR@FPR10', 'F1', 'Precision', 'Recall']

_TEXT_BOX = dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.75, edgecolor='none')


# ── helpers ────────────────────────────────────────────────────────────────────

def negate_pae(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in PAE_COLS:
        if col in df.columns:
            df[col] = -df[col]
    return df


def find_threshold(df, col):
    fpr, tpr, thresholds = roc_curve(df['label'], df[col])
    idx = np.argmin(np.sqrt(fpr**2 + (1 - tpr)**2))
    return float(thresholds[idx])


def compute_per_allele_metrics(df, col, tau, fold_label):
    """
    Returns a list of dicts, one per allele.
    AUC/AP/TPR@FPR10 are NaN when only one class is present.
    """
    rows = []
    for allele, grp in df.groupby('allele'):
        labels = grp['label'].values
        scores = grp[col].values
        preds  = (scores >= tau).astype(int)

        n_b  = int((labels == 1).sum())
        n_nb = int((labels == 0).sum())
        n    = len(labels)
        has_both = (n_b > 0 and n_nb > 0)

        if has_both:
            auc = float(roc_auc_score(labels, scores))
            ap  = float(average_precision_score(labels, scores))
            fpr_arr, tpr_arr, _ = roc_curve(labels, scores)
            mask         = fpr_arr <= 0.10
            tpr_at_fpr10 = float(tpr_arr[mask][-1]) if mask.any() else 0.0
        else:
            auc = ap = tpr_at_fpr10 = np.nan

        rows.append({
            'fold':         fold_label,
            'allele':       allele,
            'n_total':      n,
            'n_binders':    n_b,
            'n_nonbinders': n_nb,
            'single_class': not has_both,
            'threshold':    tau,
            'AUC':          auc,
            'AP':           ap,
            'TPR@FPR10':    tpr_at_fpr10,
            'F1':           float(f1_score(labels, preds, zero_division=0)),
            'Precision':    float(precision_score(labels, preds, zero_division=0)),
            'Recall':       float(recall_score(labels, preds, zero_division=0)),
        })
    return rows


def load_fold_part(root, fold, part):
    path = f'{root}/fold_{fold}/{part}.parquet'
    return negate_pae(pd.read_parquet(path)) if os.path.exists(path) else None


def load_global_test(root, flat=False):
    """flat=True: test.parquet directly in root; False: root/test/test.parquet"""
    path = f'{root}/test.parquet' if flat else f'{root}/test/test.parquet'
    return negate_pae(pd.read_parquet(path)) if os.path.exists(path) else None


# ── split configs ──────────────────────────────────────────────────────────────
CONFIGS = []

if args.split_type in ('combined', 'both'):
    CONFIGS += [
        dict(
            split_label      = 'combined_hla',
            root             = f'{COMBINED_SPLITS_ROOT}/hla',
            fold_range       = range(1, N_FOLDS_COMBINED + 1),
            has_global_test  = True,
            global_test_flat = True,
        ),
        dict(
            split_label      = 'combined_anchor',
            root             = f'{COMBINED_SPLITS_ROOT}/anchor',
            fold_range       = range(1, N_FOLDS_COMBINED + 1),
            has_global_test  = False,
            global_test_flat = False,
        ),
    ]

if args.split_type in ('dendrogram', 'both'):
    for thr in args.thresholds:
        CONFIGS.append(dict(
            split_label      = f'dendrogram_{thr}',
            root             = f'{DENDROGRAM_SPLITS_BASE}/{thr}/splits',
            fold_range       = range(0, N_FOLDS_DENDROGRAM),
            has_global_test  = True,
            global_test_flat = False,
        ))


# ── main loop ──────────────────────────────────────────────────────────────────
all_rows = []

for cfg in CONFIGS:
    split_label = cfg['split_label']
    root        = cfg['root']
    fold_range  = cfg['fold_range']
    print(f'\n{"="*60}')
    print(f'Split: {split_label}')

    # ── per-fold val metrics ──
    print('  Computing per-fold val metrics...')
    fold_taus_by_feat = {name: [] for name in FEATURES}

    for fold in fold_range:
        df_train = load_fold_part(root, fold, 'train')
        df_val   = load_fold_part(root, fold, 'val')
        if df_train is None or df_val is None:
            print(f'    WARNING: fold {fold} missing train or val — skipping')
            continue

        for feat_name, col in FEATURES.items():
            if col not in df_train.columns:
                continue
            tau = find_threshold(df_train, col)
            fold_taus_by_feat[feat_name].append(tau)

            feat_rows = compute_per_allele_metrics(df_val, col, tau, fold_label=fold)
            for r in feat_rows:
                all_rows.append({
                    'split':    split_label,
                    'eval_set': 'val',
                    'feature':  feat_name,
                    **r,
                })

    # ── global test metrics ──
    if cfg['has_global_test']:
        print('  Computing global test metrics (median τ)...')
        df_test = load_global_test(root, flat=cfg['global_test_flat'])
        if df_test is not None:
            for feat_name, col in FEATURES.items():
                taus = fold_taus_by_feat[feat_name]
                if not taus:
                    continue
                tau = float(np.median(taus))
                print(f'    {feat_name}: fold τs={[round(t,3) for t in taus]}  median τ={tau:.3f}')
                feat_rows = compute_per_allele_metrics(df_test, col, tau, fold_label='global')
                for r in feat_rows:
                    all_rows.append({
                        'split':    split_label,
                        'eval_set': 'global_test',
                        'feature':  feat_name,
                        **r,
                    })

# ── save CSV ───────────────────────────────────────────────────────────────────
results_df = pd.DataFrame(all_rows)
csv_path   = os.path.join(OUT_DIR, 'per_allele_metrics.csv')
results_df.to_csv(csv_path, index=False)
print(f'\nSaved: {csv_path}  ({len(results_df):,} rows)')

# summary
print('\nSummary (median AUC per split/eval/feature, alleles with both classes):')
val_df = results_df[(results_df['eval_set'] == 'val') & ~results_df['single_class']]
if not val_df.empty:
    summary = (
        val_df.groupby(['split', 'feature'])['AUC']
        .agg(['count', 'median', 'mean', 'std'])
        .round(3)
    )
    print(summary.to_string())


# ── plotting ───────────────────────────────────────────────────────────────────
print('\nGenerating per-allele box plots...')
rng = np.random.default_rng(42)
palette = plt.rcParams['axes.prop_cycle'].by_key()['color']


def _allele_order(val_sub, feat_name):
    """Alleles present in val data for this feature, sorted by median AUC desc."""
    feat_val = val_sub[val_sub['feature'] == feat_name]
    valid    = feat_val[~feat_val['single_class']]
    if valid.empty:
        return feat_val['allele'].unique().tolist()
    med_auc = valid.groupby('allele')['AUC'].median().sort_values(ascending=False)
    # include single-class alleles at the end
    single_alleles = [a for a in feat_val['allele'].unique() if a not in med_auc.index]
    return med_auc.index.tolist() + single_alleles


def _draw_boxes(ax, fold_data, alleles, metric, color):
    for j, allele in enumerate(alleles):
        vals = pd.to_numeric(
            fold_data[fold_data['allele'] == allele][metric]
        ).dropna().values
        if len(vals) == 0:
            continue
        ax.boxplot(
            vals, positions=[j], widths=0.6,
            patch_artist=True,
            boxprops=dict(facecolor=color, alpha=0.35),
            medianprops=dict(color='black', lw=2),
            whiskerprops=dict(color=color, lw=1.2),
            capprops=dict(color=color, lw=1.2),
            showfliers=False,
        )
        jitter = rng.uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(j + jitter, vals, color=color, s=25, zorder=5,
                   alpha=0.8, edgecolors='white', linewidths=0.3)
        median  = np.median(vals)
        y_label = np.max(vals) + 0.015
        ax.text(j, y_label, f'{median:.2f}', ha='center', va='bottom',
                fontsize=6, fontweight='bold', color='black', zorder=7,
                bbox=_TEXT_BOX)


def _draw_diamonds(ax, test_data, alleles, metric, color):
    for j, allele in enumerate(alleles):
        row = test_data[test_data['allele'] == allele]
        if row.empty:
            continue
        val = pd.to_numeric(row[metric]).values[0]
        if np.isnan(val):
            continue
        ax.scatter([j], [val], marker='D', color=color, s=50, zorder=6,
                   edgecolors='black', linewidths=0.7)


def _style_allele_ax(ax, alleles, show_xlabels=True):
    ax.set_xticks(range(len(alleles)))
    ax.set_xlim(-0.7, len(alleles) - 0.3)
    ax.tick_params(axis='y', labelsize=8)
    ax.grid(axis='y', lw=0.5, alpha=0.4, linestyle='--')
    ax.spines[['top', 'right']].set_visible(False)
    if show_xlabels:
        ax.set_xticklabels(alleles, fontsize=7, rotation=60, ha='right')
    else:
        ax.set_xticklabels([])
        ax.tick_params(bottom=False)


for cfg in CONFIGS:
    split_label = cfg['split_label']
    split_df    = results_df[results_df['split'] == split_label]
    val_df      = split_df[split_df['eval_set'] == 'val']
    test_df     = split_df[split_df['eval_set'] == 'global_test']

    for feat_name, col in FEATURES.items():
        feat_color = palette[list(FEATURES.keys()).index(feat_name) % len(palette)]
        feat_val   = val_df[val_df['feature'] == feat_name]
        feat_test  = test_df[test_df['feature'] == feat_name]

        if feat_val.empty:
            continue

        alleles   = _allele_order(val_df, feat_name)
        n_alleles = len(alleles)

        fig_width = max(12, 0.55 * n_alleles)
        fig, axes = plt.subplots(
            len(METRICS), 1,
            figsize=(fig_width, 3.5 * len(METRICS)),
            squeeze=False,
        )
        fig.subplots_adjust(hspace=0.55, left=0.06, right=0.98, top=0.94, bottom=0.08)

        for row_i, metric in enumerate(METRICS):
            ax = axes[row_i, 0]
            _draw_boxes(ax, feat_val, alleles, metric, feat_color)
            if not feat_test.empty:
                _draw_diamonds(ax, feat_test, alleles, metric, feat_color)

            all_vals = pd.to_numeric(feat_val[metric]).dropna().values
            test_vals = pd.to_numeric(feat_test[metric]).dropna().values if not feat_test.empty else np.array([])
            all_shown = np.concatenate([all_vals, test_vals]) if len(test_vals) > 0 else all_vals
            if len(all_shown) > 0:
                ax.set_ylim(max(0.0, all_shown.min() - 0.05),
                            min(1.05, all_shown.max() + 0.12))

            _style_allele_ax(ax, alleles, show_xlabels=(row_i == len(METRICS) - 1))
            ax.set_ylabel(metric, fontsize=9)

        # legend
        diamond_h = mlines.Line2D([], [], marker='D', color=feat_color,
                                   linestyle='None', markersize=6,
                                   markeredgecolor='black', label='test')
        box_h = plt.Rectangle((0, 0), 1, 1, fc=feat_color, alpha=0.4, label='val (CV folds)')
        fig.legend(handles=[box_h, diamond_h], loc='lower right',
                   fontsize=9, frameon=False, bbox_to_anchor=(0.99, 0.01))

        n_single = feat_val[feat_val['single_class']]['allele'].nunique()
        title = (f'Per-allele metrics — {split_label} | {feat_name}'
                 f'\n(alleles sorted by median AUC; {n_single} single-class alleles at end)')
        fig.suptitle(title, fontsize=11, fontweight='bold')

        safe_feat = feat_name.replace(' ', '_').replace('(', '').replace(')', '')
        png_path  = os.path.join(OUT_DIR, f'per_allele_boxplots_{safe_feat}_{split_label}.png')
        plt.savefig(png_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f'Saved: {png_path}')

print('\nDone.')