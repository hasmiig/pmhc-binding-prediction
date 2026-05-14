import pandas as pd
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt
import numpy as np

# load HLA test splits
binder_test    = pd.read_csv('/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/data/pmgen_outputs/binder/plddt_means_binder_best.csv')
nonbinder_test = pd.read_csv('/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/data/pmgen_outputs/non_binder/plddt_means_nonbinder_best.csv')

binder_test['label']    = 1
nonbinder_test['label'] = 0

df = pd.concat([binder_test, nonbinder_test], ignore_index=True)

# ROC curves
features = {
    'peptide mean pLDDT': 'pep_mean_plddt',
    'anchor mean pLDDT':  'anchor_mean_plddt',
}

fig, ax = plt.subplots(figsize=(5, 5))
for name, col in features.items():
    fpr, tpr, thresholds = roc_curve(df['label'], df[col])
    auc = roc_auc_score(df['label'], df[col])
    
    # optimal threshold — closest point to top-left corner
    distances = np.sqrt(fpr**2 + (1 - tpr)**2)
    optimal_idx = np.argmin(distances)
    optimal_threshold = thresholds[optimal_idx]

    line = ax.plot(fpr, tpr, label=f'{name} (AUC={auc:.3f})', lw=2)
    ax.scatter(fpr[optimal_idx], tpr[optimal_idx], color=line[0].get_color(), s=50, edgecolor='k', zorder=5, label=f'  τ={thresholds[optimal_idx]:.1f}, TPR={tpr[optimal_idx]:.3f}, FPR={fpr[optimal_idx]:.3f}')
    
    print(f'{name}:')
    print(f'  optimal threshold: {optimal_threshold:.2f}')
    print(f'  TPR: {tpr[optimal_idx]:.3f}')
    print(f'  FPR: {fpr[optimal_idx]:.3f}')
    print()

ax.plot([0,1],[0,1], 'k--', lw=0.8, label='random')
ax.set_xlabel('false positive rate')
ax.set_ylabel('true positive rate')
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig('stage1_roc.png', dpi=150)
print('saved stage1_roc.png')