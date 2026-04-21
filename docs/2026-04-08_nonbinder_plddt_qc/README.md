# pLDDT Quality Control — Non-binder Structures

## Date
2026-04-08

## Analyst
Hasmig Aintablian

## What we did
Performed quality control on PMGen predicted non-binder peptide-MHC structures 
using per-residue pLDDT scores. For each structure, extracted the peptide region 
and anchor residues and computed their mean pLDDT. Up to 2 structures were predicted 
per peptide-MHC pair (seeds 0 and 1).

## Why we did it
pLDDT is a proxy for structural confidence. A threshold of 80 is commonly used.
For non-binders we expect a slight shift of the pLDDT distribution to the left 
compared to binders, but not too much — AlphaFold2 has no capacity to distinguish 
binders from non-binders, so both should be reasonably well predicted structurally.
Both structures per peptide are kept for now as both may be useful as negatives.

## Results

### All structures (both seeds)
- Peptide mean pLDDT ≥ 80: 5.41%
- Anchor mean pLDDT ≥ 80: 8.39%

### Best structures (one per allele-peptide pair)
- Peptide mean pLDDT ≥ 80: 9.95%
- Anchor mean pLDDT ≥ 80: 14.69%

## Conclusion
As expected, non-binder structures show substantially lower pLDDT scores compared 
to binders. This shift to the left is expected since non-binders do not 
naturally fit the MHC groove, making AlphaFold2 less confident about their 
positioning. However, since AlphaFold2 cannot distinguish binders from non-binders, 
both structures per peptide are retained for downstream use as negative examples. 
Selection of best structures will be revisited if needed.

## Files
- `plots/` — pLDDT distribution plots (all structures)
- `plddt_means_nonbinder.csv` — pLDDT means for all structures
- `plddt_means_nonbinder_best.csv` — pLDDT means for best structures only