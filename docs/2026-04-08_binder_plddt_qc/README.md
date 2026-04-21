# pLDDT Quality Control — Binder Structures

## Date
2026-04-08

## Analyst
Hasmig Aintablian

## What we did
Performed quality control on PMGen predicted binder peptide-MHC structures 
using per-residue pLDDT scores. For each structure, extracted the peptide region 
and anchor residues and computed their mean pLDDT. Up to 2 structures were predicted 
per peptide-MHC pair (seeds 0 and 1). Selected the best structure per allele-peptide 
pair based on the higher of the two scores (peptide or anchor pLDDT).

## Why we did it
pLDDT is a proxy for structural confidence. A threshold of 80 is commonly used —
structures below this are considered unreliable for downstream analysis.

## Results

### All structures (both seeds)
- Total structures: ~328,000
- Peptide mean pLDDT ≥ 80: 29.13%
- Anchor mean pLDDT ≥ 80: 43.11%

### Best structures (one per allele-peptide pair)
- Peptide mean pLDDT ≥ 80: 52.89%
- Anchor mean pLDDT ≥ 80: 70.39%

## Conclusion
Selecting the best structure per allele-peptide pair substantially improved confidence,
with anchor pLDDT ≥ 80 rising from 43% to 70%. The best structure CSV is used for 
all downstream analyses.

## Files
- `plots/` — pLDDT distribution plots (all and best structures)
- `plddt_means_binder.csv` — pLDDT means for all structures
- `plddt_means_binder_best.csv` — pLDDT means for best structures only