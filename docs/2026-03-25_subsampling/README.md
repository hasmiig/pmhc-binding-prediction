# Subsampling — pMHC Binder and Non-binder Datasets

## Date
2026-03-25

## Analyst
Hasmig Aintablian

## What we did
Applied a two-phase iterative median sampling pipeline (`pmhc_sampling.py`) to balance 
the raw pMHC class I dataset from Amirreza's database (PMDb_2025_11_18_class1.parquet).

- **Binders**: Phase 1 only (MHC allele balancing). Anchor dominance in binders is 
  biological and should not be removed.
- **Non-binders**: Both phases (MHC allele balancing + anchor residue diversity balancing).

## Why we did it
The raw dataset is heavily imbalanced — some alleles have hundreds of thousands of 
peptides while others have only a handful. This would introduce bias into any downstream 
model. Phase 2 additionally balances anchor residue diversity within each allele for 
non-binders, since their anchor composition should not be biologically constrained.

## Results

### Binders
- Raw: 1,137,039 peptides across 459 alleles
- Post-Phase 1: 163,948 peptides across 459 alleles
- Cap: 1000 peptides per allele
- 128 alleles at cap, 331 below cap (too few peptides)

### Non-binders
- Raw: 42,727,479 peptides across 266 alleles
- Post-Phase 1: 215,108 peptides across 266 alleles
- Post-Phase 2: 120,123 peptides across 266 alleles
- Cap: 1000 peptides per allele
- Phase 2 reduced median from 1000 to 491 per allele, improving anchor diversity

## Conclusion
Binder dataset reduced from 1.1M to 163k peptides. Non-binder dataset reduced from 
42.7M to 120k peptides. Both datasets are now balanced across alleles and ready for 
PMGen input preparation.

## GitHub
Pipeline and full documentation available in the scripts/ directory.