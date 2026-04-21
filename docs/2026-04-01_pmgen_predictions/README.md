# Structure Predictions via PMGen

## Date
2026-03-26 to 2026-04-01

## Analyst
Hasmig Aintablian (binders), Amir (non-binders)

## What we did
Ran AlphaFold2 structure predictions for all binder and non-binder peptide-MHC 
complexes using PMGen on GPU nodes (A100). Each chunk was submitted as a separate 
SLURM job on the `scc-gpu` partition. Up to 2 structures were predicted per 
peptide-MHC pair (seeds 0 and 1).

## Why we did it
3D structures of peptide-MHC complexes are needed for downstream structural analysis 
and model training. PMGen wraps AlphaFold2 for pMHC-specific structure prediction.

## Results
- Binders: 41 chunks predicted by Hasmig
- Non-binders: 31 chunks predicted by Amir
- Each ID produced up to 2 predicted structures (_0 and _1)
- Output files per structure: .pdb, plddt.npy, predicted_aligned_error.npy, ptm.npy

## Conclusion
All structures successfully predicted. Quality control (pLDDT QC) performed as 
next step to assess structural confidence.

## Storage
- Binder outputs: `/projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig/outputs/binder/`
- Non-binder outputs: `/projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/Amir/PMGen/outputs/hasmig/non_binders/chunk_outputs/`