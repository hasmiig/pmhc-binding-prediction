# pMHC Binding Prediction — Data Methodology

## Overview

This document describes the data curation methodology for an early-stage pMHC class I **binding classifier** — a model that predicts whether a given peptide binds a given MHC allele. Unlike a sequence design model, a classifier requires both positive (binder) and negative (non-binder) examples. The pipeline prepares balanced, structure-predicted, QC-filtered datasets for training and evaluating this classifier.

---

## Data Preparation

### Data Source

We downloaded eluted ligand and binding affinity data from the **Immune Epitope Database (IEDB)**. We restrict to **MHC class I** data.

We retain:
- **Binders** — eluted ligands (`ELlabel = 1`) or high-affinity measurements (BA ≤ 500 nM), labeled `assigned_label = 1.0`
- **Non-binders** — peptides confirmed as non-binding (`assigned_label = 0.0`)

Both classes are kept and processed independently. Non-binders are essential as negative examples for the classifier.

### MHC-Allele Balancing

Raw pMHC datasets are heavily imbalanced — some MHC alleles have thousands of peptides while others have only a handful. This imbalance introduces inductive bias into any downstream model.

#### Iterative Median Sampling Strategy

We applied an **iterative median sampling** approach to balance allele representation within each class independently:

1. Count the number of pMHC pairs observed for each MHC allele
2. Calculate the **median** of these counts across all alleles
3. Randomly sample that number of pMHC pairs for each allele
4. Remove all alleles with counts below the median from the next iteration
5. Continue with a **sampling cap of 1,000 samples per allele**
   - Once an allele reaches 1,000 samples, it is excluded from further iterations
   - This ensures highly abundant alleles are capped while rarer alleles are fully represented

**Result:** A well-distributed dataset with balanced allele representation.

### Anchor Residue Analysis and Balancing

We analyzed anchor residue diversity at position **P2** (anchor 1) and the **C-terminal** position (anchor 2) to assess whether further allele-specific balancing was needed.

**Binders — Phase 1 only:**  
Anchor dominance in binders reflects true allele-specific binding preferences (e.g. HLA-A*02:01 strongly prefers leucine at P2). Forcing artificial anchor diversity would destroy this biological signal. We apply allele balancing only and preserve anchor profiles.

**Non-binders — both phases:**  
Anchor dominance in non-binders is less biologically constrained and more likely to reflect data collection artifacts or sampling biases. We additionally apply Phase 2 (anchor residue balancing within each allele) to reduce this noise.

---

## Structure Prediction with PMGen

### Input Preparation

We processed the allele-balanced samples using the **PMGen pipeline** in initial guess mode.

For each peptide–MHC pair, anchor combinations were:
- **Enumerated without relying on NetMHCpan** (to reduce external dependencies)
- **Restricted to known MHC-I constraints:**
  - First anchor at position P1 or P2
  - Last anchor at the C-terminal residue

### Structure Generation

For each peptide–MHC pair:
1. Generated **two structures** per pair
2. Retained the structure with **higher mean peptide pLDDT**
3. Excluded samples with **mean peptide pLDDT < 80** (confidence threshold)

**Initial yield:** 87,187 high-confidence structures spanning 426 MHC-I alleles (binders)

### Re-balancing After Prediction

Structure prediction introduces residual allele bias (some alleles fail more often). We re-apply iterative median sampling to the predicted structures.

**Final binder dataset:** 63,817 structures (73% retention)

Non-binder structures follow the same pipeline and QC independently.

---

## Data Splitting for Evaluation

We evaluate generalizability using **two complementary splitting strategies**. Both are applied to the combined (binder + non-binder) dataset.

### Strategy 1: Allele-based Splitting

- Rank MHC alleles by frequency
- Hold out the **rarest 20%** (85 alleles out of 426) as an independent test set
- Divide remaining alleles into **5 folds** for cross-validation

**Use case:** Evaluate generalization to unseen MHC alleles.

### Strategy 2: Anchor-Combination-based Splitting

- Divide folds based on **unique anchor combinations** (P2 × C-terminal pairs)
- Each observed anchor combination is held out exactly once across 5 rotating folds:
  - **Test set:** one chunk of combinations
  - **Validation set:** next chunk
  - **Training set:** remaining three chunks

**Use case:** Evaluate generalization to unseen anchor combinations within known alleles.

---

## Summary Statistics

| Metric | Value |
|--------|-------|
| Initial IEDB MHC-I pairs (binders) | ~250,000 |
| After allele balancing (binders) | ~100,000 |
| After structure prediction, pLDDT ≥ 80 (binders) | 87,187 |
| After final re-balancing (binders) | **63,817** |
| Unique MHC-I alleles | 426 |
| Sampling cap per allele | 1,000 |
| Train/Val/Test split | 5-fold (60% / 20% / 20%) |

---

## Quality Control

### pLDDT Confidence

Mean peptide pLDDT is the primary quality metric:
- **Threshold:** ≥ 80 — high-confidence structure
- **Rationale:** Peptide pLDDT is a strong proxy for sampling quality and stable local geometry

### Allele-Specific Anchor Profiles

We validated that training data preserved known allele-specific binding preferences through per-allele anchor analysis. This confirmed that balancing did not over-homogenize the binder data.

For non-binders, we separately confirmed that Phase 2 anchor balancing reduced artifactual anchor dominance without eliminating biological variation.

---

## Downstream Use

The final curated dataset (balanced binders + non-binders, pLDDT-filtered, structure-annotated) is used to train and evaluate a **binary binding classifier** that predicts whether a peptide binds a given MHC-I allele. This is an early-stage project; model architecture and training details are outside the scope of this document.

All data curation scripts are provided in `scripts/`. The pipeline is modular and can be adapted for:
- MHC-II peptides (adjust filtering criteria)
- Different allele balancing caps (adjust `SAMPLE_CAP`)
- Alternative pLDDT thresholds (adjust `--plddt_threshold`)

---

## References

- IEDB: Vita, R., et al. "The immune epitope database (IEDB): 2018 update." *Nucleic Acids Research* (2019)
- PMGen: structure prediction pipeline for pMHC complexes
- AlphaFold: Jumper, J., et al. "Highly accurate protein structure prediction with AlphaFold2." *Nature* (2021)
