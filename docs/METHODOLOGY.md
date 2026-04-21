# Fine-tuning ProteinMPNN for pMHC Binding Prediction

## Overview

ProteinMPNN is the most widely used model for protein design tasks. However, similar to AlphaFold, it lacks sufficient observations of stable **peptide–MHC structures** in its training data. This document describes our methodology for fine-tuning ProteinMPNN on a large, high-quality dataset of predicted pMHC structures.

---

## Data Preparation

### Data Source

We downloaded eluted ligand and binding affinity data from the **Immune Epitope Database (IEDB)**. To reduce training noise from anchor complexity, we used **MHC-I data exclusively** for finetuning; however, the same methodology is also applicable to MHC-II.

We removed all pMHC pairs labeled as:
- Non-binders (ELlabel = 0)
- Low-affinity peptides (BA > 500 nM)

### MHC-Allele Balancing

Raw pMHC datasets are heavily imbalanced—some MHC alleles have thousands of peptides while others have only a handful. This imbalance introduces inductive bias into any downstream model trained on the data.

#### Iterative Median Sampling Strategy

We applied an **iterative median sampling** approach to balance allele representation:

1. Count the number of pMHC pairs observed for each MHC allele
2. Calculate the **median** of these counts across all sampled MHC alleles
3. Randomly sample that number of pMHC pairs for each allele
4. Remove all alleles with counts below the median from the next iteration
5. Continue iterations with a **sampling cap of 1,000 samples per allele**
   - Once an allele reaches 1,000 samples, it is excluded from subsequent sampling iterations
   - This ensures highly abundant alleles are capped while rarer alleles are fully represented

**Result:** A well-distributed MHC dataset with 426 MHC-I alleles and balanced representation across the allele frequency spectrum.

### Anchor Residue Analysis

We additionally analyzed anchor residue diversity as a diagnostic to assess whether further allele-specific balancing was needed.

For each peptide–MHC pair, we defined:
- **Anchor 1 (P2):** The second residue of the peptide
- **Anchor 2 (C-terminal):** The last residue of the peptide

For each MHC allele independently, we:
1. Calculated the frequency of each amino acid at both anchor positions
2. Applied the same iterative median sampling strategy with an **amino acid cap of 1,000**
3. Analyzed anchor combination diversity using 20×20 heatmaps (P2 × C-terminal)

**Finding:** No significant change in the total global anchor combination distribution after balancing. However, stratified analysis per MHC allele revealed **variation in anchor combination profiles, consistent with known allele-specific binding preferences**. 

**Decision:** Anchor residue balancing was **not applied** to the final dataset, as the allele-level balance already captured this biological variation.

---

## Structure Prediction with PMGen

### Input Preparation

We processed the MHC-allele-balanced samples using the **PMGen pipeline** in initial guess mode. 

For each peptide–MHC pair, anchor combinations were:
- **Enumerated without relying on NetMHCpan** (to reduce external dependencies)
- **Restricted to known MHC-I constraints:**
  - First anchor at position P1 or P2
  - Last anchor at the C-terminal residue
  - Consistent with established MHC-I binding constraints

### Structure Generation

For each peptide–MHC pair:
1. Generated **two structures** per pair
2. Retained the structure with **higher mean peptide pLDDT**
3. Excluded samples with **mean peptide pLDDT < 80** (confidence threshold)

**Initial yield:** 87,187 high-confidence structures spanning 426 MHC-I alleles

### Re-balancing After Prediction

To correct residual MHC allele bias introduced during structure prediction:
- Re-applied iterative median sampling to the predicted structures
- **Final dataset:** 63,817 structures (73% retention)

---

## Data Splitting for Evaluation

We evaluated generalizability using **two complementary splitting strategies**:

### Strategy 1: Allele-based Splitting

- Rank MHC alleles by frequency
- Hold out the **rarest 20%** (85 alleles out of 426) as an independent test set
- Divide remaining alleles into **5 folds** for cross-validation

**Use case:** Evaluate model ability to generalize to unseen MHC alleles.

### Strategy 2: Anchor-Combination-based Splitting

- Divide folds based on **unique anchor combinations** (P2 × C-terminal pairs)
- Each observed anchor combination is held out exactly once across 5 rotating folds:
  - **Test set:** One chunk of combinations
  - **Validation set:** Next chunk of combinations
  - **Training set:** Remaining three chunks

**Use case:** Evaluate model ability to generalize to unseen anchor residue combinations within known alleles.

---

## Summary Statistics

| Metric | Value |
|--------|-------|
| Initial IEDB MHC-I pairs | ~500,000 |
| After binding filtering | ~250,000 |
| After allele balancing | ~100,000 |
| After structure prediction (pLDDT ≥ 80) | 87,187 |
| After final re-balancing | **63,817** |
| Unique MHC-I alleles | 426 |
| Sampling cap per allele | 1,000 |
| Anchor AA cap per allele | 1,000 |
| Train/Val/Test split | 5-fold (60% / 20% / 20%) |

---

## Quality Control Metrics

### pLDDT Confidence

Mean peptide pLDDT was used as the primary quality metric:
- **Threshold:** ≥ 80 indicates a high-confidence prediction
- **Rationale:** Peptide pLDDT is a strong proxy for sampling quality and stable local geometry

### Allele-Specific Anchor Profiles

We validated that training data preserved known allele-specific MHC-peptide binding preferences through per-allele anchor analysis. This confirmed that balancing did not over-homogenize the data in a way that would lose biological signal.

---

## Downstream Tasks

### Fine-tuning ProteinMPNN

The final curated dataset of 63,817 structures spanning 426 MHC-I alleles is used to fine-tune ProteinMPNN via:

1. **Sequence generation:** Generate sequences conditioned on the pMHC complex backbone
2. **Binding prediction:** Fine-tune a classification head to predict binding affinity or probability
3. **Structure-guided design:** Improve designs that stabilize the pMHC interface

### Reproducibility

All scripts used for data curation are provided in the `scripts/` directory. The complete data preparation pipeline is modular and can be adapted for:
- MHC-II peptides (change filtering criteria)
- Different allele balancing caps (adjust `SAMPLE_CAP`)
- Alternative quality thresholds (adjust `pLDDT_MIN`)

---

## References

- IEDB: Vita, R., et al. "The immune epitope database (IEDB): 2018 update." Nucleic Acids Research (2019)
- AlphaFold: Jumper, J., et al. "Highly accurate protein structure prediction with AlphaFold2." Nature (2021)
- ProteinMPNN: Dauparas, J., et al. "Robust deep learning-based protein sequence design using proteinmpnn." bioRxiv (2022)
