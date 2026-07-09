# Endothelial Instruction Competence Score (EICS)

## Project Overview

This project aims to develop the **Endothelial Instruction Competence Score (EICS)**, a computational metric that quantifies the developmental instructive capacity of endothelial cells in cardiac organoids.

Rather than evaluating vascularization solely by endothelial cell abundance, EICS integrates three functional dimensions of endothelial biology:

- Paracrine instruction
- Mechanical instruction
- Metabolic instruction

The score will be used to benchmark vascularized cardiac organoids against human fetal heart reference datasets and validated using spatial transcriptomics.

---

## Objectives

1. Construct the Endothelial Instruction Competence Score (EICS).
2. Benchmark publicly available cardiac organoid datasets.
3. Investigate endothelial–cardiomyocyte communication.
4. Validate EICS using human fetal heart spatial transcriptomics.

---

## Public datasets

### Human fetal heart
- Human fetal heart scRNA-seq (reference)
- Kanemaru et al., Nature (2023) (Spatial multiomics)

### Cardiac organoids
- Lewis-Israeli et al. (2021)
- Hofbauer et al. (2021)
- Drakhlis et al. (2021)
- Rossi et al. (2021)
- Volmert et al. (2023)
- Abilez et al. (2025)
- Pocock et al. (2025)

---

## Repository structure

```
data/
    raw/
    processed/

scripts/
    01_reference/
    02_EICS/
    03_organoid/
    04_cellchat/
    05_spatial/

gene_sets/
figures/
results/
manuscript/
```

---

## Project status

- [x] Repository initialized
- [x] RStudio Project created
- [ ] Download reference datasets
- [ ] Construct EICS gene sets
- [ ] Calculate EICS
- [ ] Cell-cell communication analysis
- [ ] Spatial validation