# Code and Data repository

**Associated manuscript:**  
Eriksson & Schiller et al. – *Variations in the latitudinal diversity gradients of the ocean microbiome*

---

## Overview

This repository contains the code and processed data required to reproduce the analyses and figures presented in the associated manuscript.

The study integrates global marine metagenomic datasets, environmental climatologies, and habitat modeling to investigate patterns and drivers of microbial diversity across ocean basins, seasons, and taxonomic groups.

The repository includes:

- **R scripts** used for statistical analyses and figure generation
- **Processed data files** used in the analyses and visualizations
- **Sample metadata and contextual information**

Detailed descriptions of individual scripts and datasets are provided in the corresponding subdirectories.

---

## Repository Structure

The repository is organized by major analysis steps:

```
Code/
│
├── 1_Alpha_diversity_computation
├── 2_Extract_model_output
├── 3_Compute_ensembles_and_uncertainties
├── 4_Hotspots_diversity
├── 5_Wilcox_test
├── 6_Coastal_influence
├── Figures
└── Functions

Data/
│
├── samples_contextual.xlsx
├── species_diversity_taxonomic_rank.xlsx
└── Data_generated/
```

### Code

The `Code` directory contains R scripts used to:

- Compute **alpha diversity metrics** from metagenomic taxonomic profiles
- Extract and process **species distribution model (SDM) outputs**
- Calculate **ensemble means and uncertainty estimates**
- Quantify **latitudinal diversity gradients (LDGs)**
- Identify **diversity hotspots and associated environmental conditions**
- Perform **statistical tests** and sensitivity analyses
- Generate **all figures included in the manuscript**

### Data

The `Data` directory contains:

- **Sample metadata and contextual information** --> samples_contextual.xlsx
- **Species diversity (Richness, Shannon, Chao1) across taxonomic ranks** --> species_diversity_taxonomic_rank.xlsx
- **Derived datasets used for visualization and statistical analyses** --> directory "Data_generated"

Raw external datasets used in this study (e.g., metagenomic profiles and environmental climatologies) are referenced but not redistributed in this repository.

---

## Contact

For questions regarding the code or data, please contact:

**Jonas Schiller** - jonas.schiller@embl.de
**Dominic Eriksson** - dominic.eriksson@usys.ethz.ch
