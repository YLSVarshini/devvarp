# DevVarP: Developmental Variant Prioritization

## Pathogenic Non-Coding Variants in Congenital-Onset Pediatric Rare Diseases Are Depleted from Annotated Fetal Regulatory Elements: Evidence for a Temporal Annotation Gap

**Author:** Lakshmi Sai Varshini Yedavalli
**Contact:** YLSVarshini@gmail.com

---

## Overview

This repository contains the complete computational pipeline for the DevVarP study — the first systematic investigation of fetal regulatory element enrichment in pathogenic non-coding variants stratified by pediatric developmental onset.

### Key Finding

Pathogenic non-coding variants in congenital-onset pediatric rare diseases are significantly depleted from annotated fetal chromatin accessibility sites (OR=0.686, FDR-corrected p=0.0096, permutation p=0.0030). A monotonic enrichment gradient exists across developmental onset groups (ratio: 0.74 to 1.43). Disease system analysis reveals that congenital brain disorders show the strongest depletion (OR=0.640), suggesting that neurodevelopmental regulatory programs active in the earliest stages of embryogenesis remain systematically underannotated in current fetal regulatory databases.

---

## Repository Structure

    devvarp_project/
    01_download_data.py                  Downloads ClinVar, HPO, fetal atlas
    02_curate_variants.py                Filters ClinVar to pediatric non-coding variants
    03_fetal_regulatory_annotation.py    Overlaps variants with fetal chromatin atlas
    04_permutation_validation.py         Permutation test to validate main finding
    05_disease_system_analysis.py        Stratifies G1 by congenital disease system
    06_generate_figures.py               Generates all publication-grade figures
    data/processed/                      Curated variant datasets and results
    figures/                             All publication-grade figures (PNG, 300 DPI)

---

## Data Sources

| Dataset | Source | Version |
|---|---|---|
| ClinVar variant summary | NCBI ClinVar | May 2026 |
| Human fetal chromatin accessibility | Cao et al. 2020, Science (GSE149683) | v1 |
| HPO gene-phenotype annotations | Human Phenotype Ontology | v2026-02-16 |
| HPO ontology OBO | Human Phenotype Ontology | v2026-02-16 |

---

## Methods Summary

**Variant Curation:** 8,980,556 ClinVar variants were filtered to 65,124 high-confidence pathogenic and benign non-coding variants in pediatric-onset rare diseases using HPO age-of-onset annotations. Only variants with 2-star or higher review status were retained.

**Onset Stratification:** Variants were assigned to four developmental onset groups — Congenital (G1, n=26,770), Neonatal (G2, n=8,047), Infantile (G3, n=25,674), and Childhood (G4, n=4,633).

**Regulatory Annotation:** All variants were overlapped with 1,048,564 human fetal chromatin accessibility sites from the Cao et al. 2020 Science atlas using PyRanges genomic interval analysis.

**Statistical Analysis:** Fisher's exact tests with Benjamini-Hochberg FDR correction were applied across all onset groups. A permutation test with 1,000 iterations validated the primary finding.

**Disease System Analysis:** G1 Congenital variants were classified into four organ systems — Heart, Brain, Metabolic, Skeletal — using curated gene lists and tested independently.

---

## Results Summary

| Onset Group | n Pathogenic | Fetal Element % | Enrichment Ratio | FDR p |
|---|---|---|---|---|
| G1 Congenital | 479 | 15.9% | 0.736 | 0.0096 |
| G2 Neonatal | 154 | 16.2% | 0.825 | 0.307 |
| G3 Infantile | 591 | 18.9% | 0.881 | 0.189 |
| G4 Childhood | 102 | 28.4% | 1.431 | 0.089 |

Permutation test (G1 Congenital): observed ratio=0.736, permutation mean=1.000, p=0.003

---

## Installation

    git clone https://github.com/YLSVarshini/devvarp.git
    cd devvarp
    conda create -n devvarp python=3.10 -y
    conda activate devvarp
    conda install pandas numpy scipy scikit-learn matplotlib seaborn -y
    pip install pyranges statsmodels requests pyBigWig

## Usage

Run scripts in order:

    python 01_download_data.py
    python 02_curate_variants.py
    python 03_fetal_regulatory_annotation.py
    python 04_permutation_validation.py
    python 05_disease_system_analysis.py
    python 06_generate_figures.py

---

## Citation

Yedavalli, Lakshmi Sai Varshini. "Pathogenic Non-Coding Variants in Congenital-Onset Pediatric Rare Diseases Are Depleted from Annotated Fetal Regulatory Elements: Evidence for a Temporal Annotation Gap." Global Research Challenge, 2026.

---

## License

MIT License — free to use, modify, and distribute with attribution.
