# GLP1R Pharmacogenetics Analysis

Code repository for:

**"Rare GLP1R variants selectively attenuate weight loss with semaglutide versus tirzepatide"**

*Pat. Pending, US Pat App. 19/002,020.*

---

## Overview

This repository contains the analysis code for a pharmacogenetic study of rare coding variants in *GLP1R* and *GIPR* and their effect on weight loss during GLP-1 receptor agonist treatment (semaglutide and tirzepatide) in the Helix Research Network (HRN).

Key analyses:
- Rare variant carrier identification from exome sequencing (Hail / Helix Research Network)
- Treatment cohort assembly from OMOP CDM EHR data
- 1:10 nearest-neighbour matching with Mahalanobis distance
- Weighted linear regression for the primary outcome (minimum % weight change at 6–12 months)
- Drug × carrier interaction test (semaglutide vs tirzepatide)
- Regenie v2.2.4 burden test sensitivity analysis (European-ancestry subset)
- Replication in the All of Us Research Program

---

## Repository Structure

```
glp1r-pharmacogenetics/
├── config.R                          # Shared parameters, paths, colours
├── scripts/
│   ├── 01_variant_classification.py  # Hail: extract GLP1R/GIPR carriers from HRN exome WGS
│   ├── 02_cohort_selection.R         # OMOP CDM: assemble sema/tirz treatment cohorts
│   ├── 03_outcome_definition.R       # Variant classification; primary outcome definition
│   ├── 04_matching.R                 # 1:10 nearest-neighbour matching; covariate balance
│   ├── 05_primary_analysis.R         # Weighted regression; drug × carrier interaction
│   ├── 06_subgroup_analyses.R        # Variant class and clinical subgroup analyses
│   ├── 07_regenie_sensitivity.R      # Prepare regenie inputs; parse output
│   ├── 07_regenie_sensitivity.sh     # Shell: run regenie v2.2.4 (HRN infrastructure)
│   ├── 08_allofus_replication.R      # Replication in All of Us (BigQuery / AoU WGS)
│   └── 09_figures.R                  # All manuscript figures and Table 1
├── demo/
│   ├── generate_demo_data.R          # Generate fully synthetic demo input data
│   └── run_demo.R                    # Run pipeline end-to-end on demo data
└── output/                           # Demo input data committed; generated results not tracked
```

---

## Data Access

This analysis was conducted on two controlled-access platforms:

**Helix Research Network (HRN)**
Scripts `01`–`07` run on HRN infrastructure. Paths and credentials marked `[HELIX INTERNAL]` require approved access to HRN genotype and EHR data. See the Methods section of the paper for data access details.

**All of Us Research Program**
Script `08` runs on the All of Us Researcher Workbench (Google Cloud / BigQuery). Items marked `[AOU INTERNAL]` require an approved All of Us data access tier. See [https://www.researchallofus.org](https://www.researchallofus.org) for access information.

---

## System Requirements

**Operating system:** Linux (analysis performed on Debian GNU/Linux 10); demo tested on macOS and Linux. Windows is not supported for scripts requiring Hail or bash (`01`, `07`), but the demo (scripts `03`–`06`, `09`) is platform-independent.

**Hardware:** No special hardware required for the demo. Full analysis (scripts `01`–`07`) was run on the Helix Research Network HPC cluster with access to exome WGS data.

**Software:**
- R ≥ 4.2.3 (tested on 4.2.3)
- Python 3.7.11 (script `01` only; requires Hail 0.2.x on Apache Spark)
- regenie v2.2.4 (script `07` only)

---

## Demo (No Data Access Required)

A fully synthetic demo dataset is included in the repository (`output/` directory) so the pipeline can be run immediately without access to HRN or All of Us data.

```r
# From the glp1r-pharmacogenetics/ directory:
source("demo/run_demo.R")
```

This runs scripts `03`–`06` and `09` on the pre-committed synthetic data. Scripts `01`, `02`, `07`, and `08` require infrastructure access and are excluded from the demo. To regenerate the synthetic input data from scratch (e.g. with a different random seed), run `source("demo/generate_demo_data.R")` before `run_demo.R`.

**Expected output:** Running the demo produces the following files in `output/figures/`:
- `fig1a_forest_glp1r.png` — Fig 1a: GLP1R forest plot (response ratio by variant class, overall and by drug)
- `fig1b_forest_gipr.png` — Fig 1b: GIPR forest plot (response ratio by variant class, overall and by drug)
- `fig1c_drug_comparison.png` — Fig 1c: GLP1R violin/box plots of % weight change for carriers vs controls, stratified by semaglutide vs tirzepatide
- `fig1d_slopegraph_glp1r.png` — Fig 1d: GLP1R individual carrier weight trajectories by semaglutide vs tirzepatide
- `ed_fig_subgroup_glp1r.png` — Extended Data: GLP1R variant class subgroups and treatment-naive sensitivity analysis (forest plot)
- `ed_fig_subgroup_gipr.png` — Extended Data: GIPR variant class subgroups (forest plot)
- `supp_fig_gipr_violin.png` — Supplementary: GIPR violin/box plots by semaglutide vs tirzepatide (mirrors Fig 1c for GIPR)
- `supp_fig_gipr_slopegraph.png` — Supplementary: GIPR individual carrier weight trajectories (mirrors Fig 1d for GIPR)
- `table1_glp1r.csv` — Supplementary table: GLP1R cohort baseline characteristics (pre- and post-matching with standardised mean differences); the post-matching portion corresponds to manuscript Table 1
- `table1_gipr.csv` — Supplementary table: GIPR cohort baseline characteristics (pre- and post-matching with standardised mean differences); the post-matching portion corresponds to manuscript Table 1

**Expected runtime:** < 10 minutes on a standard laptop (synthetic dataset of ~15,600 individuals).

**Required R packages:**
```r
install.packages(c("tidyverse", "MatchIt", "sandwich", "lmtest", "RNOmni", "patchwork"))
```

**Typical install time:** ~5 minutes on a standard laptop (depends on internet speed and whether packages are already cached).

---

## Analysis Pipeline

### Step 1 — Variant extraction (`01_variant_classification.py`)
Extracts rare coding variant carriers for *GLP1R* and *GIPR* from pre-filtered HRN exome sequencing data using Hail 0.2.x on Apache Spark. Variants are filtered to MAF < 0.1% in gnomAD v4.1. Outputs one TSV per gene (one row per carrier–variant pair).

### Step 2 — Cohort selection (`02_cohort_selection.R`)
Assembles semaglutide and tirzepatide treatment cohorts from the HRN OMOP CDM. Applies inclusion/exclusion criteria, computes baseline BMI and follow-up duration, and builds longitudinal weight trajectories.

### Step 3 — Outcome definition (`03_outcome_definition.R`)
Classifies GLP1R and GIPR rare variants into a five-tier hierarchy (pLoF > damaging missense > rs146868158 > uncertain missense > benign missense) using REVEL and AlphaMissense scores. Defines the primary outcome: minimum % weight change from baseline in the 6–12 month window.

### Step 4 — Matching (`04_matching.R`)
Performs 1:10 nearest-neighbour matching using Mahalanobis distance (age, BMI, follow-up) with exact matching on drug, dose category, prior GLP-1 RA use, genetic similarity group, sex, and T2D status, plus BMI (±3 kg/m²) and follow-up (±30 day) calipers.

### Step 5 — Primary analysis (`05_primary_analysis.R`)
Weighted linear regression (ATT weights; HC1 SEs) for the primary outcome. Drug-stratified analyses and drug × carrier interaction test.

### Step 6 — Subgroup analyses (`06_subgroup_analyses.R`)
Variant class subgroups (pLoF+damaging, rs146868158, uncertain missense), clinical subgroup interactions (sex, T2D), and treatment-naive sensitivity analysis.

### Step 7 — Regenie sensitivity (`07_regenie_sensitivity.R` / `.sh`)
Gene-burden and SKAT-O tests in European-ancestry individuals using regenie v2.2.4. Three burden masks mirror the primary analysis variant groupings. Outcome is rank-normalized (RINT).

### Step 8 — All of Us replication (`08_allofus_replication.R`)
Independent replication in the All of Us Research Program using WGS data and BigQuery OMOP CDM. Same analytical specification as the primary HRN analysis.

### Step 9 — Figures (`09_figures.R`)
Generates all manuscript figures (forest plots, violin/box plots, covariate balance love plot, variant annotation summary, regenie results, AoU replication) and Table 1. Outputs PDF and PNG for each figure.

---

## Variant Classification Thresholds

| Category | Criterion |
|---|---|
| pLoF | LOFTEE HC, or stop_gained / frameshift / splice donor / splice acceptor |
| rs146868158 | chr6:39085942:C>T (GRCh38); queried separately (MAF slightly exceeds 0.1% threshold in Finnish; functionally validated) |
| Damaging missense | REVEL ≥ 0.75 **or** AlphaMissense ≥ 0.564 |
| Uncertain missense | (REVEL ≥ 0.25 **or** AlphaMissense ≥ 0.340) and not damaging |
| Benign missense | AlphaMissense < 0.34 **and** REVEL < 0.25 (scored but below both uncertain thresholds) |

Thresholds are defined in `config.R` (`PARAMS` list) and are not modified in any downstream script.

---

## Matching Specification

| Parameter | Value |
|---|---|
| Method | Nearest-neighbour without replacement (MatchIt) |
| Distance | Mahalanobis (age, BMI, follow-up duration) |
| Exact match | drug, dose category, prior sema/tirz, prior other GLP-1 RA, ancestry group, sex, T2D |
| Calipers | BMI ±3 kg/m², follow-up ±30 days (raw, not standardised) |
| Ratio | 1 carrier : 10 controls |
| Weights | Carrier = 1; each control = 1 / n controls in matched set |

---

## Citation

If you use this code, please cite the associated paper (citation to be added upon publication).

---

## License

CC-BY-NC-SA 4.0 — see `LICENSE`.

---

## Contact

For questions about data access, please refer to the Methods section of the paper.
For questions about the code, please open an issue in this repository.
