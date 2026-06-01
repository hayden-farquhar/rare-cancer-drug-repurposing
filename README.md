# Pharmacogenomic Transfer Learning for Rare-Cancer Drug Repurposing

Replication code and results for:

> Farquhar H. Computational nomination of candidate repurposing drugs for five rare cancers via GDSC-to-TCGA pharmacogenomic transfer learning. 2026. Zenodo preprint. https://doi.org/10.5281/zenodo.19546768
>
> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19535800.svg)](https://doi.org/10.5281/zenodo.19535800)

## Overview

This repository contains the analysis pipeline for nominating drug-repurposing candidates across five TCGA rare cancer types (mesothelioma [MESO], adrenocortical carcinoma [ACC], uveal melanoma [UVM], cholangiocarcinoma [CHOL], thymoma [THYM]) using transfer learning from the GDSC cell-line pharmacogenomic dataset.

**Interpretation caveat (important).** The associations reported here are **prognostic, not predictive**: stratifying patients by imputed drug sensitivity and observing a survival difference shows that the expression programme the model reads is associated with prognosis, not that administering the drug would change outcome. The output is a prioritised, mechanistically interpretable **hypothesis set** for experimental follow-up — strongest for ACC, intermediate for MESO, most tentative for UVM — not evidence of clinical efficacy.

## Quick start

The core pipeline runs in a single self-contained Jupyter notebook on [Google Colab](https://colab.research.google.com/) (free tier, CPU-only, ~90 minutes):

1. Open `rare_cancer_pipeline.ipynb` in Google Colab
2. Runtime > Run all
3. Results are saved and downloaded automatically as a zip file

No local installation, GPU, or institutional compute is required. All data is downloaded automatically from public APIs during execution.

> **Scope of the notebook vs. extension analyses.** `rare_cancer_pipeline.ipynb` runs the base pipeline: data acquisition, ridge-regression imputation, per-drug survival association (log-rank + Cox, FDR-corrected), IDWAS, GSEA, SHAP, permutation testing, and split-half validation. Four robustness/validation analyses added during peer-review revision were run as **separate steps** after the base notebook and are documented under [Extension analyses](#extension-analyses) below. Their output tables are included in `results/`; the drug-specificity null is provided as a runnable script.

## Repository structure

```
├── rare_cancer_pipeline.ipynb   # Core analysis pipeline (Colab-ready)
├── scripts/
│   └── drug_specificity_null.py # PC1 prognostic-axis conditioning (extension; backs Table S6)
├── results/                     # Pre-computed results
│   ├── ridge_cv_scores.csv          # CV performance for all 286 drugs
│   ├── imputed_drug_responses.parquet   # Imputed IC50 for 401 patients x 278 drugs
│   ├── tcga_clinical.parquet        # Per-patient overall-survival time + event
│   ├── all_cancer_survival.csv      # Survival analysis across 5 cancer types
│   ├── meso_survival_results.csv    # Mesothelioma-specific survival results
│   ├── gsea_all_cancers.csv         # GSEA pathway enrichment
│   ├── idwas_all_cancers.csv.gz     # Gene-drug association results
│   ├── shap_biomarkers_corrected.csv    # SHAP feature importance (corrected direction)
│   ├── shap_biomarkers.csv          # SHAP (original; superseded — see note)
│   ├── firth_cox_headline.csv       # Firth-penalised per-SD HRs (extension; Table 1)
│   ├── splithalf_repeated.csv       # 100-split repeated reproducibility (extension)
│   ├── splithalf_r_values.csv       # Per-split correlations
│   ├── meso_perm_counts.csv         # Permutation null counts (MESO)
│   ├── validation_extension.csv     # Per-cancer permutation + split-half summary
│   ├── ctrp_bet_meso_crosscheck.csv # CTRPv2 BET-inhibitor cross-check (extension; Table S5)
│   ├── ctrp_bet_percell.csv         # CTRPv2 per-cell-line AUCs
│   └── specificity_null.csv         # Drug-specificity PC1 null (extension; Table S6)
├── figures/                     # All figures from the analysis
├── requirements.txt             # Python dependencies
├── LICENSE
└── README.md
```

> **SHAP correction.** `shap_biomarkers.csv` derived feature direction from the mean signed SHAP value, which is ~zero for every feature by construction and therefore non-informative. `shap_biomarkers_corrected.csv` assigns direction from the sign of the standardised ridge coefficient and is the version used in the manuscript (Figure 5, Table S4).

## Data sources

All data is publicly available. The core notebook downloads its inputs automatically.

| Dataset | Source | Access | Used in |
|---------|--------|--------|---------|
| GDSC2 drug sensitivity (IC50) | [cancerrxgene.org](https://www.cancerrxgene.org/) | Free download | Core pipeline |
| Cell-line gene expression | [DepMap portal](https://depmap.org/portal/) (release 24Q4) | Free API | Core pipeline |
| TCGA RNA-seq + clinical | [GDC API](https://portal.gdc.cancer.gov/) | Open access | Core pipeline |
| CTRPv2 drug sensitivity (AUC) | [CTD² / CTRPv2](https://ctd2-data.nci.nih.gov/Public/Broad/) | Open access | External validation (extension) |

## Methods summary (core pipeline)

1. **Data acquisition:** GDSC2 IC50 for 286 drugs; DepMap RNA-seq (log2 TPM+1) for matched cell lines; TCGA RNA-seq (log2 FPKM+1) and clinical data for 401 patients across 5 rare cancers
2. **Preprocessing:** gene harmonisation (12,512 common protein-coding genes), quantile normalisation (TCGA to GDSC reference), variance filtering (top 5,000 genes)
3. **Modelling:** ridge regression per drug (5-fold CV, alpha from 20 values); 278/286 drugs predictable at CV R² > 0.1 (171 at R² > 0.3, 237 at R² > 0.2 — the headline proportion is threshold-dependent)
4. **Imputation:** trained models applied to TCGA expression to predict drug sensitivity per patient
5. **Survival analysis:** log-rank + Cox PH per drug per cancer, FDR-corrected (nominal within-cancer Benjamini-Hochberg)
6. **Interpretation:** IDWAS pharmacogenomic associations, GSEA pathway enrichment, SHAP biomarkers
7. **Robustness:** label-permutation testing, split-half replication

## Extension analyses

Run after the base notebook during peer-review revision. Output tables are included in `results/`; for the analyses whose Colab source is not bundled here, the tables are provided as the reproducible record and the methods are specified in the manuscript.

- **Firth-penalised per-SD Cox** (`results/firth_cox_headline.csv`, manuscript Table 1). Unstandardised HRs in these small cohorts are inflated by near-complete separation; per-standard-deviation Firth-penalised estimation gives bounded effect sizes: MESO Remodelin 1.72 (1.32–2.22), ACC Methotrexate 0.40 (0.27–0.57), UVM Serdemetan 0.48 (0.32–0.69).
- **Repeated split-half resampling** (`results/splithalf_repeated.csv`, `splithalf_r_values.csv`). Reproducibility over 100 random partitions per cancer: mean between-half r = 0.24 in ACC (positive in 95% of splits), 0.19 in MESO (88%), 0.14 in UVM (78%).
- **CTRPv2 external validation** (`results/ctrp_bet_meso_crosscheck.csv`, `ctrp_bet_percell.csv`, manuscript Figure 7, Table S5). The leading mesothelioma BET/bromodomain theme did **not** replicate in the independent CTRPv2 cell-line screen: mesothelioma lines showed no preferential BET-class sensitivity (median AUC 0.84 vs 0.83; one-sided Mann-Whitney p = 0.76), even though the same screen recovered the canonical haematopoietic BET dependency (p = 9×10⁻⁵¹), confirming the test had power.
- **Drug-specificity (prognostic-axis) conditioning** (`scripts/drug_specificity_null.py` → `results/specificity_null.csv`, manuscript Table S6). Because the 278 imputed profiles are strongly correlated (mean pairwise r = 0.78), PC1 of the within-cancer standardised matrix is a data-driven proxy for the shared prognostic axis. Each headline association persists or strengthens after adjusting for PC1 (ACC 0.40→0.29, MESO 1.71→1.71, UVM 0.48→0.34) and sits in the top ~1–2% of all 278 drugs by univariate OS association. Because PC1 captures only about half the matrix variance (47–55%), this conditions out the single dominant correlated direction but not residual structure in higher principal components — so the signals are specific *beyond that dominant axis* rather than fully independent of all shared structure. This script is runnable: `python scripts/drug_specificity_null.py`.

## Key results

- **114 nominal FDR-significant drug-survival associations** across 3 cancer types (ACC: 59, MESO: 39, UVM: 16; none in CHOL or THYM). Because imputed profiles are strongly correlated, these counts index a **ranking**, not independent discoveries.
- Per-cancer aggregate signal exceeds its own correlation-preserving permutation null (empirical p = 0.01–0.02; the smallest p, in UVM, reflects that cancer's especially thin null rather than a stronger signal — UVM was the least reproducible on resampling).
- Top candidates are pathway-concordant: bromodomain/BET inhibitors in BAP1-mutant MESO, antiproliferatives (Methotrexate, Axitinib) in ACC, MDM2 inhibition (Serdemetan) in p53-intact UVM; SHAP recovers BCL2L1.
- Bounded Firth per-SD HRs and a drug-specificity null support the three headline drugs; the MESO BET theme did not externally replicate (CTRPv2 p = 0.76).
- Reproducibility is strongest in ACC, modest in MESO, weakest in UVM (see Extension analyses).

## Requirements

- Python >= 3.10
- Google Colab (free tier) or a local environment with the packages in `requirements.txt`
- Internet connection (for data download on first run)
- ~6 GB RAM, ~500 MB disk

## Citation

If you use this code or data, please cite:

```bibtex
@misc{farquhar2026rarecancer,
  title={Computational nomination of candidate repurposing drugs for five rare cancers via GDSC-to-TCGA pharmacogenomic transfer learning},
  author={Farquhar, Hayden},
  year={2026},
  note={Zenodo preprint. \url{https://doi.org/10.5281/zenodo.19546768}},
  doi={10.5281/zenodo.19535800}
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

The underlying datasets (GDSC, DepMap, TCGA, CTRPv2) are subject to their own data use agreements. See their respective websites for terms.
