# Pharmacogenomic Drug Repurposing for Rare Cancers

Replication code and results for:

> Farquhar H. Pharmacogenomic imputation identifies epigenetic drug candidates for mesothelioma and actionable therapies for four other rare cancers. *[Journal]*, 2026.

## Overview

This repository contains the complete analysis pipeline for identifying drug repurposing candidates across five TCGA rare cancer types (mesothelioma, adrenocortical carcinoma, uveal melanoma, cholangiocarcinoma, thymoma) using transfer learning from the GDSC cell line pharmacogenomic dataset.

## Quick start

The entire analysis runs in a single self-contained Jupyter notebook on [Google Colab](https://colab.research.google.com/) (free tier, CPU-only, ~90 minutes total):

1. Open `rare_cancer_pipeline.ipynb` in Google Colab
2. Runtime > Run all
3. Results are saved and downloaded automatically as a zip file

No local installation, GPU, or institutional compute is required. All data is downloaded automatically from public APIs during execution.

## Repository structure

```
├── rare_cancer_pipeline.ipynb   # Complete analysis pipeline (Colab-ready)
├── results/                     # Pre-computed results
│   ├── ridge_cv_scores.csv          # CV performance for all 286 drugs
│   ├── imputed_drug_responses.parquet   # Imputed IC50 for 401 patients x 278 drugs
│   ├── all_cancer_survival.csv      # Survival analysis across 5 cancer types
│   ├── meso_survival_results.csv    # Mesothelioma-specific survival results
│   ├── gsea_all_cancers.csv         # GSEA pathway enrichment
│   ├── idwas_all_cancers.csv        # Gene-drug association results
│   └── shap_biomarkers.csv          # SHAP feature importance
├── figures/                     # All figures from the analysis
│   ├── ridge_cv_performance.png
│   ├── imputed_heatmap.png
│   ├── meso_km_curves.png
│   ├── gsea_barplots.png
│   ├── shap_biomarkers.png
│   ├── meso_volcano.png
│   ├── all_cancer_km.png
│   ├── permutation_test.png
│   ├── split_half.png
│   ├── drug_drug_clustering.png
│   └── normalisation.png
├── requirements.txt             # Python dependencies
├── LICENSE
└── README.md
```

## Data sources

All data is publicly available and downloaded automatically by the notebook:

| Dataset | Source | Access |
|---------|--------|--------|
| GDSC2 drug sensitivity (IC50) | [cancerrxgene.org](https://www.cancerrxgene.org/) | Free download |
| Cell line gene expression | [DepMap portal](https://depmap.org/portal/) (release 24Q4) | Free API |
| TCGA RNA-seq + clinical | [GDC API](https://portal.gdc.cancer.gov/) | Open access |

## Methods summary

1. **Data acquisition:** GDSC2 IC50 values for 286 drugs across 969 cell lines; DepMap RNA-seq expression (log2 TPM+1) for 717 matched cell lines; TCGA RNA-seq (log2 FPKM+1) and clinical data for 401 patients across 5 rare cancers
2. **Preprocessing:** Gene harmonisation (12,512 common protein-coding genes), quantile normalisation (TCGA to GDSC reference), variance filtering (top 5,000 genes)
3. **Modelling:** Ridge regression per drug (5-fold CV, alpha selection from 20 values); 278/286 drugs predictable (CV R² > 0.1)
4. **Imputation:** Trained models applied to TCGA expression to predict drug sensitivity for all patients
5. **Survival analysis:** Log-rank test + Cox PH per drug per cancer type, FDR-corrected
6. **Validation:** IDWAS pharmacogenomic associations, GSEA pathway enrichment, SHAP biomarker analysis, permutation testing (100 permutations), split-half replication

## Key results

- **114 FDR-significant drug-survival associations** across 3 cancer types (ACC: 59, MESO: 39, UVM: 16)
- Top mesothelioma candidates: Remodelin, LCL161, OF-1, CPI-637 (epigenetic drugs clustering consistent with BAP1 loss biology)
- Top ACC candidates: Methotrexate, Axitinib (FDA-approved for RCC)
- Top UVM candidate: Serdemetan (MDM2 inhibitor; consistent with intact p53 in UVM)
- Permutation test confirms signal (p = 0.01); split-half validates replicability (r = 0.357, p = 9.1e-10)

## Requirements

- Python >= 3.10
- Google Colab (free tier) or local environment with packages in `requirements.txt`
- Internet connection (for data download on first run)
- ~6 GB RAM, ~500 MB disk

## Citation

If you use this code or data, please cite:

```bibtex
@article{farquhar2026rare,
  title={Pharmacogenomic imputation identifies epigenetic drug candidates for mesothelioma and actionable therapies for four other rare cancers},
  author={Farquhar, Hayden},
  journal={[Journal]},
  year={2026}
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

The underlying datasets (GDSC, DepMap, TCGA) are subject to their own data use agreements. See their respective websites for terms.
