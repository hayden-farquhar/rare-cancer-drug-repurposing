"""Drug-specificity (prognostic-axis) conditioning for the three headline drugs.

Extension analysis (run after the main pipeline notebook). Tests whether each
headline drug-survival association is drug-specific or merely a re-expression of
the dominant correlated prognostic axis of the imputed-sensitivity matrix.

The 278 imputed drug-sensitivity profiles are strongly correlated (mean pairwise
r=0.78), so the first principal component (PC1) of the within-cancer standardised
matrix is a data-driven proxy for that shared axis. We (1) test PC1's own overall
-survival association, and (2) refit each headline drug's Cox model adjusting for
PC1. If the drug term collapses toward HR=1 / loses significance, the association
is largely the generic axis; if it persists, it is drug-specific. We also run a
drug-identity null: among all 278 drugs, how many reach a univariate OS
association at least as strong (by p-value) as the headline drug?

Inputs (in ../results, relative to this script):
    imputed_drug_responses.parquet   imputed per-patient sensitivity, 401 x 278
    tcga_clinical.parquet            os_days, os_event per patient
Output:
    ../results/specificity_null.csv  one row per headline drug (backs Table S6)

Run:  python scripts/drug_specificity_null.py
"""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter

RESULTS = Path(__file__).resolve().parent.parent / "results"
OUT = RESULTS / "specificity_null.csv"

imp = pd.read_parquet(RESULTS / "imputed_drug_responses.parquet")
clin = pd.read_parquet(RESULTS / "tcga_clinical.parquet")

drug_cols = [c for c in imp.columns if c != "project_id"]
HEADLINE = {"TCGA-ACC": "Methotrexate", "TCGA-MESO": "Remodelin", "TCGA-UVM": "Serdemetan"}


def cox1(df, cov):
    cph = CoxPHFitter().fit(df[["os_days", "os_event", cov]], "os_days", "os_event")
    return cph.summary.loc[cov, ["exp(coef)", "p"]]


rows = []
for cancer, drug in HEADLINE.items():
    pts = imp.index[imp["project_id"] == cancer]
    sub = imp.loc[pts, drug_cols].join(clin[["os_days", "os_event"]], how="inner")
    sub = sub[sub["os_days"].notna() & (sub["os_days"] > 0) & sub["os_event"].notna()]
    surv = sub[["os_days", "os_event"]].astype(float)
    X = StandardScaler().fit_transform(sub[drug_cols].values)
    Xdf = pd.DataFrame(X, columns=drug_cols, index=sub.index)

    pca = PCA(n_components=1).fit(X)
    pc1 = StandardScaler().fit_transform(pca.transform(X))[:, 0]
    var1 = pca.explained_variance_ratio_[0]
    drug_z = Xdf[drug].values
    # align PC1 sign with the headline drug for interpretable HRs
    if np.corrcoef(pc1, drug_z)[0, 1] < 0:
        pc1 = -pc1
    r_drug_pc1 = np.corrcoef(pc1, drug_z)[0, 1]

    a = cox1(surv.assign(drug=drug_z), "drug")
    b = cox1(surv.assign(pc1=pc1), "pc1")
    joint = surv.assign(drug=drug_z, pc1=pc1)
    cj = CoxPHFitter().fit(joint[["os_days", "os_event", "drug", "pc1"]],
                           "os_days", "os_event").summary

    # drug-identity null: univariate per-SD Cox for every drug
    pvals = {}
    for d in drug_cols:
        try:
            pvals[d] = cox1(surv.assign(x=Xdf[d].values), "x")["p"]
        except Exception:
            pvals[d] = np.nan
    pser = pd.Series(pvals).dropna()
    head_p = pser[drug]
    n_at_least = int((pser <= head_p).sum())

    rows.append(dict(
        cancer=cancer.replace("TCGA-", ""), drug=drug, n=len(sub),
        events=int(surv["os_event"].sum()),
        HR_drug_uni=round(a["exp(coef)"], 3), p_drug_uni=a["p"],
        HR_pc1_uni=round(b["exp(coef)"], 3), p_pc1_uni=b["p"],
        HR_drug_adj=round(cj.loc["drug", "exp(coef)"], 3), p_drug_adj=cj.loc["drug", "p"],
        HR_pc1_adj=round(cj.loc["pc1", "exp(coef)"], 3), p_pc1_adj=cj.loc["pc1", "p"],
        pc1_var_explained=round(var1, 3), r_drug_pc1=round(r_drug_pc1, 3),
        n_drugs_ge=n_at_least, frac_drugs_ge=round(n_at_least / len(pser), 3),
    ))

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)
pd.set_option("display.width", 200, "display.max_columns", 30)
print(out.to_string(index=False))
print("\nsaved", OUT.name)
