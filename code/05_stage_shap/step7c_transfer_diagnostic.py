#!/usr/bin/env python3
"""
Step 7c - Diagnostic: is the partial non-replication a methodological artefact
of the pseudo-label recalibration, or a genuine biological cohort difference?

Tests THREE transfer methods for computing stage-stratified importance on
GSE39582, to see whether the gene-level pattern is robust to the transfer choice:

  M1 (direct linear attribution): apply the TCGA model's linear coefficients
      directly to z-scored GSE39582 panel; per-sample contribution of gene g is
      coef_g * z_g. Stage importance = mean|contribution|. NO recalibration, NO
      pseudo-labels - the cleanest transfer (linear SHAP == coef*value for a
      linear model, exactly).
  M2 (variance-weighted): importance = |coef_g| * SD_stage(z_g), capturing how
      much each gene actually varies within each stage subgroup.
  M3 (pseudo-label recalibration): the original Step-7b approach, for comparison.

If M1/M2 (clean, assumption-free) agree with each other but differ from M3, the
non-replication was partly an M3 artefact. If all three agree, the gene-level
non-replication is genuine biology.

License: MIT | Random seed: 42
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import pickle
from sklearn.preprocessing import StandardScaler

SEED = 42
np.random.seed(SEED)
PROC = Path("data/processed"); RES = Path("results/step7")
PANEL = ['CDH3','OTOP2','DHRS7C','AADACL2','KRT80','ETV4','ESM1','TMEFF2','LGI1','KHDRBS2']
STAGES = ["I","II","III","IV"]


def load_gse():
    expr = pd.read_csv(PROC/"GSE39582_expression_common.csv", index_col=0)
    meta = pd.read_csv("data/raw/GSE39582_metadata.csv").set_index("gsm")
    X = expr.loc[PANEL].T
    common = [s for s in X.index if s in meta.index]
    X = X.loc[common]
    stage = meta.loc[common,"tnm.stage"].map({1.:"I",2.:"II",3.:"III",4.:"IV"})
    keep = stage.notna().values
    return X.loc[keep], stage[keep]


def importance_by_stage(contrib, stage):
    """contrib: samples x genes signed contributions. importance=mean|contrib|."""
    out = {}
    for s in STAGES:
        m = (stage.values==s)
        out[s] = np.abs(contrib[m]).mean(axis=0)
    return out


def concord_matrix(imp):
    C = {}
    for i,a in enumerate(STAGES):
        for b in STAGES[i+1:]:
            rho,_ = spearmanr(imp[a], imp[b]); C[f"{a}_vs_{b}"]=float(rho)
    return C


def main():
    d = pickle.load(open("results/step6/final_model.pkl","rb"))
    coef = d["model"].coef_[0]
    Xpanel, stage = load_gse()
    Xz = StandardScaler().fit_transform(Xpanel.values)

    # M1: direct linear attribution (coef * z) - exact SHAP for linear model
    contrib_m1 = Xz * coef
    imp_m1 = importance_by_stage(contrib_m1, stage)

    # M2: variance-weighted |coef| * SD within stage
    imp_m2 = {}
    for s in STAGES:
        m=(stage.values==s)
        imp_m2[s] = np.abs(coef) * Xz[m].std(axis=0)

    # M3: load original pseudo-label result
    ext = json.load(open(RES/"external_validation_gse39582.json"))

    # TCGA reference
    tcga = pd.read_csv(RES/"shap_importance_by_stage.csv").set_index("gene").loc[PANEL]

    # --- compare concordance matrices ---
    print("=== Stage I-vs-IV concordance by transfer method ===", flush=True)
    print(f"  TCGA (reference):        {spearmanr(tcga['stage_I'],tcga['stage_IV'])[0]:.3f}", flush=True)
    print(f"  GSE39582 M1 (direct):    {concord_matrix(imp_m1)['I_vs_IV']:.3f}", flush=True)
    print(f"  GSE39582 M2 (var-wt):    {concord_matrix(imp_m2)['I_vs_IV']:.3f}", flush=True)
    print(f"  GSE39582 M3 (pseudo):    {ext['gse39582_stage_I_vs_IV_rho']:.3f}", flush=True)

    # --- CDH3 / ETV4 / ESM1 trajectories per method ---
    def ranks(imp, gene):
        gi = PANEL.index(gene)
        return {s:int(np.argsort(imp[s])[::-1].tolist().index(gi))+1 for s in STAGES}

    print("\n=== Key gene trajectories (rank; 1=most important) ===", flush=True)
    for gene in ["CDH3","ETV4","ESM1","KRT80"]:
        tcga_r = {s:int(np.argsort(tcga[f'stage_{s}'].values)[::-1].tolist().index(PANEL.index(gene)))+1 for s in STAGES}
        print(f"\n{gene}:", flush=True)
        print(f"  TCGA:        {tcga_r}", flush=True)
        print(f"  GSE39582 M1: {ranks(imp_m1, gene)}", flush=True)
        print(f"  GSE39582 M2: {ranks(imp_m2, gene)}", flush=True)

    # --- cross-cohort per-stage agreement for M1 (cleanest) ---
    print("\n=== Cross-cohort per-stage concordance (TCGA vs GSE39582 M1) ===", flush=True)
    cross={}
    for s in STAGES:
        rho,p = spearmanr(tcga[f"stage_{s}"].values, imp_m1[s])
        cross[s]={"rho":float(rho),"p":float(p)}
        print(f"  Stage {s}: rho={rho:.3f} p={p:.3f}", flush=True)

    # --- agreement BETWEEN methods (are M1/M2 consistent?) ---
    print("\n=== Agreement between transfer methods (mean per-stage rho) ===", flush=True)
    m1m2 = np.mean([spearmanr(imp_m1[s], imp_m2[s])[0] for s in STAGES])
    print(f"  M1 vs M2: {m1m2:.3f}", flush=True)

    result = {
        "tcga_I_vs_IV": float(spearmanr(tcga['stage_I'],tcga['stage_IV'])[0]),
        "gse_I_vs_IV_M1_direct": concord_matrix(imp_m1)['I_vs_IV'],
        "gse_I_vs_IV_M2_varwt": concord_matrix(imp_m2)['I_vs_IV'],
        "gse_I_vs_IV_M3_pseudo": ext['gse39582_stage_I_vs_IV_rho'],
        "cross_cohort_M1": cross,
        "M1_vs_M2_agreement": float(m1m2),
        "cdh3_M1": ranks(imp_m1,"CDH3"),
        "etv4_M1": ranks(imp_m1,"ETV4"),
        "esm1_M1": ranks(imp_m1,"ESM1"),
    }
    json.dump(result, open(RES/"transfer_diagnostic.json","w"), indent=2)
    # save M1 importances for potential figure
    pd.DataFrame({"gene":PANEL, **{f"stage_{s}":imp_m1[s] for s in STAGES}}).to_csv(
        RES/"gse39582_M1_importance.csv", index=False)
    print("\nDONE", flush=True)


if __name__ == "__main__":
    main()
