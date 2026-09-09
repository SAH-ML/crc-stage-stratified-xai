#!/usr/bin/env python3
"""
Step 10 - Three-tier external validation (Methods 2.7.2), TRIPOD+AI-aligned.

TIER 1 (same-platform, RNA-seq) - GSE50760:
  Classification replication ONLY (cohort is all Stage IV; no per-stage/survival).
  The Step-6 tumour-vs-normal Elastic Net (10-gene panel) is applied to GSE50760
  primary-tumour vs normal samples via the signature-transfer protocol (per-gene
  z-score standardisation within GSE50760, then the trained model's linear
  decision function). Metric: AUC + AUPRC + bootstrap 95% CI.

TIER 2 (cross-platform, microarray) - GSE39582 + GSE17536:
  SURVIVAL validation via signature transfer. The 10-gene panel is z-scored within
  each cohort; a signature risk score is formed from the panel and its association
  with overall survival assessed by Harrell's C-index with bootstrap 95% CI.
  EPV is reported per cohort; per-stage survival is exploratory only.

All cohorts are external (never used in model development). Random seed 42.
Honesty: results reported regardless of direction; weak signals stated plainly.

License: MIT   |   Random seed: 42
"""
from pathlib import Path
import json, warnings
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score
from lifelines.utils import concordance_index
import pickle

warnings.filterwarnings("ignore")
SEED = 42
rng = np.random.RandomState(SEED)
np.random.seed(SEED)

PROC = Path("data/processed")
RES = Path("results/step10")
RES.mkdir(parents=True, exist_ok=True)
PANEL = ['CDH3','OTOP2','DHRS7C','AADACL2','KRT80','ETV4','ESM1','TMEFF2','LGI1','KHDRBS2']
N_BOOT = 1000
N_FEATURES = len(PANEL)


def boot_ci(metric_fn, *args, n_boot=N_BOOT, need_both_classes=False):
    """Generic percentile bootstrap CI."""
    base = metric_fn(*args)
    n = len(args[0])
    boots = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        resampled = [a[idx] for a in args]
        if need_both_classes and len(np.unique(resampled[0])) < 2:
            continue
        try:
            boots.append(metric_fn(*resampled))
        except Exception:
            pass
    if not boots:
        return base, np.nan, np.nan, 0
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return base, lo, hi, len(boots)


# ---------------------------------------------------------------------------
# TIER 1 - GSE50760 classification replication
# ---------------------------------------------------------------------------
def tier1():
    print("=== TIER 1: GSE50760 classification replication ===", flush=True)
    d = pickle.load(open("results/step6/final_model.pkl","rb"))
    model = d["model"]; coef = model.coef_[0]; intercept = model.intercept_[0]

    expr = pd.read_csv(PROC/"GSE50760_expression_common.csv", index_col=0)
    meta = pd.read_csv("data/raw/GSE50760_metadata.csv")
    tissue = meta.set_index("gsm")["tissue"]

    # primary tumour vs normal (exclude metastatic for the clean tumour-vs-normal test)
    samples = [s for s in expr.columns if s in tissue.index]
    label = tissue.loc[samples].map(
        lambda x: 1 if "primary colorectal cancer" in str(x)
        else (0 if "normal" in str(x) else np.nan))
    keep = label.notna().values
    samples = np.array(samples)[keep]; y = label[keep].values.astype(int)

    X = expr.loc[PANEL, samples].T.values
    Xz = StandardScaler().fit_transform(X)            # signature-transfer z-score
    score = Xz @ coef + intercept                     # linear decision function
    prob = 1/(1+np.exp(-score))

    auc, auc_lo, auc_hi, nb = boot_ci(roc_auc_score, y, prob, need_both_classes=True)
    auprc, ap_lo, ap_hi, _ = boot_ci(average_precision_score, y, prob, need_both_classes=True)
    print(f"  n={len(y)} ({y.sum()} tumour / {(y==0).sum()} normal)", flush=True)
    print(f"  AUC={auc:.3f} [{auc_lo:.3f}-{auc_hi:.3f}]", flush=True)
    print(f"  AUPRC={auprc:.3f} [{ap_lo:.3f}-{ap_hi:.3f}]", flush=True)
    out = {"tier":"1","cohort":"GSE50760","task":"primary tumour vs normal (RNA-seq)",
           "n":int(len(y)),"n_tumor":int(y.sum()),"n_normal":int((y==0).sum()),
           "auc":round(auc,4),"auc_ci":[round(auc_lo,4),round(auc_hi,4)],
           "auprc":round(auprc,4),"auprc_ci":[round(ap_lo,4),round(ap_hi,4)],
           "n_bootstrap":nb}
    return out


# ---------------------------------------------------------------------------
# TIER 2 - survival validation via signature transfer
# ---------------------------------------------------------------------------
def signature_risk_score(expr_cohort):
    """z-score panel within cohort, weight by TCGA model coefficients -> risk."""
    d = pickle.load(open("results/step6/final_model.pkl","rb"))
    coef = d["model"].coef_[0]
    X = expr_cohort.loc[PANEL].T.values
    Xz = StandardScaler().fit_transform(X)
    # risk score = linear combination (higher = more 'tumour-like'); test survival assoc
    return Xz @ coef, expr_cohort.columns.tolist()


def tier2_survival(cohort_name, expr_path, meta_path, time_col, event_col,
                   event_map=None, stage_col=None, stage_map=None, time_unit="months"):
    print(f"\n=== TIER 2: {cohort_name} survival validation ===", flush=True)
    expr = pd.read_csv(expr_path, index_col=0)
    meta = pd.read_csv(meta_path)
    idcol = "gsm" if "gsm" in meta.columns else meta.columns[0]
    meta = meta.set_index(idcol)

    risk, samples = signature_risk_score(expr)
    risk = pd.Series(risk, index=samples)
    common = [s for s in samples if s in meta.index]
    risk = risk.loc[common]

    t = pd.to_numeric(meta.loc[common, time_col], errors="coerce")
    if event_map:
        e = meta.loc[common, event_col].map(event_map)
    else:
        e = pd.to_numeric(meta.loc[common, event_col], errors="coerce")
    ok = t.notna() & e.notna()
    t, e, risk = t[ok].values.astype(float), e[ok].values.astype(int), risk[ok].values
    if time_unit == "months":
        t = t * 30.44

    events = int(e.sum())
    epv = events / N_FEATURES
    ci, lo, hi, nb = boot_ci(lambda tt,ee,rr: concordance_index(tt,-rr,ee),
                             t, e, risk)
    print(f"  n={len(t)}, events={events}, EPV={epv:.1f}", flush=True)
    print(f"  C-index={ci:.3f} [{lo:.3f}-{hi:.3f}]", flush=True)
    out = {"tier":"2","cohort":cohort_name,"task":"survival (signature transfer)",
           "n":int(len(t)),"events":events,"epv":round(epv,2),
           "epv_meets_threshold": bool(epv>=10),
           "cindex":round(ci,4),"cindex_ci":[round(lo,4),round(hi,4)],
           "n_bootstrap":nb}
    return out


def main():
    results = {}
    results["tier1_gse50760"] = tier1()

    results["tier2_gse39582"] = tier2_survival(
        "GSE39582", PROC/"GSE39582_expression_common.csv",
        "data/raw/GSE39582_metadata.csv",
        time_col="os.delay (months)", event_col="os.event", time_unit="months")

    results["tier2_gse17536"] = tier2_survival(
        "GSE17536", PROC/"GSE17536_expression_common.csv",
        "data/raw/GSE17536_metadata.csv",
        time_col="overall survival follow-up time",
        event_col="overall_event (death from any cause)",
        event_map={"death":1,"no death":0}, time_unit="months")

    json.dump(results, open(RES/"external_validation_summary.json","w"), indent=2)

    # tidy table
    rows = []
    t1 = results["tier1_gse50760"]
    rows.append({"tier":"1","cohort":"GSE50760","task":"classification",
                 "metric":"AUC","value":t1["auc"],
                 "ci":f"[{t1['auc_ci'][0]}-{t1['auc_ci'][1]}]","n":t1["n"]})
    for k in ["tier2_gse39582","tier2_gse17536"]:
        r = results[k]
        rows.append({"tier":"2","cohort":r["cohort"],"task":"survival",
                     "metric":"C-index","value":r["cindex"],
                     "ci":f"[{r['cindex_ci'][0]}-{r['cindex_ci'][1]}]","n":r["n"]})
    pd.DataFrame(rows).to_csv(RES/"external_validation_table.csv", index=False)

    print("\n=== EXTERNAL VALIDATION SUMMARY ===", flush=True)
    print(pd.DataFrame(rows).to_string(index=False), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
