#!/usr/bin/env python3
"""
Step 9 - Survival benchmark on the 10-gene SHAP panel (TCGA-COAD overall survival).

Benchmarks three survival architectures (Methods 2.7.1) on the SHAP-prioritised
panel using TCGA-COAD overall survival (OS):
  * Cox Proportional Hazards (lifelines)            - linear baseline
  * Random Survival Forest (scikit-survival)        - non-parametric ensemble
  * DeepSurv (pycox CoxPH neural net)               - deep-learning comparator

EVENTS-PER-VARIABLE (EPV) HONESTY (Methods 2.7.1, Table 6 Issue 3, pre-registered)
--------------------------------------------------------------------------------
The 10-gene panel needs EPV>=10 (>=100 events) for a confirmatory model.
Actual TCGA OS events: 69 overall (EPV=6.9); per stage only 3/22/21/19
(EPV 0.3-2.2). THEREFORE:
  - the overall/pooled TCGA model is reported with an explicit EPV=6.9 caveat;
  - EVERY per-stage survival result is labelled EXPLORATORY, never confirmatory.
The EPV table is written out first; nothing is hidden.

Evaluation: Harrell's C-index + time-dependent AUC (3y, 5y) under 5-fold CV,
1,000-iteration percentile bootstrap 95% CIs, and a bootstrap comparison of
C-index between models. Leakage control: standardisation fit on each training
fold only.

License: MIT   |   Random seed: 42
"""
from pathlib import Path
import json, warnings, gc
import numpy as np
import pandas as pd

from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import cumulative_dynamic_auc
from sksurv.util import Surv

warnings.filterwarnings("ignore")
SEED = 42
np.random.seed(SEED)
rng = np.random.RandomState(SEED)

PROC = Path("data/processed")
RES = Path("results/step9")
RES.mkdir(parents=True, exist_ok=True)

PANEL = ['CDH3','OTOP2','DHRS7C','AADACL2','KRT80','ETV4','ESM1','TMEFF2','LGI1','KHDRBS2']
N_BOOT = 1000
N_FOLDS = 5
HORIZONS = {"3yr": 3*365.25, "5yr": 5*365.25}
EPV_THRESHOLD = 10
N_FEATURES = len(PANEL)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_tcga_survival():
    expr = pd.read_csv(PROC/"TCGA_expression_common.csv", index_col=0)
    surv = pd.read_csv(PROC/"TCGA_survival.csv", index_col=0)
    labels = pd.read_csv(PROC/"TCGA_labels.csv", index_col=0)
    tumor = labels[labels["tumor_vs_normal"]=="tumor"]
    common = [s for s in tumor.index if s in surv.index and s in expr.columns]
    X = expr.loc[PANEL, common].T                      # samples x genes
    y_time = surv.loc[common, "OS.time"].values.astype(float)
    y_event = surv.loc[common, "OS"].values.astype(int)
    stage = tumor.loc[common, "stage"].values
    # drop any rows with missing time
    ok = ~np.isnan(y_time)
    return X.loc[ok], y_time[ok], y_event[ok], stage[ok]


# ---------------------------------------------------------------------------
# EPV table (written first)
# ---------------------------------------------------------------------------
def write_epv_table(y_event, stage):
    rows = []
    for stg in ["I","II","III","IV"]:
        m = (stage==stg)
        ev = int(y_event[m].sum())
        rows.append({"subgroup": f"Stage {stg}", "n": int(m.sum()), "events": ev,
                     "EPV": round(ev/N_FEATURES,2),
                     "meets_EPV10": ev/N_FEATURES >= EPV_THRESHOLD,
                     "designation": "exploratory"})
    ev_all = int(y_event.sum())
    rows.append({"subgroup":"Overall (pooled)","n":int(len(y_event)),"events":ev_all,
                 "EPV":round(ev_all/N_FEATURES,2),
                 "meets_EPV10": ev_all/N_FEATURES >= EPV_THRESHOLD,
                 "designation": "primary (with EPV caveat)" if ev_all/N_FEATURES < EPV_THRESHOLD
                                else "confirmatory"})
    df = pd.DataFrame(rows)
    df.to_csv(RES/"epv_table.csv", index=False)
    print("=== EPV table (TCGA OS) ===", flush=True)
    print(df.to_string(index=False), flush=True)
    return df


# ---------------------------------------------------------------------------
# Survival models - fit + risk-score prediction
# ---------------------------------------------------------------------------
def fit_predict_cox(Xtr, ttr, etr, Xte):
    dftr = Xtr.copy(); dftr["T"]=ttr; dftr["E"]=etr
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(dftr, "T", "E")
    # higher partial hazard = higher risk
    risk = cph.predict_partial_hazard(Xte).values.ravel()
    return risk

def fit_predict_rsf(Xtr, ttr, etr, Xte):
    ytr = Surv.from_arrays(event=etr.astype(bool), time=ttr)
    rsf = RandomSurvivalForest(n_estimators=200, min_samples_leaf=10,
                               max_features="sqrt", random_state=SEED, n_jobs=-1)
    rsf.fit(Xtr.values, ytr)
    # risk score = predicted cumulative hazard (higher = higher risk)
    risk = rsf.predict(Xte.values)
    return risk

def fit_predict_deepsurv(Xtr, ttr, etr, Xte):
    import torch
    import torchtuples as tt
    from pycox.models import CoxPH as DeepSurv
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    xtr = Xtr.values.astype("float32")
    xte = Xte.values.astype("float32")
    ytr = (ttr.astype("float32"), etr.astype("float32"))
    net = tt.practical.MLPVanilla(
        in_features=xtr.shape[1], num_nodes=[16, 16], out_features=1,
        batch_norm=True, dropout=0.2, output_bias=False)
    model = DeepSurv(net, tt.optim.Adam(0.01))
    model.fit(xtr, ytr, batch_size=64, epochs=100, verbose=False)
    # risk = model output (log hazard); higher = higher risk
    risk = model.predict(xte).ravel()
    return risk


MODELS = {"CoxPH": fit_predict_cox, "RSF": fit_predict_rsf, "DeepSurv": fit_predict_deepsurv}


# ---------------------------------------------------------------------------
# Cross-validated evaluation
# ---------------------------------------------------------------------------
def cv_evaluate(X, t, e):
    """5-fold CV; returns pooled out-of-fold risk scores per model + fold c-indices."""
    kf = KFold(N_FOLDS, shuffle=True, random_state=SEED)
    oof_risk = {m: np.full(len(t), np.nan) for m in MODELS}
    fold_cindex = {m: [] for m in MODELS}
    for fi,(tri,tei) in enumerate(kf.split(X),1):
        sc = StandardScaler().fit(X.iloc[tri].values)
        Xtr = pd.DataFrame(sc.transform(X.iloc[tri].values), columns=X.columns, index=X.index[tri])
        Xte = pd.DataFrame(sc.transform(X.iloc[tei].values), columns=X.columns, index=X.index[tei])
        for m, fn in MODELS.items():
            try:
                risk = fn(Xtr, t[tri], e[tri], Xte)
                oof_risk[m][tei] = risk
                ci = concordance_index(t[tei], -risk, e[tei])  # -risk: higher risk -> shorter time
                fold_cindex[m].append(ci)
                print(f"  fold {fi} {m:9s} C-index={ci:.3f}", flush=True)
            except Exception as ex:
                print(f"  fold {fi} {m:9s} FAILED: {ex}", flush=True)
            gc.collect()
    return oof_risk, fold_cindex


def bootstrap_cindex_ci(t, e, risk, n_boot=N_BOOT):
    valid = ~np.isnan(risk)
    t, e, risk = t[valid], e[valid], risk[valid]
    point = concordance_index(t, -risk, e)
    boots = []
    n = len(t)
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        if e[idx].sum() < 2:
            continue
        try:
            boots.append(concordance_index(t[idx], -risk[idx], e[idx]))
        except Exception:
            pass
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return point, lo, hi, len(boots)


def time_dependent_auc(t, e, risk, train_t, train_e):
    """sksurv cumulative_dynamic_auc at 3y, 5y horizons."""
    valid = ~np.isnan(risk)
    t, e, risk = t[valid], e[valid], risk[valid]
    surv_tr = Surv.from_arrays(event=train_e.astype(bool), time=train_t)
    surv_te = Surv.from_arrays(event=e.astype(bool), time=t)
    out = {}
    for name, hz in HORIZONS.items():
        # only evaluate at horizons within follow-up range
        if hz >= t.max() or hz <= t.min():
            out[name] = None
            continue
        try:
            auc, _ = cumulative_dynamic_auc(surv_tr, surv_te, risk, [hz])
            out[name] = float(auc[0])
        except Exception:
            out[name] = None
    return out


def main():
    X, t, e, stage = load_tcga_survival()
    print(f"TCGA survival cohort: {len(t)} tumours, {int(e.sum())} events\n", flush=True)
    epv = write_epv_table(e, stage)

    print("\n=== Overall (pooled) survival benchmark — PRIMARY (EPV caveat) ===", flush=True)
    oof_risk, fold_cindex = cv_evaluate(X, t, e)

    # aggregate metrics with bootstrap CIs
    results = []
    for m in MODELS:
        risk = oof_risk[m]
        if np.all(np.isnan(risk)):
            continue
        pt, lo, hi, nb = bootstrap_cindex_ci(t, e, risk)
        tdauc = time_dependent_auc(t, e, risk, t, e)
        results.append({
            "model": m,
            "cindex": round(pt,4), "cindex_lo": round(lo,4), "cindex_hi": round(hi,4),
            "cindex_ci": f"{pt:.3f} [{lo:.3f}-{hi:.3f}]",
            "fold_cindex_mean": round(np.mean(fold_cindex[m]),4),
            "fold_cindex_sd": round(np.std(fold_cindex[m]),4),
            "auc_3yr": tdauc.get("3yr"), "auc_5yr": tdauc.get("5yr"),
            "n_bootstrap": nb,
        })
    res_df = pd.DataFrame(results).sort_values("cindex", ascending=False)
    res_df.to_csv(RES/"survival_benchmark_overall.csv", index=False)
    print("\n=== Overall survival benchmark results ===", flush=True)
    print(res_df.to_string(index=False), flush=True)

    # save oof risk scores (for KM plot using best model) + median split
    best_model = res_df.iloc[0]["model"]
    best_risk = oof_risk[best_model]
    km = pd.DataFrame({"time": t, "event": e, "stage": stage,
                       "risk": best_risk,
                       "risk_group": np.where(best_risk > np.nanmedian(best_risk),
                                              "high", "low")})
    km.to_csv(RES/"risk_scores_overall.csv", index=False)

    # ---- exploratory per-stage C-index (clearly labelled) ----
    print("\n=== Per-stage C-index — EXPLORATORY ONLY (all EPV<10) ===", flush=True)
    perstage = []
    for stg in ["I","II","III","IV"]:
        m = (stage==stg)
        if m.sum() < 10 or e[m].sum() < 3:
            perstage.append({"stage":stg,"n":int(m.sum()),"events":int(e[m].sum()),
                             "cindex_best_model":None,"note":"too few events"})
            continue
        risk_s = best_risk[m]
        valid = ~np.isnan(risk_s)
        if valid.sum() < 5 or e[m][valid].sum() < 2:
            ci = None
        else:
            ci = concordance_index(t[m][valid], -risk_s[valid], e[m][valid])
        perstage.append({"stage":stg,"n":int(m.sum()),"events":int(e[m].sum()),
                         "cindex_best_model": round(ci,4) if ci else None,
                         "note":"EXPLORATORY (EPV<10)"})
    ps_df = pd.DataFrame(perstage)
    ps_df.to_csv(RES/"survival_perstage_exploratory.csv", index=False)
    print(ps_df.to_string(index=False), flush=True)

    summary = {
        "best_model_overall": best_model,
        "overall_epv": float(e.sum()/N_FEATURES),
        "overall_epv_meets_threshold": bool(e.sum()/N_FEATURES >= EPV_THRESHOLD),
        "interpretation": ("Overall pooled model reported with EPV caveat "
                           "(EPV=%.1f < 10); all per-stage survival is exploratory."
                           % (e.sum()/N_FEATURES)),
    }
    json.dump(summary, open(RES/"survival_summary.json","w"), indent=2)
    print(f"\nBest overall model: {best_model}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
