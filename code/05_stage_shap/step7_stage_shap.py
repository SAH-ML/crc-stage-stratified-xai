#!/usr/bin/env python3
"""
Step 7 - Stage-stratified SHAP attribution (the central novel analysis).

The Step-6 winning model (Elastic Net logistic regression on the 10-gene stable
panel, tumour-vs-normal) is explained with SHAP, computed:
  (a) GLOBALLY (all tumours pooled, vs the normal reference) -> stage-agnostic
      importance ranking, replicating standard CRC-SHAP practice; and
  (b) SEPARATELY WITHIN each AJCC stage subgroup (I, II, III, IV) using the SAME
      trained model and SAME background -> four stage-resolved importance
      rankings (the novel contribution).

The model is NOT retrained per stage; SHAP is a post-hoc stratification of one
full-cohort model (this is what makes the small Stage-I subgroup acceptable -
no per-subgroup training occurs).

Robustness / reviewer-proofing:
  * Exact explainer for the linear model (shap.LinearExplainer).
  * Bootstrap stability of each stage's top-gene ranking (1,000 iters):
    a stage ranking is "stable" if its top-k Spearman rho > 0.7 across resamples.
  * Concordance between global and per-stage rankings, and between stage pairs,
    via Spearman correlation with permutation-test p-values.
  * Pre-registered fallback (Methods 2.5.4): if any stage fails the stability
    bar, results are additionally reported at the early (I+II) vs late (III+IV)
    resolution. BOTH outcomes are reported; the decision is recorded.

License: MIT   |   Random seed: 42
"""
from pathlib import Path
import json, pickle, warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import shap

warnings.filterwarnings("ignore")
SEED = 42
rng = np.random.RandomState(SEED)
np.random.seed(SEED)

PROC = Path("data/processed")
RES = Path("results/step7")
RES.mkdir(parents=True, exist_ok=True)

N_BOOT = 1000
STABILITY_RHO = 0.7
TOPK = 10                # full panel is 10 genes


def mean_abs_shap(sv):
    """Mean absolute SHAP per feature -> importance vector."""
    return np.abs(sv).mean(axis=0)


def ranking(importance, genes):
    order = np.argsort(importance)[::-1]
    return [genes[i] for i in order], importance[order]


def bootstrap_stability(explainer, Xsub, genes, n_boot=N_BOOT):
    """Resample the subgroup's samples; recompute mean|SHAP| ranking each time;
    measure Spearman rho of the top-k importance vector vs the full-subgroup
    ranking. Returns mean rho and fraction of resamples with rho>threshold."""
    base_imp = mean_abs_shap(explainer.shap_values(Xsub))
    base_order = np.argsort(base_imp)[::-1]
    n = Xsub.shape[0]
    rhos = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        imp_b = mean_abs_shap(explainer.shap_values(Xsub[idx]))
        # Spearman between base and bootstrap importance (all genes)
        rho, _ = spearmanr(base_imp, imp_b)
        rhos.append(rho)
    rhos = np.array(rhos)
    return {
        "mean_rho": float(np.nanmean(rhos)),
        "sd_rho": float(np.nanstd(rhos)),
        "frac_rho_gt_thr": float(np.mean(rhos > STABILITY_RHO)),
        "stable": bool(np.nanmean(rhos) > STABILITY_RHO),
    }


def permutation_concordance(imp_a, imp_b, n_perm=1000):
    """Spearman rho between two importance vectors + permutation p-value."""
    rho_obs, _ = spearmanr(imp_a, imp_b)
    null = []
    b = imp_b.copy()
    for _ in range(n_perm):
        rng.shuffle(b)
        r, _ = spearmanr(imp_a, b)
        null.append(r)
    null = np.array(null)
    p = (np.sum(np.abs(null) >= abs(rho_obs)) + 1) / (n_perm + 1)
    return float(rho_obs), float(p)


def main():
    d = pickle.load(open("results/step6/final_model.pkl", "rb"))
    model, scaler, genes, idx = d["model"], d["scaler"], d["genes"], d["gene_indices"]
    print(f"Model: {d['winner']} | panel: {genes}", flush=True)

    expr = pd.read_csv(PROC / "TCGA_expression_common.csv", index_col=0)
    labels = pd.read_csv(PROC / "TCGA_labels.csv", index_col=0).loc[expr.columns]
    X_all = scaler.transform(expr.T.values[:, idx])     # scaled panel matrix
    tn = labels["tumor_vs_normal"].values
    stage = labels["stage"].values

    # background = full training distribution (consistent across all subgroups)
    explainer = shap.LinearExplainer(model, X_all)

    # ---- GLOBAL (all tumours) ----
    tumor_mask = tn == "tumor"
    Xt = X_all[tumor_mask]
    global_imp = mean_abs_shap(explainer.shap_values(Xt))
    g_order, g_vals = ranking(global_imp, genes)
    print(f"\nGlobal importance ranking: {g_order}", flush=True)

    # ---- PER-STAGE ----
    stage_imp = {}
    stage_stability = {}
    stage_counts = {}
    for stg in ["I", "II", "III", "IV"]:
        m = (stage == stg) & tumor_mask
        Xs = X_all[m]
        stage_counts[stg] = int(m.sum())
        imp = mean_abs_shap(explainer.shap_values(Xs))
        stage_imp[stg] = imp
        stab = bootstrap_stability(explainer, Xs, genes)
        stage_stability[stg] = stab
        order, _ = ranking(imp, genes)
        print(f"\nStage {stg} (n={m.sum()}): top genes {order[:5]}", flush=True)
        print(f"  stability mean_rho={stab['mean_rho']:.3f} stable={stab['stable']}",
              flush=True)

    # ---- concordance: global vs each stage, and stage pairs ----
    concord = {}
    for stg in ["I", "II", "III", "IV"]:
        rho, p = permutation_concordance(global_imp, stage_imp[stg])
        concord[f"global_vs_{stg}"] = {"rho": rho, "p": p}
    for a, b in [("I","II"),("I","III"),("I","IV"),
                 ("II","III"),("II","IV"),("III","IV")]:
        rho, p = permutation_concordance(stage_imp[a], stage_imp[b])
        concord[f"{a}_vs_{b}"] = {"rho": rho, "p": p}

    # ---- fallback decision (Methods 2.5.4) ----
    all_stable = all(stage_stability[s]["stable"] for s in ["I","II","III","IV"])
    min_n = min(stage_counts.values())
    decision = {
        "all_stages_meet_n30": bool(min_n >= 30),
        "all_stages_stable": bool(all_stable),
        "four_stage_primary": bool(all_stable and min_n >= 30),
        "fallback_triggered": bool(not (all_stable and min_n >= 30)),
        "stage_counts": stage_counts,
    }

    # ---- if fallback: also compute early (I+II) vs late (III+IV) ----
    early_late = {}
    for grp, stgs in [("early", ["I","II"]), ("late", ["III","IV"])]:
        m = np.isin(stage, stgs) & tumor_mask
        imp = mean_abs_shap(explainer.shap_values(X_all[m]))
        early_late[grp] = {"importance": imp.tolist(),
                           "ranking": ranking(imp, genes)[0],
                           "n": int(m.sum())}
    rho_el, p_el = permutation_concordance(
        np.array(early_late["early"]["importance"]),
        np.array(early_late["late"]["importance"]))
    early_late["early_vs_late_concordance"] = {"rho": rho_el, "p": p_el}

    # ---- save everything ----
    imp_df = pd.DataFrame({"gene": genes, "global": global_imp})
    for stg in ["I","II","III","IV"]:
        imp_df[f"stage_{stg}"] = stage_imp[stg]
    imp_df.to_csv(RES / "shap_importance_by_stage.csv", index=False)

    json.dump({k: {kk: (vv.tolist() if hasattr(vv,'tolist') else vv)
                   for kk,vv in v.items()} if isinstance(v,dict) else v
               for k,v in stage_stability.items()},
              open(RES / "stage_stability.json","w"), indent=2)
    json.dump(concord, open(RES / "concordance.json","w"), indent=2)
    json.dump(decision, open(RES / "fallback_decision.json","w"), indent=2)
    json.dump(early_late, open(RES / "early_late.json","w"), indent=2)

    # save raw shap arrays for figures
    shap_arrays = {"global": explainer.shap_values(Xt).tolist()}
    for stg in ["I","II","III","IV"]:
        m = (stage == stg) & tumor_mask
        shap_arrays[stg] = explainer.shap_values(X_all[m]).tolist()
    json.dump({"genes": genes, "shap": shap_arrays},
              open(RES / "shap_values.json","w"))

    print("\n=== DECISION ===", flush=True)
    print(json.dumps(decision, indent=2), flush=True)
    print("\n=== CONCORDANCE (global vs stage) ===", flush=True)
    for stg in ["I","II","III","IV"]:
        c = concord[f"global_vs_{stg}"]
        print(f"  global vs {stg}: rho={c['rho']:.3f} p={c['p']:.3f}", flush=True)
    print("\n=== STAGE-PAIR CONCORDANCE ===", flush=True)
    for a,b in [("I","IV"),("II","III")]:
        c = concord[f"{a}_vs_{b}"]
        print(f"  {a} vs {b}: rho={c['rho']:.3f} p={c['p']:.3f}", flush=True)
    print(f"\nArtefacts saved to {RES}/", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
