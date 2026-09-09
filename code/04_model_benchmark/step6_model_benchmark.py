#!/usr/bin/env python3
"""
Step 6 - Six-model benchmark under leak-free nested cross-validation, with
         multi-metric evaluation, class-imbalance diagnosis, DeLong-test based
         model selection, and feature-selection stability.

Rationale for the strengthened design
--------------------------------------
The cohort is imbalanced (286 tumour vs 41 normal, ~7:1). AUC alone is
insensitive to imbalance, so model screening additionally reports metrics that
ARE imbalance-sensitive: AUPRC, Matthews correlation coefficient (MCC),
balanced accuracy, and per-class sensitivity/specificity. Model selection is
made principled: models are first tested for statistical equivalence in
discrimination with DeLong's test on pooled out-of-fold predictions; among
models that are statistically indistinguishable, the most parsimonious is
chosen (Occam's razor), with the choice justified explicitly rather than by
sort order. A class-weighted variant (class_weight='balanced') is benchmarked
alongside the unweighted models to test whether imbalance correction changes
anything; SMOTE is deliberately NOT used (unreliable in ~16k-dimensional
expression space and leakage-prone), consistent with a measure-first,
correct-only-if-needed, prefer-safest-correction policy.

Leakage control: feature selection + scaling fit on each outer training fold
only; held-out fold transformed without refitting.

License: MIT   |   Random seed: 42
"""
from pathlib import Path
import json, pickle, warnings, gc
import numpy as np
import pandas as pd
from itertools import combinations
from scipy.stats import loguniform, norm
from sklearn.model_selection import StratifiedKFold, ParameterSampler
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.feature_selection import f_classif
from sklearn.metrics import (roc_auc_score, accuracy_score, recall_score,
                             f1_score, confusion_matrix, average_precision_score,
                             matthews_corrcoef, balanced_accuracy_score)
import xgboost as xgb
import lightgbm as lgb

warnings.filterwarnings("ignore")
SEED = 42
np.random.seed(SEED)

PROC = Path("data/processed")
RES = Path("results/step6")
RES.mkdir(parents=True, exist_ok=True)

PREFILTER_K = 2000
MAX_PANEL = 200
N_SEARCH = 12
N_OUTER = 5
N_INNER = 3


# --------------------------------------------------------------------------
# DeLong's test for two correlated ROC-AUCs (Sun & Xu 2014 fast implementation)
# --------------------------------------------------------------------------
def _compute_midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T
    return T2


def _fast_delong(preds_sorted, m):
    n = preds_sorted.shape[1] - m
    pos = preds_sorted[:, :m]
    neg = preds_sorted[:, m:]
    k = preds_sorted.shape[0]
    tx = np.empty([k, m]); ty = np.empty([k, n]); tz = np.empty([k, m + n])
    for r in range(k):
        tx[r] = _compute_midrank(pos[r]); ty[r] = _compute_midrank(neg[r])
        tz[r] = _compute_midrank(preds_sorted[r])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    sx = np.cov(v01); sy = np.cov(v10)
    delongcov = sx / m + sy / n
    return aucs, delongcov


def delong_test(y_true, p1, p2):
    """Two-sided p-value for H0: AUC(p1) == AUC(p2). Returns (p_value, auc1, auc2)."""
    y_true = np.asarray(y_true)
    order = np.argsort(-y_true)            # positives first
    label_1 = y_true[order]
    m = int(label_1.sum())
    preds = np.vstack((p1, p2))[:, order]
    aucs, cov = _fast_delong(preds, m)
    l = np.array([[1, -1]])
    z_var = l.dot(cov).dot(l.T)
    if z_var <= 0:
        return (1.0, aucs[0], aucs[1])    # identical predictions
    z = (aucs[0] - aucs[1]) / np.sqrt(z_var)
    p = 2 * norm.sf(abs(z[0][0] if np.ndim(z) else z))
    return (float(p), float(aucs[0]), float(aucs[1]))


# --------------------------------------------------------------------------
def select_features(Xtr, ytr, genes):
    F, _ = f_classif(Xtr, ytr); F = np.nan_to_num(F, nan=0.0)
    top = np.argsort(F)[::-1][:PREFILTER_K]
    Xs = StandardScaler().fit_transform(Xtr[:, top])
    l1 = LogisticRegression(penalty="l1", solver="liblinear", C=0.1,
                            random_state=SEED, max_iter=1000).fit(Xs, ytr)
    nz = np.abs(l1.coef_[0]) > 1e-8
    sel = top[nz]
    if len(sel) > MAX_PANEL:
        sel = sel[np.argsort(np.abs(l1.coef_[0][nz]))[::-1][:MAX_PANEL]]
    if len(sel) == 0:
        sel = top[:20]
    return sel, [genes[i] for i in sel]


def grids():
    return {
        "ElasticNet": list(ParameterSampler(
            {"C": loguniform(1e-2, 10), "l1_ratio": [0.1, 0.5, 0.9]},
            n_iter=N_SEARCH, random_state=SEED)),
        "RandomForest": list(ParameterSampler(
            {"n_estimators": [100, 300, 500], "max_depth": [3, 5, 8, 12],
             "min_samples_leaf": [1, 2, 4]}, n_iter=N_SEARCH, random_state=SEED)),
        "SVM": list(ParameterSampler(
            {"C": loguniform(1e-1, 100), "gamma": loguniform(1e-4, 1)},
            n_iter=N_SEARCH, random_state=SEED)),
        "XGBoost": list(ParameterSampler(
            {"n_estimators": [100, 300, 500], "max_depth": [2, 4, 6],
             "learning_rate": loguniform(1e-2, 0.3), "subsample": [0.7, 0.9, 1.0]},
            n_iter=N_SEARCH, random_state=SEED)),
        "LightGBM": list(ParameterSampler(
            {"n_estimators": [100, 300, 500], "num_leaves": [15, 31, 63],
             "learning_rate": loguniform(1e-2, 0.3)},
            n_iter=N_SEARCH, random_state=SEED)),
        "MLP": list(ParameterSampler(
            {"hidden_layer_sizes": [(32,), (64,), (64, 32)],
             "alpha": loguniform(1e-5, 1e-1)}, n_iter=N_SEARCH, random_state=SEED)),
    }


def make(model, p, balanced=False):
    cw = "balanced" if balanced else None
    if model == "ElasticNet":
        return LogisticRegression(penalty="elasticnet", solver="saga",
                                  C=p["C"], l1_ratio=p["l1_ratio"], max_iter=1000,
                                  class_weight=cw, random_state=SEED)
    if model == "RandomForest":
        return RandomForestClassifier(random_state=SEED, n_jobs=-1,
                                      class_weight=cw, **p)
    if model == "SVM":
        return SVC(probability=True, random_state=SEED, class_weight=cw, **p)
    if model == "XGBoost":
        spw = (np.sum(_Y == 0) / np.sum(_Y == 1)) if balanced else 1
        return xgb.XGBClassifier(eval_metric="logloss", random_state=SEED,
                                 verbosity=0, scale_pos_weight=spw, **p)
    if model == "LightGBM":
        return lgb.LGBMClassifier(random_state=SEED, verbose=-1,
                                  class_weight=cw, **p)
    if model == "MLP":
        # MLP has no class_weight; imbalance handled by other models for comparison
        return MLPClassifier(max_iter=300, early_stopping=True,
                             random_state=SEED, **p)


def inner_best(model, grid, Xtr, ytr, balanced=False):
    inner = StratifiedKFold(N_INNER, shuffle=True, random_state=SEED)
    best, bp = -1, grid[0]
    for p in grid:
        sc = []
        for tri, vai in inner.split(Xtr, ytr):
            est = make(model, p, balanced).fit(Xtr[tri], ytr[tri])
            sc.append(roc_auc_score(ytr[vai], est.predict_proba(Xtr[vai])[:, 1]))
        if np.mean(sc) > best:
            best, bp = np.mean(sc), p
    return make(model, bp, balanced).fit(Xtr, ytr), bp


MODELS = ["ElasticNet", "RandomForest", "SVM", "XGBoost", "LightGBM", "MLP"]
_Y = None   # module-level for scale_pos_weight


def jac(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a | b) else 0.0


def evaluate_variant(X, y, genes, balanced):
    """Run the full nested-CV benchmark for one variant (balanced or not)."""
    tag = "weighted" if balanced else "unweighted"
    outer = StratifiedKFold(N_OUTER, shuffle=True, random_state=SEED)
    fold_metrics = {m: [] for m in MODELS}
    oof = {m: {"y": [], "p": []} for m in MODELS}   # pooled out-of-fold preds
    fold_sets, counter = [], {}

    for fi, (tri, tei) in enumerate(outer.split(X, y), 1):
        sel_idx, sel_genes = select_features(X[tri], y[tri], genes)
        if not balanced:                       # count features once (same selection)
            fold_sets.append(set(sel_genes))
            for g in sel_genes:
                counter[g] = counter.get(g, 0) + 1
        sc = StandardScaler().fit(X[tri][:, sel_idx])
        Xtr, Xte = sc.transform(X[tri][:, sel_idx]), sc.transform(X[tei][:, sel_idx])
        ytr, yte = y[tri], y[tei]
        for m in MODELS:
            global _Y; _Y = ytr
            est, _ = inner_best(m, GRIDS[m], Xtr, ytr, balanced)
            proba = est.predict_proba(Xte)[:, 1]
            pred = (proba >= 0.5).astype(int)
            tn, fp, fn, tp = confusion_matrix(yte, pred, labels=[0, 1]).ravel()
            fold_metrics[m].append({
                "fold": fi,
                "auc": roc_auc_score(yte, proba),
                "auprc": average_precision_score(yte, proba),
                "mcc": matthews_corrcoef(yte, pred),
                "balanced_acc": balanced_accuracy_score(yte, pred),
                "sensitivity": recall_score(yte, pred, pos_label=1, zero_division=0),
                "specificity": tn / (tn + fp) if (tn + fp) else np.nan,
                "f1": f1_score(yte, pred, zero_division=0),
            })
            oof[m]["y"] += yte.tolist(); oof[m]["p"] += proba.tolist()
            del est; gc.collect()
        print(f"  [{tag}] fold {fi} done", flush=True)

    rows = []
    for m in MODELS:
        d = pd.DataFrame(fold_metrics[m])
        rows.append({"model": m, "variant": tag,
                     "auc_mean": d.auc.mean(), "auc_sd": d.auc.std(),
                     "auprc_mean": d.auprc.mean(), "mcc_mean": d.mcc.mean(),
                     "balanced_acc_mean": d.balanced_acc.mean(),
                     "sens_mean": d.sensitivity.mean(),
                     "spec_mean": d.specificity.mean(), "f1_mean": d.f1.mean()})
    return pd.DataFrame(rows), oof, fold_metrics, fold_sets, counter


def main():
    expr = pd.read_csv(PROC / "TCGA_expression_common.csv", index_col=0)
    labels = pd.read_csv(PROC / "TCGA_labels.csv", index_col=0).loc[expr.columns]
    genes = expr.index.tolist()
    X = expr.T.values
    y = (labels["tumor_vs_normal"].values == "tumor").astype(int)
    print(f"Data: {X.shape[0]} samples x {X.shape[1]} genes | "
          f"tumour={y.sum()}, normal={(y==0).sum()} (ratio {y.sum()/(y==0).sum():.1f}:1)",
          flush=True)

    global GRIDS; GRIDS = grids()

    print("\n=== Unweighted variant ===", flush=True)
    sum_u, oof_u, fm_u, fold_sets, counter = evaluate_variant(X, y, genes, False)
    print("\n=== Class-weighted variant (imbalance correction) ===", flush=True)
    sum_w, oof_w, fm_w, _, _ = evaluate_variant(X, y, genes, True)

    summary = pd.concat([sum_u, sum_w], ignore_index=True)
    summary.to_csv(RES / "model_benchmark_summary.csv", index=False)
    print("\n=== Multi-metric benchmark ===", flush=True)
    print(summary.to_string(index=False), flush=True)

    # ---- principled model selection on the UNWEIGHTED variant ----
    # 1) rank by a composite of imbalance-robust metrics
    su = sum_u.copy()
    su["rank_score"] = su[["auc_mean", "auprc_mean", "mcc_mean",
                           "balanced_acc_mean", "f1_mean"]].mean(axis=1)
    su = su.sort_values("rank_score", ascending=False).reset_index(drop=True)
    top_model = su.iloc[0]["model"]

    # 2) DeLong test: which models are statistically indistinguishable from top?
    delong_rows = []
    yb = np.array(oof_u[top_model]["y"])
    pb = np.array(oof_u[top_model]["p"])
    for m in MODELS:
        if m == top_model:
            delong_rows.append({"model": m, "vs_top_auc_p": np.nan,
                                "statistically_tied_with_top": True})
            continue
        pm = np.array(oof_u[m]["p"])
        pval, a1, a2 = delong_test(yb, pb, pm)
        delong_rows.append({"model": m, "vs_top_auc_p": pval,
                            "statistically_tied_with_top": (pval > 0.05)})
    delong_df = pd.DataFrame(delong_rows)
    delong_df.to_csv(RES / "delong_vs_top.csv", index=False)

    # 3) among statistically-tied models, pick the most parsimonious
    PARSIMONY_RANK = {"ElasticNet": 1, "SVM": 2, "LightGBM": 3,
                      "RandomForest": 4, "XGBoost": 5, "MLP": 6}
    tied = delong_df.loc[delong_df["statistically_tied_with_top"], "model"].tolist()
    winner = min(tied, key=lambda m: PARSIMONY_RANK.get(m, 99))

    print("\n=== Model selection ===", flush=True)
    print(f"  Top by composite metric: {top_model}", flush=True)
    print(f"  Statistically tied with top (DeLong p>0.05): {tied}", flush=True)
    print(f"  Selected (most parsimonious among tied): {winner}", flush=True)

    # ---- feature stability (from unweighted selection) ----
    jpairs = [jac(a, b) for a, b in combinations(fold_sets, 2)]
    freq = (pd.Series(counter, name="folds_selected")
            .sort_values(ascending=False).to_frame())
    freq["selection_frequency"] = freq.folds_selected / N_OUTER
    freq.to_csv(RES / "feature_selection_frequency.csv")
    stable = freq.index[freq.folds_selected >= (N_OUTER // 2 + 1)].tolist()
    stab = {"mean_pairwise_jaccard": float(np.mean(jpairs)),
            "sd_pairwise_jaccard": float(np.std(jpairs)),
            "n_genes_selected_any_fold": len(counter),
            "n_genes_in_all_folds": int((freq.folds_selected == N_OUTER).sum()),
            "n_stable_genes_majority": len(stable)}
    json.dump(stab, open(RES / "feature_stability.json", "w"), indent=2)
    print("\n=== Feature stability ===", flush=True)
    for k, v in stab.items():
        print(f"  {k}: {v}", flush=True)

    # ---- selection decision record ----
    decision = {
        "winner": winner,
        "top_by_composite": top_model,
        "statistically_tied_models": tied,
        "selection_rule": ("Among models statistically indistinguishable from the "
                           "top performer by DeLong's test (p>0.05) across pooled "
                           "out-of-fold predictions, the most parsimonious model was "
                           "selected (Occam's razor), favouring interpretability for "
                           "the downstream SHAP analysis."),
        "imbalance_ratio": float(y.sum() / (y == 0).sum()),
        "imbalance_handling": ("Class-weighted variant benchmarked alongside "
                               "unweighted; SMOTE deliberately avoided "
                               "(high-dimensional, leakage-prone)."),
    }
    json.dump(decision, open(RES / "selection_decision.json", "w"), indent=2)

    # ---- refit final model on FULL cohort with stable panel ----
    idx = [genes.index(g) for g in stable]
    scf = StandardScaler().fit(X[:, idx])
    Xf = scf.transform(X[:, idx])
    global _Y; _Y = y
    final, best_p = inner_best(winner, GRIDS[winner], Xf, y, balanced=False)
    pickle.dump({"model": final, "winner": winner, "scaler": scf,
                 "genes": stable, "gene_indices": idx, "best_params": best_p},
                open(RES / "final_model.pkl", "wb"))
    pd.Series(stable, name="gene").to_csv(RES / "final_feature_panel.csv", index=False)
    pd.concat([pd.DataFrame(fm_u[m]).assign(model=m, variant="unweighted") for m in MODELS] +
              [pd.DataFrame(fm_w[m]).assign(model=m, variant="weighted") for m in MODELS],
              ignore_index=True).to_csv(RES / "per_fold_metrics.csv", index=False)
    json.dump({m: oof_u[m] for m in MODELS}, open(RES / "roc_data.json", "w"))
    print(f"\nWinner: {winner} | stable panel: {len(stable)} genes", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
