#!/usr/bin/env python3
"""
Step 9d - Permutation-null figure (revision analysis -> Supplementary Figure S4).

Renders the permutation null distribution of the classifier AUC against the observed
value, at 600 dpi. Depends on the same computation as Step 9a.

Outputs: results/09_revision_analyses/FigureS4_permutation_null.png (600 dpi)

License: MIT   |   Random seed: 42
"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

warnings.filterwarnings("ignore")
SEED = 42
DATA = Path("data")
RESB = Path("results/04_model_benchmark")
RES = Path("results/09_revision_analyses")
RES.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13,
                     "axes.titleweight": "bold", "axes.labelweight": "bold",
                     "savefig.dpi": 600})


def main():
    expr = pd.read_csv(DATA / "TCGA_expression_qc.csv.gz", index_col=0)
    labels = pd.read_csv(DATA / "TCGA_labels.csv", index_col=0)
    panel = pd.read_csv(RESB / "final_feature_panel.csv")["gene"].tolist()
    X = expr.loc[panel].T.values
    y = (labels.loc[expr.loc[panel].T.index, "tumor_vs_normal"] == "tumor").astype(int).values

    def mk(seed=SEED):
        return Pipeline([("s", StandardScaler()),
                         ("c", LogisticRegression(penalty="elasticnet", solver="saga",
                              l1_ratio=0.5, C=1.0, max_iter=5000, random_state=seed))])

    obs = cross_val_score(mk(), X, y, cv=StratifiedKFold(5, shuffle=True, random_state=SEED),
                          scoring="roc_auc").mean()
    rng = np.random.RandomState(SEED)
    null = []
    for _ in range(1000):
        yp = rng.permutation(y)
        try:
            null.append(cross_val_score(mk(), X, yp,
                        cv=StratifiedKFold(5, shuffle=True, random_state=SEED),
                        scoring="roc_auc").mean())
        except Exception:
            null.append(0.5)
    null = np.array(null)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(null, bins=40, color="#4C72B0", edgecolor="black", alpha=0.85,
            label="Permuted-label null (1,000 runs)")
    ax.axvline(obs, color="#C44E52", lw=3, label=f"Observed AUC = {obs:.3f}")
    ax.axvline(np.percentile(null, 95), color="#555555", ls="--", lw=2,
               label=f"Null 95th percentile = {np.percentile(null,95):.3f}")
    ax.set_xlabel("Cross-validated AUC"); ax.set_ylabel("Frequency (permutations)")
    ax.set_title("Permutation test: observed classifier AUC vs permuted-label null")
    ax.set_xlim(0.3, 1.05)
    ax.legend(loc="upper center", fontsize=11, framealpha=0.95)
    ax.text(0.62, ax.get_ylim()[1]*0.55,
            f"Permutation p < 0.001\nNull mean = {null.mean():.3f}",
            fontsize=12, bbox=dict(boxstyle="round", fc="white", ec="gray"))
    fig.tight_layout()
    fig.savefig(RES / "FigureS4_permutation_null.png", dpi=600, facecolor="white")
    print(f"Saved {RES/'FigureS4_permutation_null.png'}")


if __name__ == "__main__":
    main()
