# How to Host This Repository on Your GitHub Account

Two methods. **Method A (website, no commands)** is easiest if you are not familiar with
git. **Method B (command line)** is faster if you have git installed.

Your GitHub: **https://github.com/SAH-ML**  (suggested repo name: `crc-stage-stratified-xai`)

---

## Method A — Upload via the GitHub website (no git needed)

1. Go to https://github.com/new (sign in as SAH-ML).
2. **Repository name:** `crc-stage-stratified-xai`
3. **Description:** "Stage-stratified explainable AI for colorectal cancer transcriptomics — code, data, and reproducibility for the Diagnostics manuscript."
4. Choose **Public**. Do NOT tick "Add a README" (we already have one). Click **Create repository**.
5. On the new empty repo page, click **"uploading an existing file"** (the link in the
   "Quick setup" box).
6. Drag the **entire contents** of this folder (the `code/`, `data/`, `results/`,
   `figures/`, `docs/` folders plus `README.md`, `requirements.txt`, `LICENSE`,
   `.gitignore`) into the upload area.
   - If the browser struggles with many files at once, upload folder-by-folder:
     first `code/`, commit; then `results/`, commit; then `figures/`, etc.
7. At the bottom, write a commit message ("Initial release") and click **Commit changes**.
8. Done — your repo is live at `https://github.com/SAH-ML/crc-stage-stratified-xai`.

> **Large files:** the 600-dpi figures total ~70 MB. GitHub's per-file limit is 100 MB,
> so all individual figures are fine. If the total repo feels large, you can instead
> attach the `figures/` folder as a **Release asset** (see below) and keep the repo lean.

## Method B — Command line (git)

```bash
# from inside this folder:
git init
git add .
git commit -m "Initial release: code, data pointers, results, figures"
git branch -M main
git remote add origin https://github.com/SAH-ML/crc-stage-stratified-xai.git
git push -u origin main
```
If prompted, authenticate with a GitHub Personal Access Token (Settings → Developer
settings → Personal access tokens).

---

## Recommended: attach large data as a Release (keeps the repo small)

The QC'd expression matrix (`TCGA_expression_qc.csv.gz`, ~17 MB) and optionally the
figure bundle are best attached to a **GitHub Release** rather than committed:

1. On your repo page → **Releases** → **Create a new release**.
2. Tag: `v1.0`. Title: "Manuscript release v1.0".
3. Drag `TCGA_expression_qc.csv.gz` (and optionally a zip of `figures/`) into the
   **"Attach binaries"** area.
4. **Publish release.** Users then download these from the Releases tab.
5. (Optional) update `data/README.md` to point to the Release asset URL.

---

## After uploading — finishing touches

- Add **topics** (repo home → ⚙ next to About): `colorectal-cancer`, `explainable-ai`,
  `shap`, `transcriptomics`, `machine-learning`, `bioinformatics`.
- Confirm the **LICENSE** shows as "MIT" in the repo sidebar.
- Put the final repository URL into the manuscript's **Data Availability Statement**
  before submission.
- (Optional) enable **Zenodo** integration to mint a DOI for the repo, which journals
  like to see in the Data Availability statement.
