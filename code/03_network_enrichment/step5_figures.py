#!/usr/bin/env python3
"""Step 5 figures - PPI network + triangulated enrichment + two-level. 600 DPI, caption-free."""
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

RES = Path("results/step5"); FIG = Path("figures"); FIG.mkdir(exist_ok=True)
plt.rcParams.update({"savefig.dpi":600,"font.family":"DejaVu Sans","font.size":12,
    "axes.titlesize":14,"axes.titleweight":"bold","axes.labelsize":12.5,
    "axes.labelweight":"bold","xtick.labelsize":10.5,"ytick.labelsize":10,
    "legend.fontsize":10,"axes.linewidth":1.1,"pdf.fonttype":42,"ps.fonttype":42})

def main():
    robust = pd.read_csv(RES/"triangulated_enrichment_robust.csv")
    panel_pw = pd.read_csv(RES/"panel_genes_in_robust_pathways.csv")
    twolevel = pd.read_csv(RES/"two_level_panel_in_network.csv")
    saudi = pd.read_csv(RES/"step8c_saudi_pathway_crossref.csv")

    fig = plt.figure(figsize=(17,13))
    gs = gridspec.GridSpec(2,2,figure=fig,hspace=0.42,wspace=0.32)

    axA = fig.add_subplot(gs[0,0])
    top3 = robust[robust['n_platforms']>=3].copy()
    keep = ['Wnt signaling pathway','Pathways in cancer','Cytokine-cytokine receptor interaction',
            'Calcium signaling pathway','Cell adhesion molecule (CAM) interaction','cAMP signaling pathway',
            'Chemokine signaling pathway','Complement and coagulation cascades','Protein digestion and absorption',
            'Neuroactive ligand-receptor interaction','extracellular matrix organization','Hemostasis']
    sub = top3[top3['pathway'].isin(keep)].drop_duplicates('pathway')
    if len(sub)<8: sub = top3.head(12)
    sub = sub.sort_values('n_platforms')
    axA.barh(range(len(sub)), sub['n_platforms'], color="#0072B2", edgecolor="black", linewidth=0.8)
    axA.set_yticks(range(len(sub))); axA.set_yticklabels(sub['pathway'], fontsize=9)
    axA.set_xlabel("# platforms significant (of 4)")
    axA.set_xlim(0,4.3); axA.axvline(2,ls='--',color='red',lw=1.3,label='Robust threshold (\u22652)')
    axA.set_title("(A) Triangulated robust pathways", loc="left")
    axA.legend(loc='lower right'); axA.grid(axis='x',alpha=0.3)

    axB = fig.add_subplot(gs[0,1])
    pp = panel_pw.sort_values('n_robust_pathways')
    cols = ['#009E73' if n>0 else '#BBBBBB' for n in pp['n_robust_pathways']]
    axB.barh(range(len(pp)), pp['n_robust_pathways'], color=cols, edgecolor="black", linewidth=0.8)
    axB.set_yticks(range(len(pp))); axB.set_yticklabels(pp['gene'])
    axB.set_xlabel("# robust enriched pathways containing gene")
    axB.set_title("(B) Two-level: panel genes in enriched pathways", loc="left")
    axB.grid(axis='x',alpha=0.3)
    for i,(_,r) in enumerate(pp.iterrows()):
        if r['gene']=='ESM1': axB.text(r['n_robust_pathways']+0.1,i,"angiogenesis",fontsize=8,va='center',style='italic')
        if r['gene']=='CDH3': axB.text(r['n_robust_pathways']+0.1,i,"Wnt/adhesion",fontsize=8,va='center',style='italic')

    axC = fig.add_subplot(gs[1,0])
    piv = twolevel.pivot_table(index='gene',columns='confidence',values='in_network',aggfunc='first')
    genes = piv.index.tolist()
    present_09 = [(1 if piv.loc[g,0.9]=='yes' else 0) for g in genes]
    present_07 = [(1 if piv.loc[g,0.7]=='yes' else 0) for g in genes]
    x = np.arange(len(genes)); w=0.38
    axC.bar(x-w/2, present_09, w, label='conf>0.9', color="#D55E00", edgecolor='black')
    axC.bar(x+w/2, present_07, w, label='conf>0.7', color="#E69F00", edgecolor='black')
    axC.set_xticks(x); axC.set_xticklabels(genes, rotation=45, ha='right', fontsize=9)
    axC.set_ylabel("Present in PPI network (1=yes)")
    axC.set_yticks([0,1]); axC.set_ylim(0,1.3)
    axC.set_title("(C) Panel genes peripheral to PPI core (sensitivity)", loc="left")
    axC.legend(); axC.grid(axis='y',alpha=0.3)

    axD = fig.add_subplot(gs[1,1])
    sd = saudi.copy()
    conv = (sd['convergent_with_our_enrichment']=='yes').astype(int)
    cols = ['#009E73' if c else '#BBBBBB' for c in conv]
    short = [s.split('/')[0].split(' signaling')[0] for s in sd['saudi_pathway']]
    axD.barh(range(len(sd)), conv, color=cols, edgecolor='black', linewidth=0.8)
    axD.set_yticks(range(len(sd))); axD.set_yticklabels(short, fontsize=10)
    axD.set_xlim(0,1.3); axD.set_xticks([0,1]); axD.set_xticklabels(['no','yes'])
    axD.set_xlabel("Convergent with our robust enrichment")
    axD.set_title("(D) Saudi CRC pathway cross-reference", loc="left")
    axD.grid(axis='x',alpha=0.3)

    fig.savefig(FIG/"Figure9_enrichment_network.png",dpi=600,bbox_inches="tight",facecolor="white")
    fig.savefig(FIG/"Figure9_enrichment_network.jpg",dpi=600,bbox_inches="tight",
                facecolor="white",pil_kwargs={"quality":95})
    print("Saved Figure9_enrichment_network.png/.jpg at 600 DPI")

if __name__=="__main__":
    main()
