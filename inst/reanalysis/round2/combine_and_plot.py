#!/usr/bin/env python
"""Combine three analysis tiers and plot every cell pair across relative tumor zones."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ZONES = ["Very Low", "Low", "Intermediate", "High", "Very High"]
NETWORKS = ["Shared"] + ZONES
TAGS = {q: q.replace(" ", "_") for q in NETWORKS}
LABELS = {
    "CD8_T":"CD8 T", "CD4_T":"CD4 T", "T_cell":"other T",
    "B_cell":"B", "Other_immune":"other immune",
    "Tumor_HER2pos":"HER2 tumor", "Tumor_basal":"basal tumor",
    "Tumor_luminal":"luminal tumor", "Tumor_other":"other tumor",
    "CD8_T_cell":"CD8 T", "CD4_T_cell":"CD4 T",
    "Plasma_cell":"plasma", "Myoepithelial":"myoepithelial",
}


def read_pcor(folder, tag):
    return pd.read_csv(folder / "factor_fit" / f"pcor_{tag}.csv", index_col=0)


def run(root, budgets, dataset):
    tiers = [root / f"tier_{n}" for n in budgets]
    out = root / "combined"
    out.mkdir(parents=True, exist_ok=True)
    matrices = {}
    sd_matrices = {}
    for network in NETWORKS:
        mats = [read_pcor(t, TAGS[network]) for t in tiers]
        labels = mats[0].index.tolist()
        if any(m.index.tolist() != labels or m.columns.tolist() != labels for m in mats):
            raise ValueError("Cell-type order differs across tiers")
        stack = np.stack([m.to_numpy() for m in mats])
        z = np.arctanh(np.clip(stack, -0.999999, 0.999999))
        pcor = np.tanh(z.mean(0))
        np.fill_diagonal(pcor, 1)
        sd = stack.std(0, ddof=1)
        np.fill_diagonal(sd, 0)
        matrices[network] = pcor
        sd_matrices[network] = sd
        pd.DataFrame(pcor,index=labels,columns=labels).to_csv(out / f"pcor_combined_{TAGS[network]}.csv")
        pd.DataFrame(sd,index=labels,columns=labels).to_csv(out / f"pcor_sd_{TAGS[network]}.csv")
    g = len(labels)
    pairs = [(i,j) for i in range(g) for j in range(i+1,g)]
    effects = np.array([[matrices[q][i,j] for q in ZONES] for i,j in pairs])
    signs = np.zeros_like(effects)
    ranges = np.zeros_like(effects)
    for zidx, network in enumerate(ZONES):
        stack = np.stack([read_pcor(t,TAGS[network]).to_numpy() for t in tiers])
        for pidx,(i,j) in enumerate(pairs):
            vals = stack[:,i,j]
            signs[pidx,zidx] = max(np.sum(vals>0),np.sum(vals<0))
            ranges[pidx,zidx] = vals.max()-vals.min()
    pair_names = [f"{LABELS.get(labels[i],labels[i])} : {LABELS.get(labels[j],labels[j])}" for i,j in pairs]
    order = np.argsort(-np.max(np.abs(effects),axis=1))
    effects = effects[order]; signs = signs[order]; ranges = ranges[order]
    pair_names = [pair_names[i] for i in order]
    pd.DataFrame({"pair":pair_names,**{f"rho_{q}":effects[:,k] for k,q in enumerate(ZONES)},
                  **{f"sign_tiers_{q}":signs[:,k].astype(int) for k,q in enumerate(ZONES)},
                  **{f"range_{q}":ranges[:,k] for k,q in enumerate(ZONES)}}).to_csv(
        out / "edge_effects_and_tier_stability.csv",index=False)
    height = max(11,0.35*len(pairs)+2.2)
    fig,(ax,ax2) = plt.subplots(1,2,figsize=(14,height),
                               gridspec_kw={"width_ratios":[5,2.2]},
                               constrained_layout=True)
    im = ax.imshow(effects,aspect="auto",cmap="RdBu",norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1))
    ax.set_xticks(range(5),ZONES,rotation=35,ha="right")
    ax.set_yticks(range(len(pairs)),pair_names,fontsize=7.3 if dataset=="crc" else 6.1)
    ax.tick_params(length=0)
    ax.set_title(f"{dataset.upper()} conditional partial correlations across relative tumor-density zones")
    for i in range(len(pairs)):
        for j in range(5):
            val = effects[i,j]
            ax.text(j,i,f"{val:+.2f}",ha="center",va="center",
                    fontsize=7 if dataset=="crc" else 5.5,
                    color="white" if abs(val)>0.55 else "black")
    ax2.imshow(signs,aspect="auto",cmap="Greys",vmin=1,vmax=3)
    ax2.set_xticks(range(5),ZONES,rotation=35,ha="right")
    ax2.set_yticks([])
    ax2.tick_params(length=0)
    ax2.set_title("Same sign in tiers")
    for i in range(len(pairs)):
        for j in range(5):
            ax2.text(j,i,f"{int(signs[i,j])}/3",ha="center",va="center",
                     fontsize=7 if dataset=="crc" else 5.5,
                     color="white" if signs[i,j]==3 else "black")
    fig.colorbar(im,ax=ax,label="Partial correlation",fraction=0.015,pad=0.01)
    fig.savefig(out / f"figure_{dataset}_pair_heatmap.pdf")
    plt.close(fig)
    print(dataset,"combined",len(pairs),"pairs; median sign agreement",
          np.median(signs),"3/3 fraction",np.mean(signs==3))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,required=True)
    p.add_argument("--budgets",type=int,nargs=3,required=True)
    p.add_argument("--dataset",choices=["crc","bc"],required=True)
    a=p.parse_args()
    run(a.root,a.budgets,a.dataset)
