#!/usr/bin/env python
"""Full-pair zone heatmaps with tier and section-block stability."""
import argparse
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from combine_and_plot import ZONES,LABELS
def run(dataset,root_arg):
 root=Path(root_arg)/dataset/"combined"
 boot=pd.read_csv(root/"section_block_resampling.csv")
 edge=pd.read_csv(root/"edge_effects_and_tier_stability.csv")
 mats=[pd.read_csv(root/f"pcor_combined_{z.replace(' ','_')}.csv",index_col=0) for z in ZONES]
 labels=mats[0].index.tolist()
 pairs=[(i,j) for i in range(len(labels)) for j in range(i+1,len(labels))]
 effects=np.array([[m.iloc[i,j] for m in mats] for i,j in pairs])
 names=[f"{LABELS.get(labels[i],labels[i])} : {LABELS.get(labels[j],labels[j])}" for i,j in pairs]
 order=np.argsort(-np.max(np.abs(effects),axis=1))
 effects=effects[order];pairs=[pairs[k] for k in order];names=[names[k] for k in order]
 tier=edge.set_index("pair").loc[names,[f"sign_tiers_{z}" for z in ZONES]].to_numpy()
 block=np.zeros_like(effects)
 for k,(i,j) in enumerate(pairs):
  for l,z in enumerate(ZONES):
   row=boot[(boot.zone==z)&(boot.cell_1==labels[i])&(boot.cell_2==labels[j])]
   if len(row)!=1:raise ValueError((labels[i],labels[j],z,len(row)))
   block[k,l]=row.iloc[0].same_sign_fraction
 height=max(11,.35*len(pairs)+2.2)
 fig,(ax,ax2,ax3)=plt.subplots(1,3,figsize=(18,height),
  gridspec_kw={"width_ratios":[5,2.0,2.4]},constrained_layout=True)
 im=ax.imshow(effects,aspect="auto",cmap="RdBu",norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1))
 ax.set_xticks(range(5),ZONES,rotation=32,ha="right")
 ax.set_yticks(range(len(pairs)),names,fontsize=7 if dataset=="crc" else 5.8)
 ax.tick_params(length=0);ax.set_title("Partial correlation by relative density zone")
 for i in range(len(pairs)):
  for j in range(5):
   val=effects[i,j]
   ax.text(j,i,f"{val:+.2f}",ha="center",va="center",
    fontsize=7 if dataset=="crc" else 5.5,
    color="white" if abs(val)>.55 else "black")
 ax2.imshow(tier,aspect="auto",cmap="Greys",vmin=1,vmax=3)
 ax2.set_xticks(range(5),ZONES,rotation=32,ha="right")
 ax2.set_yticks([]);ax2.tick_params(length=0);ax2.set_title("Matching sign\nacross three budgets")
 for i in range(len(pairs)):
  for j in range(5):
   ax2.text(j,i,f"{int(tier[i,j])}/3",ha="center",va="center",
    fontsize=7 if dataset=="crc" else 5.5,
    color="white" if tier[i,j]==3 else "black")
 ax3.imshow(block,aspect="auto",cmap="YlGnBu",vmin=0,vmax=1)
 ax3.set_xticks(range(5),ZONES,rotation=32,ha="right")
 ax3.set_yticks([]);ax3.tick_params(length=0)
 ax3.set_title("Same sign across 100\nsection-block resamples")
 for i in range(len(pairs)):
  for j in range(5):
   ax3.text(j,i,f"{100*block[i,j]:.0f}%",ha="center",va="center",
    fontsize=7 if dataset=="crc" else 5.5,
    color="white" if block[i,j]>.68 else "black")
 fig.colorbar(im,ax=ax,label="Partial correlation",fraction=.015,pad=.01)
 fig.savefig(root/f"figure_{dataset}_pair_heatmap.pdf");plt.close(fig)
 print(dataset,"all",len(pairs),"pairs with bootstrap stability",flush=True)
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--dataset",choices=["crc","bc"],required=True)
 p.add_argument('--root',type=Path,required=True)
 a=p.parse_args();run(a.dataset,a.root)
