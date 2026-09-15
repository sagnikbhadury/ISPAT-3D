#!/usr/bin/env python
"""Plot matched 3D minus section-wise planar partial correlations."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from combine_and_plot import ZONES,LABELS
def run(dataset,budget,root_arg):
 root=Path(root_arg)/dataset
 base=root/f"sectionwise_2d_min10_tier_{budget}"
 out=root/"combined";out.mkdir(exist_ok=True)
 mats=[]
 records=[]
 for zone in ZONES:
  tag=zone.replace(" ","_")
  a=pd.read_csv(out/f"pcor_combined_{tag}.csv",index_col=0)
  b=pd.read_csv(base/"factor_fit"/f"pcor_{tag}.csv",index_col=0).loc[a.index,a.columns]
  d=a-b
  d.to_csv(out/f"pcor_3d_minus_sectionwise_2d_{tag}.csv")
  labels=[x for x in a.index if not (dataset=="bc" and x=="Tumor_luminal")]
  mats.append(d.loc[labels,labels].to_numpy())
 g=len(labels)
 pairs=[(i,j) for i in range(g) for j in range(i+1,g)]
 effect=np.array([[m[i,j] for m in mats] for i,j in pairs])
 names=[f"{LABELS.get(labels[i],labels[i])} : {LABELS.get(labels[j],labels[j])}" for i,j in pairs]
 order=np.argsort(-np.max(np.abs(effect),axis=1))
 effect=effect[order];names=[names[i] for i in order]
 pd.DataFrame({"pair":names,**{f"delta_{z}":effect[:,k] for k,z in enumerate(ZONES)}}).to_csv(
  out/"three_d_minus_matched_sectionwise_2d_edges.csv",index=False)
 vmax=max(0.35,float(np.ceil(np.max(np.abs(effect))*10)/10))
 fig,ax=plt.subplots(figsize=(10,max(10,len(pairs)*0.35+2)),constrained_layout=True)
 im=ax.imshow(effect,aspect="auto",cmap="RdBu",norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax))
 ax.set_yticks(range(len(pairs)),names,fontsize=7 if dataset=="crc" else 6)
 ax.set_xticks(range(5),ZONES,rotation=30,ha="right")
 ax.set_title(f"{dataset.upper()} volumetric minus matched section-wise planar estimate")
 ax.tick_params(length=0)
 for i in range(len(pairs)):
  for j in range(5):
   val=effect[i,j]
   ax.text(j,i,f"{val:+.2f}",ha="center",va="center",
    fontsize=7 if dataset=="crc" else 5.5,
    color="white" if abs(val)>vmax*.55 else "black")
 fig.colorbar(im,ax=ax,label="Difference in partial correlation")
 fig.savefig(out/f"figure_{dataset}_3d_vs_2d_matched.pdf");plt.close(fig)
 report={"dataset":dataset,"baseline_budget":budget,"kind":"matched selected cells, section-wise XY GP vs pooled XYZ GP",
         "n_pairs":len(pairs),"max_abs_delta":float(np.max(np.abs(effect))),
         "median_abs_delta":float(np.median(np.abs(effect))),
         "top_pairs":[{"pair":names[i],"deltas":dict(zip(ZONES,effect[i].tolist()))} for i in range(min(10,len(pairs)))]}
 (out/"three_d_vs_2d_summary.json").write_text(json.dumps(report,indent=2))
 print(dataset,"max delta",report["max_abs_delta"],"median",report["median_abs_delta"],
       "top",names[0],flush=True)
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--dataset",choices=["crc","bc"],required=True)
 p.add_argument("--budget",type=int,required=True);p.add_argument("--root",type=Path,required=True);a=p.parse_args();run(a.dataset,a.budget,a.root)
