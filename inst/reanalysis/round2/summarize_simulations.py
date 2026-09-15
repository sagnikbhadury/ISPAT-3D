#!/usr/bin/env python
"""Refit MSFA on saved simulation GP covariance summaries and remake core figures."""
from pathlib import Path
import sys
import json
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fit_msfa_covariances import fit_one

ROOT=Path(sys.argv[1])
ROWS=[]
def rv(a,b):
 a=np.asarray(a);b=np.asarray(b)
 return float(np.sum(a*b)/np.sqrt(np.sum(a*a)*np.sum(b*b)))
def edges(a,b,t=0.1):
 x=np.abs(a[np.triu_indices_from(a,1)])>t
 y=np.abs(b[np.triu_indices_from(b,1)])>t
 tp=int(np.sum(x&y));fp=int(np.sum(~x&y));fn=int(np.sum(x&~y))
 return (tp/(tp+fn) if tp+fn else np.nan,fp/(fp+tp) if fp+tp else np.nan)
for G in (15,20,25):
 for folder in sorted((ROOT/f"G_{G}").glob("*")):
  manifest=folder/"manifest.json"
  if not manifest.exists():continue
  meta=json.loads((folder/"provenance.json").read_text())
  if len(meta["config"])!=3:raise ValueError(folder)
  zones=[]
  for q in (1,2,3):
   cov=pd.read_csv(folder/f"cov_Zone_{q}.csv",index_col=0).to_numpy()
   zones.append({"name":f"Zone_{q}","n":meta["config"][q-1],"cov":cov})
  fitted=fit_one(zones,min(math.ceil(2*math.log(G)),G-1),steps=350,seed=2026)
  allnets=[("Shared",fitted["shared"],pd.read_csv(folder/"truth_Shared.csv",index_col=0).to_numpy())]
  allnets += [(f"Zone_{q}",fitted["full"][q-1],
               pd.read_csv(folder/f"truth_Zone_{q}.csv",index_col=0).to_numpy()) for q in (1,2,3)]
  for name,estimated,truth in allnets:
   power,fdr=edges(truth,estimated)
   ROWS.append({"G":G,"Kernel":meta["kernel"],"N1":meta["config"][0],
    "N2":meta["config"][1],"N3":meta["config"][2],"N_total":sum(meta["config"]),
    "Replicate":meta["replicate"],"Network":name,"RV":rv(truth,estimated),
    "Power":power,"FDR":fdr,"factor_calls":fitted["calls"],
    "GP_fallbacks":sum(str(s["status"]).startswith("fallback") for s in meta["statuses"])})
  (folder/"fit_metrics.json").write_text(json.dumps({"loss":fitted["loss"],
     "calls":fitted["calls"]},indent=2))
  print("FIT",G,meta["kernel"],meta["config"],meta["replicate"],flush=True)
df=pd.DataFrame(ROWS)
if len(df)!=720:raise ValueError(f"Expected 720 network rows; got {len(df)}")
df.to_csv(ROOT/"simulation_replicate_metrics.csv",index=False)
agg=df.groupby(["G","Kernel","N1","N2","N3","Network"],as_index=False).agg(
 RV_mean=("RV","mean"),RV_sd=("RV","std"),Power_mean=("Power","mean"),
 FDR_mean=("FDR","mean"),GP_fallbacks=("GP_fallbacks","max"))
agg.to_csv(ROOT/"simulation_config_metrics.csv",index=False)
colors={"Matern":"#2873a6","RBF":"#d17436"}
nets=["Shared","Zone_1","Zone_2","Zone_3"]
fig,axes=plt.subplots(3,4,figsize=(13,9),sharey=True,constrained_layout=True)
for gi,G in enumerate((15,20,25)):
 for ni,network in enumerate(nets):
  ax=axes[gi,ni]
  for kernel in ("Matern","RBF"):
   d=agg[(agg.G==G)&(agg.Network==network)&(agg.Kernel==kernel)].sort_values(
    ["N1","N2","N3"])
   x=np.arange(1,len(d)+1)+(0.09 if kernel=="RBF" else -0.09)
   ax.errorbar(x,d.RV_mean,yerr=d.RV_sd,fmt="o-",ms=3,capsize=2,
    color=colors[kernel],label=kernel)
  ax.set_title(f"G={G} {network}")
  ax.set_ylim(0,1.05);ax.set_xticks(range(1,7),[str(i) for i in range(1,7)])
  ax.grid(alpha=.2)
  if ni==0:ax.set_ylabel("Covariance RV")
  if gi==2:ax.set_xlabel("Cell-count configuration")
axes[0,0].legend()
fig.savefig(ROOT/"figure03_rv_simulation.pdf");plt.close(fig)
fig,axes=plt.subplots(3,4,figsize=(13,9),constrained_layout=True)
for gi,G in enumerate((15,20,25)):
 for ni,network in enumerate(nets):
  ax=axes[gi,ni]
  for kernel in ("Matern","RBF"):
   d=agg[(agg.G==G)&(agg.Network==network)&(agg.Kernel==kernel)].sort_values(
    ["N1","N2","N3"])
   x=np.arange(1,len(d)+1)+(0.09 if kernel=="RBF" else -0.09)
   ax.plot(x,d.Power_mean,"o-",ms=3,color=colors[kernel],label=f"{kernel} power")
   ax.plot(x,d.FDR_mean,"s--",ms=3,color=colors[kernel],label=f"{kernel} FDP")
  ax.set_title(f"G={G} {network}")
  ax.set_ylim(0,1.05);ax.set_xticks(range(1,7),[str(i) for i in range(1,7)])
  ax.grid(alpha=.2)
  if ni==0:ax.set_ylabel("Rate at |covariance| > 0.1")
  if gi==2:ax.set_xlabel("Cell-count configuration")
axes[0,0].legend(fontsize=7)
fig.savefig(ROOT/"figure04_metrics_simulation.pdf");plt.close(fig)
print("ALL SIMULATION FIGURES COMPLETE")
