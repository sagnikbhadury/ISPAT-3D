#!/usr/bin/env python
"""Measure rerun GP stage wall time and factor fit time; remake runtime figure."""
from pathlib import Path
import json,math,time,sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fit_msfa_covariances import fit_one
root=Path(sys.argv[1])
rows=[]
for G in (15,20,25):
 for folder in sorted((root/f"G_{G}").glob("*")):
  manifest=folder/"manifest.json"
  if not manifest.exists():continue
  meta=json.loads((folder/"provenance.json").read_text())
  zones=[{"name":f"Zone_{q}","n":meta["config"][q-1],
   "cov":pd.read_csv(folder/f"cov_Zone_{q}.csv",index_col=0).to_numpy()} for q in (1,2,3)]
  start=time.perf_counter()
  fit_one(zones,min(math.ceil(2*math.log(G)),G-1),steps=350,seed=2026)
  factor_s=time.perf_counter()-start
  gp_s=manifest.stat().st_mtime-folder.stat().st_ctime
  if not (0<gp_s<600):raise ValueError((folder,gp_s))
  rows.append({"G":G,"Kernel":meta["kernel"],"N1":meta["config"][0],
   "N2":meta["config"][1],"N3":meta["config"][2],"Replicate":meta["replicate"],
   "GP_stage_wall_seconds":gp_s,"Factor_fit_seconds":factor_s,
   "Combined_seconds":gp_s+factor_s})
if len(rows)!=180:raise ValueError(len(rows))
df=pd.DataFrame(rows);df.to_csv(root/"simulation_runtime_observed.csv",index=False)
agg=df.groupby(["G","Kernel","N1","N2","N3"],as_index=False).agg(
 GP_mean=("GP_stage_wall_seconds","mean"),GP_sd=("GP_stage_wall_seconds","std"),
 Factor_mean=("Factor_fit_seconds","mean"),Factor_sd=("Factor_fit_seconds","std"),
 Combined_mean=("Combined_seconds","mean"),Combined_sd=("Combined_seconds","std"))
agg.to_csv(root/"simulation_runtime_config.csv",index=False)
fig,axes=plt.subplots(3,2,figsize=(11.5,8.5),sharey=True,constrained_layout=True)
for ri,G in enumerate((15,20,25)):
 for ci,kernel in enumerate(("Matern","RBF")):
  ax=axes[ri,ci];d=agg[(agg.G==G)&(agg.Kernel==kernel)].sort_values(["N1","N2","N3"])
  x=np.arange(1,len(d)+1)
  ax.errorbar(x,d.Combined_mean,yerr=d.Combined_sd,fmt="o-",color="#2873a6",capsize=3)
  ax.plot(x,d.Factor_mean,"s--",color="#d17436",label="Factor fit only")
  ax.set_title(f"G={G} | generating kernel {kernel}")
  ax.set_xticks(x,[str(i) for i in x]);ax.set_xlabel("Cell-count configuration")
  ax.grid(alpha=.2)
  if ci==0:ax.set_ylabel("Observed seconds")
axes[0,0].legend(fontsize=8)
fig.savefig(root/"figure05_runtime_simulation.pdf");plt.close(fig)
print("runtime GP",df.GP_stage_wall_seconds.min(),df.GP_stage_wall_seconds.mean(),
      df.GP_stage_wall_seconds.max(),"factor",df.Factor_fit_seconds.mean(),flush=True)
