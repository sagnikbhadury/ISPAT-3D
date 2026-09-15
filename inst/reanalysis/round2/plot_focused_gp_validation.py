#!/usr/bin/env python
from pathlib import Path
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=Path(sys.argv[1])
df=pd.read_csv(ROOT/"focused_gp_vecchia_validation.csv")
if len(df)!=640:raise ValueError(len(df))
metrics=[("Mean_abs_Moran","Residual 3D spatial autocorrelation"),
("Partial_correlation_RMSE","Partial-correlation RMSE"),
("Null_edge_FPR","False positive rate on true null pairs"),
("Signal_residual_RMSE","Adjusted residual RMSE to known signal")]
methods=["Matched 3D Vecchia","Short scale 3D Vecchia","Planar 2D Vecchia","No GP"]
colors=["#17689b","#dc7d30","#4b9a65","#a12c2d"]
styles=["-","--","-.",":"]
summary=df.groupby(["Kernel","Spatial_SD","Method"],as_index=False).agg(
 **{f"{metric}_mean":(metric,"mean") for metric,_ in metrics},
 **{f"{metric}_se":(metric,lambda x:x.std(ddof=1)/(len(x)**.5)) for metric,_ in metrics})
summary.to_csv(ROOT/"focused_gp_vecchia_summary.csv",index=False)
fig,axes=plt.subplots(4,2,figsize=(11.5,10.5),constrained_layout=True)
for ri,(metric,title) in enumerate(metrics):
 for ci,kernel in enumerate(("Matern","RBF")):
  ax=axes[ri,ci]
  for method,color,style in zip(methods,colors,styles):
   d=summary[(summary.Kernel==kernel)&(summary.Method==method)].sort_values("Spatial_SD")
   x=d.Spatial_SD.to_numpy();y=d[f"{metric}_mean"].to_numpy()
   se=d[f"{metric}_se"].to_numpy()
   ax.errorbar(x,y,yerr=se,marker="o",linestyle=style,color=color,
    linewidth=1.5,markersize=4,capsize=2,label=method)
  ax.set_title(f"{title} | {kernel}")
  ax.set_xticks([0,.5,1,1.5]);ax.set_xlabel("Shared spatial nuisance SD")
  ax.grid(alpha=.2)
  if ci==0:ax.set_ylabel("Lower is better")
axes[0,0].legend(fontsize=7)
fig.savefig(ROOT/"figure06_gp_vecchia_validation.pdf");plt.close(fig)
print(summary[summary.Spatial_SD==1.5][["Kernel","Method",
 "Mean_abs_Moran_mean","Partial_correlation_RMSE_mean","Null_edge_FPR_mean"]].to_string(index=False))
