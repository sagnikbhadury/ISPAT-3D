#!/usr/bin/env python
"""Section-block resampling of fixed GP residual covariance summaries."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from fit_msfa_covariances import fit_one,partial_corr
def run(dataset,budget,reps,root_arg):
 root=Path(root_arg)/dataset
 data=json.loads((root/f"section_block_stats_{budget}.json").read_text())
 sections=np.array(data["sections"],dtype=int);labels=data["labels"]
 zones=[]
 for zone in data["zones"]:
  records={int(x["section"]):{"n":int(x["n"]),"sum":np.array(x["sum"],float),
   "cross":np.array(x["cross"],float)} for x in zone["records"]}
  if set(records)!=set(sections):raise ValueError("Section absent from a zone")
  zones.append((zone["name"],records))
 rng=np.random.default_rng(20260915)
 active=np.array([i for i,label in enumerate(labels) if label!="Tumor_luminal"]) if dataset=="bc" else np.arange(len(labels))
 baseline={}
 for name,_ in zones:
  tag=name.replace(" ","_")
  baseline[name]=pd.read_csv(root/"combined"/f"pcor_combined_{tag}.csv",index_col=0).loc[labels,labels].to_numpy()
 boot=np.zeros((reps,5,len(labels),len(labels)))
 for r in range(reps):
  drawn=rng.choice(sections,size=len(sections),replace=True)
  counts=np.bincount(drawn,minlength=int(sections.max())+1)
  inputs=[]
  for name,records in zones:
   n=0;sy=np.zeros(len(labels));sx=np.zeros((len(labels),len(labels)))
   for sec in sections:
    weight=counts[sec]
    if weight:
     item=records[int(sec)]
     n+=weight*item["n"];sy+=weight*item["sum"];sx+=weight*item["cross"]
   cov=(sx-np.outer(sy,sy)/n)/(n-1)
   cov=(cov+cov.T)/2
   inputs.append({"name":name,"n":n,"cov":cov[np.ix_(active,active)]})
  fit=fit_one(inputs,min(5,len(active)-1),steps=350,seed=2026+r)
  for q in range(5):
   small=partial_corr(fit["full"][q])
   full=np.eye(len(labels));full[np.ix_(active,active)]=small
   boot[r,q]=full
  if (r+1)%25==0:print("BOOT",dataset,r+1,"/",reps,flush=True)
 rows=[]
 for q,(name,_) in enumerate(zones):
  for i in range(len(labels)):
   for j in range(i+1,len(labels)):
    vals=boot[:,q,i,j];value=baseline[name][i,j]
    rows.append({"dataset":dataset,"zone":name,"cell_1":labels[i],"cell_2":labels[j],
     "rho_combined":value,"rho_boot_median":float(np.median(vals)),
     "rho_boot_p025":float(np.quantile(vals,.025)),
     "rho_boot_p975":float(np.quantile(vals,.975)),
     "same_sign_fraction":float(np.mean(np.sign(vals)==np.sign(value))),
     "bootstrap_sd":float(np.std(vals,ddof=1))})
 out=root/"combined"/"section_block_resampling.csv"
 pd.DataFrame(rows).to_csv(out,index=False)
 summary={"dataset":dataset,"budget":budget,"n_sections":len(sections),"n_replicates":reps,
  "same_sign_fraction_ge_0_9":float(np.mean(pd.DataFrame(rows).same_sign_fraction>=.9)),
  "scope":"Reweights serial sections of fixed GP residuals; does not refit the GP or quantify registration and model-selection uncertainty"}
 (root/"combined"/"section_block_resampling_summary.json").write_text(json.dumps(summary,indent=2))
 print(summary,flush=True)
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--dataset",choices=["crc","bc"],required=True)
 p.add_argument("--budget",type=int,required=True);p.add_argument("--reps",type=int,default=100)
 p.add_argument('--root',type=Path,required=True)
 a=p.parse_args();run(a.dataset,a.budget,a.reps,a.root)
