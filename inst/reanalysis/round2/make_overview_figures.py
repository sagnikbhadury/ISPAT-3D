#!/usr/bin/env python
"""Revised vector overview figures for the final analytical pipeline."""
from pathlib import Path
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
root=Path(sys.argv[1]);root.mkdir(parents=True,exist_ok=True)
def box(ax,x,y,w,h,title,body,face):
 p=FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.02,rounding_size=0.03",
  linewidth=1.2,edgecolor="#42546a",facecolor=face)
 ax.add_patch(p)
 ax.text(x+w/2,y+h*.68,title,ha="center",va="center",fontsize=11,weight="bold")
 ax.text(x+w/2,y+h*.29,body,ha="center",va="center",fontsize=9,color="#263845")
def arrow(ax,a,b,y):
 ax.add_patch(FancyArrowPatch((a,y),(b,y),arrowstyle="-|>",mutation_scale=18,
  linewidth=1.5,color="#4f6575"))
fig,ax=plt.subplots(figsize=(12.5,3.2));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
items=[
 ("Registered volume","Cell labels and\nX Y Z coordinates","#e9f2f8"),
 ("Relative tumor zones","Five density strata\nwithin each specimen","#e9f2f8"),
 ("Spatial adjustment","Smooth 3D trends\nestimated by GP","#e8f3eb"),
 ("Covariance model","Shared and zone\nvarying covariance","#f4eedf"),
 ("Conditional graphs","Signed partial correlations\nand tier stability","#f7e9e9")]
w=.17;x=[.02,.22,.42,.62,.82]
for j,(title,body,color) in enumerate(items):
 box(ax,x[j],.26,w,.54,title,body,color)
 if j<4:arrow(ax,x[j]+w,x[j+1],.53)
ax.text(.5,.08,"Edges describe conditional associations among measured cell-density variables; their biological cause remains unresolved.",
 ha="center",fontsize=9,color="#42546a")
fig.savefig(root/"figure01_conceptual_workflow.pdf",bbox_inches="tight");plt.close(fig)
fig,ax=plt.subplots(figsize=(12.5,5.3));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
items=[
 (.03,.68,.27,.22,"Input","Registered serial sections\ncell-type KDE values and coordinates","#e9f2f8"),
 (.365,.68,.27,.22,"Zone assignment","Five mutually exclusive relative\ntumor-density strata","#e9f2f8"),
 (.70,.68,.27,.22,"Selected cells","CRC three per-zone budgets\nbreast three total budgets","#e9f2f8"),
 (.03,.37,.27,.22,"Anisotropic GP","Matérn 3/2 with separate\nX Y Z ranges","#e8f3eb"),
 (.365,.37,.27,.22,"Scalable fit","Spatially balanced anchors\nVecchia 15-neighbor likelihood","#e8f3eb"),
 (.70,.37,.27,.22,"Adjusted densities","Predict smooth field on selected\ncells then subtract it","#e8f3eb"),
 (.03,.06,.27,.22,"Covariance fit","Shared plus zone factors and\ndiagonal variance by Gaussian ML","#f4eedf"),
 (.365,.06,.27,.22,"Zone graphs","Invert complete zone covariance\nto obtain partial correlations","#f4eedf"),
 (.70,.06,.27,.22,"Checks","Section-wise planar GP\nEBIC graphical lasso and tiers","#f7e9e9")]
for x,y,w,h,t,b,c in items:box(ax,x,y,w,h,t,b,c)
for y in (.79,.48,.17):
 arrow(ax,.30,.365,y);arrow(ax,.635,.70,y)
ax.add_patch(FancyArrowPatch((.835,.68),(.165,.59),connectionstyle="arc3,rad=-.24",
 arrowstyle="-|>",mutation_scale=15,lw=1.3,color="#748596"))
ax.add_patch(FancyArrowPatch((.835,.37),(.165,.28),connectionstyle="arc3,rad=-.24",
 arrowstyle="-|>",mutation_scale=15,lw=1.3,color="#748596"))
fig.savefig(root/"figure02_pipeline.pdf",bbox_inches="tight");plt.close(fig)
