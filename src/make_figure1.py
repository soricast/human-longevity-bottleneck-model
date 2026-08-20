#!/usr/bin/env python3
"""
Generate the single five-panel manuscript figure from frozen article outputs.

This script intentionally uses only the frozen numerical result tables stored
in ../results/ plus the frozen central VAF-growth assumptions in
../config/frozen_config.json.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
R = ROOT / "results"
C = ROOT / "config" / "frozen_config.json"
OUT = ROOT / "figures" / "Figure1_integrated_exact.png"

cfg = json.loads(C.read_text())
primary = pd.read_csv(R / "Efimov_ModelII_CHIP_primary_results.csv")
rep = pd.read_csv(R / "Efimov_literal_reproduction_check.csv")
m3 = pd.read_csv(R / "ModelIII_MonteCarlo_check.csv")
switch = pd.read_csv(R / "Efimov_ModelII_switching_sweep.csv")
validation = pd.read_csv(R / "Final_ModelII_primary_and_validation_results.csv")

drivers = ["DNMT3A","TET2","ASXL1","JAK2"]
groups = ["No CHIP"] + drivers
ages = np.linspace(20,150,400)
V0, VMAX, onset = 0.02, 0.50, 50.0
ratio = (VMAX-V0)/V0
growth = {k:v[0] for k,v in cfg["growth_params"].items()}
vaf={}
for g in drivers:
    d=np.maximum(ages-onset,0)
    z=VMAX/(1+ratio*np.exp(-growth[g]*d))
    z[ages<onset]=0
    vaf[g]=100*z

blue="#2F6F9F"; orange="#E89A3C"; green="#6AA56A"; red="#C95B55"; dark="#2B2B2B"
plt.rcParams.update({
    "font.family":"DejaVu Sans","font.size":9,"axes.labelsize":9,
    "xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8,
    "axes.linewidth":0.9,"axes.spines.top":False,"axes.spines.right":False
})

fig=plt.figure(figsize=(11.5,8.0))
gs=GridSpec(2,6,figure=fig,height_ratios=[1,1],hspace=0.48,wspace=0.65)

# a
ax=fig.add_subplot(gs[0,0:2])
labels=["Brain","Heart","Combined","Lung","Liver"]
vals=[
    float(rep.iloc[0]["literal_reconstruction_years"]),
    float(rep.iloc[1]["literal_reconstruction_years"]),
    float(rep.iloc[2]["literal_reconstruction_years"]),
    float(m3.loc[m3.model=="lungs","observed_death_median"].iloc[0]),
    float(m3.loc[m3.model=="liver","observed_death_median"].iloc[0]),
]
bars=ax.bar(labels,vals,color=orange,edgecolor=dark,linewidth=.3)
ax.set_yscale("log"); ax.set_ylim(80,1.1e5)
ax.set_ylabel("Median model time to failure (years)")
for b,x in zip(bars,vals):
    ax.text(b.get_x()+b.get_width()/2,x*1.12,f"{x:,.0f}",ha="center",fontsize=8)
ax.text(-.20,1.08,"a",transform=ax.transAxes,fontsize=15,fontweight="bold")

# b
ax=fig.add_subplot(gs[0,2:4])
for g,c in zip(drivers,[blue,orange,green,red]):
    ax.plot(ages,vaf[g],label=g,color=c,linewidth=1.8)
ax.axhline(2,color=blue,ls="--",lw=1)
ax.text(21,3.2,"2% CHIP\nthreshold",color=blue,fontsize=8)
ax.set(xlabel="Age (years)",ylabel="Synthetic VAF (%)",xlim=(20,150),ylim=(0,52))
ax.legend(frameon=False,loc="upper left")
ax.text(-.20,1.08,"b",transform=ax.transAxes,fontsize=15,fontweight="bold")

# c
ax=fig.add_subplot(gs[0,4:6])
ph=[float(primary.loc[primary.group==g,"P_heart_first"].iloc[0]) for g in groups]
t50=[float(primary.loc[primary.group==g,"T50_years"].iloc[0]) for g in groups]
bars=ax.bar(groups,ph,color=orange,edgecolor=dark,linewidth=.3)
ax.axhline(.5,color=blue,ls="--",lw=1)
ax.set_ylabel("P(Heart-first)"); ax.set_ylim(.40,.60)
for b,x in zip(bars,t50):
    ax.text(b.get_x()+b.get_width()/2,b.get_height()+.006,f"{x:.1f} y",ha="center",fontsize=8)
ax.text(-.20,1.08,"c",transform=ax.transAxes,fontsize=15,fontweight="bold")

# d
ax=fig.add_subplot(gs[1,0:3])
ax.plot(switch["HR_at_10pct_VAF"],switch["P_heart_first_CHIP"],marker="o",color=blue,lw=1.8,ms=4.5)
ax.axhline(.5,color=blue,ls="--",lw=1)
ax.axvline(cfg["estimated_switch_HR_among_CHIP"],color=blue,ls=":",lw=1)
ax.text(1.34,.505,"Switch ≈ 1.27",fontsize=8)
ax.set(xlabel="Common cardiovascular HR at 10% VAF",
       ylabel="P(Heart-first) among CHIP carriers",ylim=(.39,.65))
ax.text(-.12,1.06,"d",transform=ax.transAxes,fontsize=15,fontweight="bold")

# e
ax=fig.add_subplot(gs[1,3:6])
colors={"Primary incident HF":blue,"Positive CHD":orange,"CVD-null":green}
for scen in ["Primary incident HF","Positive CHD","CVD-null"]:
    d=validation[validation.scenario==scen]
    ys=[float(d.loc[d.group==g,"P_heart_first"].iloc[0]) for g in groups]
    ax.plot(groups,ys,marker="o",label=scen,color=colors[scen],lw=1.8,ms=4.5)
ax.axhline(.5,color=blue,ls="--",lw=1)
ax.set_ylabel("P(Heart-first)"); ax.set_ylim(.40,.90)
ax.legend(frameon=False,loc="upper left")
ax.text(-.12,1.06,"e",transform=ax.transAxes,fontsize=15,fontweight="bold")

fig.savefig(OUT,dpi=600,bbox_inches="tight")
print(OUT)
