from pathlib import Path
import pandas as pd,numpy as np,matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'font.family':'DejaVu Sans','font.size':10,'savefig.dpi':220})
o=Path('analysis_outputs');d=pd.read_csv(o/'joint_uncertainty_draws.csv')
fig,ax=plt.subplots(1,3,figsize=(12,4),layout='constrained')
colors=['#0072B2','#E69F00']
for i,(scenario,label) in enumerate([('persistent','Persistent HR'),('stop_at_100','HR ends at 100')]):
 g=d[d.scenario==scenario]
 for a,col,title in zip(ax,['primary_P','delta_P','delta_T50'],['Absolute P(heart before brain)','Paired ΔP','Paired ΔT50 (years)']):
  v=g[col].to_numpy();q=np.quantile(v,[.025,.5,.975])
  a.errorbar(i,q[1],yerr=[[q[1]-q[0]],[q[2]-q[1]]],fmt='o',color=colors[i],markersize=7,capsize=4)
  a.set_title(title);a.grid(axis='y',alpha=.2)
for a in ax:a.set_xticks([0,1],['Persistent','End at 100']);a.tick_params(axis='x',rotation=15)
ax[0].axhline(.5,ls='--',color='0.5',lw=.8)
fig.suptitle('Joint parameter draws (200): median and central 95% simulation quantiles')
fig.savefig(o/'Figure5_joint_uncertainty.png');plt.close(fig)
