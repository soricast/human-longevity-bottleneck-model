import pandas as pd,numpy as np,matplotlib.pyplot as plt
from pathlib import Path
from matplotlib import rcParams
rcParams.update({'font.family':'DejaVu Sans','font.size':10,'savefig.dpi':220})
O=Path('analysis_outputs');x=pd.read_csv(O/'replicate_paired.csv',keep_default_na=False);h=pd.read_csv(O/'age_hazard_raw.csv')
groups=['DNMT3A','TET2','ASXL1','JAK2'];colors=['#0072B2','#E69F00','#009E73','#CC79A7']

fig,ax=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
for j,(metric,label) in enumerate([('delta_P_heart_before_brain_matched','ΔP(heart before brain) versus matched null'),('delta_T50_matched','ΔT50 versus matched null (years)')]):
 for i,g in enumerate(groups):
  a=x[(x.scenario=='cardiac')&(x.group==g)][metric].to_numpy()
  ax[j].errorbar(i,a.mean(),yerr=[[a.mean()-a.min()],[a.max()-a.mean()]],fmt='o',color=colors[i],capsize=4,markersize=7)
 ax[j].axhline(0,color='0.5',lw=.8);ax[j].set_xticks(range(4),groups);ax[j].set_ylabel(label);ax[j].grid(axis='y',alpha=.2)
fig.suptitle('Paired simulations: 10 seeds × 100,000 individuals; bars show seed range')
fig.savefig(O/'Figure1_paired_effects.png');plt.close(fig)

fig,ax=plt.subplots(figsize=(7.4,4.6),layout='constrained')
for g,c in [('No CHIP','#333333'),*zip(groups,colors)]:
 a=h[h.group==g].groupby('age').ratio_heart_brain.agg(['mean','min','max'])
 ax.plot(a.index,a['mean'],label=g,color=c,marker='o',lw=2)
 ax.fill_between(a.index,a['min'],a['max'],color=c,alpha=.12)
ax.axhline(1,color='0.5',lw=.8,ls='--');ax.set(xlabel='Age (years)',ylabel='Mean modeled heart/brain hazard ratio',title='Association-calibrated hazard contrast at observable ages')
ax.legend(ncol=2,frameon=False);ax.grid(alpha=.2);fig.savefig(O/'Figure2_age_hazards.png');plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
sc=['null','cardiac','brain_half','brain_equal','brain_double'];labels=['Uncoupled','Cardiac only','Brain + half','Brain + equal','Brain + double'];xx=np.arange(len(sc))
for i,g in enumerate(['TET2','ASXL1','JAK2']):
 a=x[x.group==g].groupby('scenario').P_heart_before_brain.mean().reindex(sc)
 ax[0].plot(xx,a,marker=['o','s','^'][i],label=g,color=colors[i+1])
ax[0].axhline(.5,color='0.4',ls='--',lw=.8)
ax[0].set(ylabel='P(heart before brain)',title='Exploratory brain coupling')
for g,c in [('Global','#333333'),('TET2',colors[1])]:
 a=x[x.group==g].groupby('scenario').P_terminal_heart.mean().reindex(['null','cardiac','background_equal'])
 ax[1].plot(range(3),a,marker='o',label=g,color=c)
ax[1].set_xticks(range(3),['Uncoupled','Cardiac only','+ unrelated hazard'],rotation=15)
ax[1].set(ylabel='P(heart is terminal cause)',title='Exploratory other-cause hazard')
for a in ax:a.grid(axis='y',alpha=.2);a.legend(frameon=False)
ax[0].set_xticks(xx,labels,rotation=20)
fig.suptitle('Additional hazards are hypothetical multiples of the imposed cardiac effect')
fig.savefig(O/'Figure3_competing_hazards.png');plt.close(fig)

import openpyxl
w=openpyxl.load_workbook('upload/Supplementary_Data_1_Master_Parameters (1).xlsx',read_only=True,data_only=True)
sens=pd.DataFrame(list(w['Sensitivity'].values));sens.columns=sens.iloc[0];sens=sens.iloc[1:]
sc=['Frozen primary','mu_brain 95% CI low','mu_brain 95% CI high','mu_heart 95% CI low','mu_heart 95% CI high']
fig,ax=plt.subplots(figsize=(9.2,4.7),layout='constrained')
for i,g in enumerate(['No CHIP',*groups]):
 a=sens[(sens['group']==g)&sens['scenario'].isin(sc)].set_index('scenario').reindex(sc)
 ax.plot(range(5),a['P_heart_first'].astype(float),label=g,marker='o',color=(['#333333']+colors)[i])
ax.axhline(.5,color='0.4',ls='--',lw=.8);ax.set_xticks(range(5),['Central','Brain μ low','Brain μ high','Heart μ low','Heart μ high'],rotation=15)
ax.set(ylabel='P(heart before brain)',title='Intrinsic-hazard limits: pre-existing single-seed Weibull sensitivity')
ax.legend(ncol=3,frameon=False);ax.grid(axis='y',alpha=.2)
fig.savefig(O/'Figure4_intrinsic_sensitivity.png');plt.close(fig)
