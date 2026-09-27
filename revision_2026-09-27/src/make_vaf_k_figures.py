import numpy as np,pandas as pd,matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path
from scipy.special import expit,logit
rcParams.update({'font.family':'DejaVu Sans','font.size':10,'savefig.dpi':220})
O=Path('analysis_outputs');O.mkdir(exist_ok=True)
N=100000;rng=np.random.default_rng(20260819)
drivers=np.array(['DNMT3A','TET2','ASXL1','JAK2']);weights=np.array([39,18,3,8],float);weights/=weights.sum()
carrier=rng.random(N)<.30;d=rng.choice(drivers,size=carrier.sum(),p=weights)
a=-4.359484788519582;b=.051171374892018846
L0=expit(a+b*20);L1=expit(a+b*80)
onset=(logit(L0+rng.random(carrier.sum())*(L1-L0))-a)/b
G={'DNMT3A':(.05,.05),'TET2':(.10,.05),'ASXL1':(.10,.05),'JAK2':(.10,.10)}
growth=np.zeros(len(d))
for g in drivers:
 mu,sd=G[g];ix=np.where(d==g)[0];v=rng.normal(mu,sd,len(ix));bad=(v<.002)|(v>.5)
 while bad.any():v[bad]=rng.normal(mu,sd,bad.sum());bad=(v<.002)|(v>.5)
 growth[ix]=v
ages=np.arange(20,111,2);rows=[];colors=['#0072B2','#E69F00','#009E73','#CC79A7']
fig,axs=plt.subplots(2,2,figsize=(10.5,7),layout='constrained',sharex=True,sharey=True)
for ax,g,c in zip(axs.flat,drivers,colors):
 ix=np.where(d==g)[0];o=onset[ix,None];r=growth[ix,None]
 delta=ages[None,:]-o
 vaf=np.where(delta>=0,50/(1+24*np.exp(-r*np.maximum(delta,0))),0)
 q=np.quantile(vaf,[.1,.5,.9],axis=0)
 for j,age in enumerate(ages):rows.append(dict(driver=g,n=len(ix),age=age,q10_VAF_pct=q[0,j],median_VAF_pct=q[1,j],q90_VAF_pct=q[2,j],mean_VAF_pct=np.mean(vaf[:,j]),fraction_above_2pct=np.mean(vaf[:,j]>=2)))
 ax.fill_between(ages,q[0],q[2],color=c,alpha=.22,label='10th–90th percentile')
 ax.plot(ages,q[1],color=c,lw=2,label='Population median')
 # Central-rate representative uses median onset; JAK2 differs in dispersion, not central growth.
 oo=np.median(onset[ix]);rr=G[g][0]
 rep=np.where(ages>=oo,50/(1+24*np.exp(-rr*np.maximum(ages-oo,0))),0)
 ax.plot(ages,rep,color=c,lw=1.4,ls='--',label='Representative central rate')
 ax.axhline(2,color='0.45',lw=.8,ls=':')
 ax.set_title(f'{g} (n={len(ix):,})');ax.grid(alpha=.15)
 ax.set(xlim=(20,110),ylim=(0,50))
for ax in axs[1]:ax.set_xlabel('Age (years)')
for ax in axs[:,0]:ax.set_ylabel('Modeled VAF (%)')
axs[0,0].legend(frameon=False,loc='upper left',fontsize=8)
fig.suptitle('Assigned-carrier VAF trajectories: population spread and representative curve')
fig.savefig(O/'Figure6_VAF_population.png');plt.close(fig)
pd.DataFrame(rows).to_csv(O/'VAF_population_quantiles.csv',index=False)

# Every curve is exactly anchored at the same HR at VAF=10%, so K only changes shape.
v=np.linspace(0,30,301);HR=1.52;kv_values=[.5,1,2,5,10]
curve=[]
fig,ax=plt.subplots(figsize=(7.6,4.6),layout='constrained')
for k,c in zip(kv_values,['#0072B2','#E69F00','#009E73','#CC79A7','#D55E00']):
 f=v/(k+v);f10=10/(k+10)
 multiplier=np.exp(np.log(HR)*f/f10)
 ax.plot(v,multiplier,label=f'K_V={k:g} pp',color=c,lw=2)
 for vi,mi in zip(v,multiplier):curve.append(dict(K_V_pp=k,VAF_pct=vi,HR_multiplier=mi))
ax.axvline(10,color='0.5',ls='--',lw=.8);ax.axhline(HR,color='0.5',ls='--',lw=.8)
ax.set(xlabel='VAF (%)',ylabel='Modeled cardiac hazard multiplier',title='K_V changes extrapolation around the shared 10% VAF anchor')
ax.legend(frameon=False,ncol=2);ax.grid(alpha=.15)
fig.savefig(O/'Figure7_KV_extrapolation.png');plt.close(fig)
pd.DataFrame(curve).to_csv(O/'KV_extrapolation_curves.csv',index=False)
