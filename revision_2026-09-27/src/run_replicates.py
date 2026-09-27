import sys, json, math, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expit,logit
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_efimov_modelII_chip_extension import efimov_model_ii_curve, PARAMS, DRIVERS, WEIGHTS, GROWTH, PRIMARY_HR, LAMBDA_BG, K_V, VAF_REF

OUT=Path('analysis_outputs');OUT.mkdir(exist_ok=True)
N=100000;DT=.25;TMAX=1000
t=np.arange(0,TMAX+DT,DT);tm=(t[:-1]+t[1:])/2
cache=OUT/'intrinsic_curves_1000.npz'
if cache.exists():
 z=np.load(cache);Sb=z['Sb'];Sh=z['Sh']
else:
 Sb=efimov_model_ii_curve('brain',t);Sh=efimov_model_ii_curve('heart',t)
 np.savez(cache,Sb=Sb,Sh=Sh)
Hb=-np.log(np.clip(Sb,1e-300,1));Hh=-np.log(np.clip(Sh,1e-300,1))
dHb=np.diff(Hb);dHh=np.diff(Hh)
def invert(e,H):
 return np.interp(e,H,t,left=0,right=TMAX+DT)

def sample(seed):
 rng=np.random.default_rng(seed);carrier=rng.random(N)<.30
 driver=np.full(N,'No CHIP',dtype=object);driver[carrier]=rng.choice(DRIVERS,size=carrier.sum(),p=WEIGHTS)
 onset=np.full(N,np.nan)
 a=-4.359484788519582;b=.051171374892018846
 L0=expit(a+b*20);L1=expit(a+b*80)
 onset[carrier]=(logit(L0+rng.random(carrier.sum())*(L1-L0))-a)/b
 growth=np.zeros(N)
 for g in DRIVERS:
  mu,sd=GROWTH[g];ix=np.where(driver==g)[0];v=rng.normal(mu,sd,len(ix));bad=(v<.002)|(v>.5)
  while bad.any():v[bad]=rng.normal(mu,sd,bad.sum());bad=(v<.002)|(v>.5)
  growth[ix]=v
 return driver,onset,growth,rng.exponential(size=N),rng.exponential(size=N),rng.exponential(size=N)

def run(seed):
 driver,onset,growth,Ebg,Eb,Eh=sample(seed)
 # Same person, clone and three exponential thresholds in every scenario.
 base_b=invert(Eb,Hb);base_h=invert(Eh,Hh);base_bg=Ebg/LAMBDA_BG
 scenarios=[('null',0,0),('cardiac',1,0),('brain_half',1,.5),('brain_equal',1,1),('brain_double',1,2)]
 times={name:[base_b.copy(),base_h.copy(),base_bg.copy()] for name,_,_ in scenarios}
 # Independent extra background-hazard stress test, anchored to cardiac excess hazard.
 times['background_equal']=[base_b.copy(),base_h.copy(),base_bg.copy()]
 ratios=[]
 for g in DRIVERS:
  ids=np.where(driver==g)[0]
  beta=math.log(PRIMARY_HR[g])/(VAF_REF/(K_V+VAF_REF))
  for start in range(0,len(ids),500):
   ix=ids[start:start+500]
   delta=tm[None,:]-onset[ix,None]
   vf=.5/(1+24*np.exp(-growth[ix,None]*np.maximum(delta,0)))
   vp=np.where(delta>=0,100*vf,0)
   logm=beta*vp/(K_V+vp)
   m=np.exp(logm)
   ch=np.cumsum(m*dHh[None,:],axis=1)
   # Initial zero at t=0, then invert row-specific cumulative hazard.
   hh=np.concatenate((np.zeros((len(ix),1)),ch),axis=1)
   for k in range(len(ix)):
    times['cardiac'][1][ix[k]]=np.interp(Eh[ix[k]],hh[k],t,right=TMAX+DT)
   times['background_equal'][1][ix]=times['cardiac'][1][ix]
   for name,_,alpha in scenarios[2:]:
    bh=np.cumsum(np.exp(alpha*logm)*dHb[None,:],axis=1)
    for k in range(len(ix)):
     times[name][0][ix[k]]=np.interp(Eb[ix[k]],np.r_[0,bh[k]],t,right=TMAX+DT)
    times[name][1][ix]=times['cardiac'][1][ix]
   # Extra unrelated competing hazard, cumulative excess equal to the cardiac excess.
   bgh=np.cumsum(LAMBDA_BG*DT+(m-1)*dHh[None,:],axis=1)
   for k in range(len(ix)):
    times['background_equal'][2][ix[k]]=np.interp(Ebg[ix[k]],np.r_[0,bgh[k]],t,right=TMAX+DT)
  # Representative age-specific hazard ratio, weighted over clone parameter draws.
  for age in [60,70,80,90,100]:
   vi=.5/(1+24*np.exp(-growth[ids]*np.maximum(age-onset[ids],0)))
   vp=np.where(age>=onset[ids],100*vi,0)
   mult=np.exp(beta*vp/(K_V+vp))
   j=round(age/DT)
   ratios.append(dict(seed=seed,group=g,age=age,ratio_heart_brain=(dHh[j]/dHb[j])*np.mean(mult),mean_multiplier=np.mean(mult),q025_multiplier=np.quantile(mult,.025),q975_multiplier=np.quantile(mult,.975)))
 for age in [60,70,80,90,100]:
  j=round(age/DT)
  ratios.append(dict(seed=seed,group='No CHIP',age=age,ratio_heart_brain=dHh[j]/dHb[j],mean_multiplier=1,q025_multiplier=1,q975_multiplier=1))
 rows=[]
 for name,(tb,th,tbg) in times.items():
  life=np.minimum.reduce([tb,th,tbg]);cause=np.argmin(np.stack([tb,th,tbg]),axis=0)
  for group in ['Global','No CHIP',*DRIVERS]:
   ix=np.ones(N,dtype=bool) if group=='Global' else driver==group
   rows.append(dict(seed=seed,scenario=name,group=group,n=int(ix.sum()),T50=float(np.median(life[ix])),P_heart_before_brain=float(np.mean(th[ix]<tb[ix])),P_terminal_heart=float(np.mean(cause[ix]==1)),P_terminal_brain=float(np.mean(cause[ix]==0)),P_terminal_background=float(np.mean(cause[ix]==2))))
 return rows,ratios

if __name__=='__main__':
 seeds=[20260819+i for i in range(10)]
 allrows=[];allrat=[]
 for seed in seeds:
  now=time.time();r,h=run(seed);allrows+=r;allrat+=h
  pd.DataFrame(allrows).to_csv(OUT/'replicate_raw.csv',index=False)
  pd.DataFrame(allrat).to_csv(OUT/'age_hazard_raw.csv',index=False)
  print(seed,'seconds',round(time.time()-now,1),flush=True)
 df=pd.DataFrame(allrows); null=df[df.scenario=='null'].set_index(['seed','group'])
 for col in ['T50','P_heart_before_brain','P_terminal_heart']:
  df['delta_'+col+'_matched']=df.apply(lambda x:x[col]-null.loc[(x.seed,x.group),col],axis=1)
 df.to_csv(OUT/'replicate_paired.csv',index=False)
 agg=df.groupby(['scenario','group'],sort=False).agg({c:['mean',lambda x:np.quantile(x,.025),lambda x:np.quantile(x,.975)] for c in ['T50','P_heart_before_brain','P_terminal_heart','delta_T50_matched','delta_P_heart_before_brain_matched','delta_P_terminal_heart_matched']})
 agg.columns=[a+'_'+b.replace('<lambda_0>','lo').replace('<lambda_1>','hi') for a,b in agg.columns]
 agg.reset_index().to_csv(OUT/'replicate_summary.csv',index=False)
 hr=pd.DataFrame(allrat).groupby(['group','age'],sort=False).agg({'ratio_heart_brain':['mean','min','max'],'mean_multiplier':['mean','min','max']})
 hr.columns=[a+'_'+b for a,b in hr.columns];hr.reset_index().to_csv(OUT/'age_hazard_summary.csv',index=False)
