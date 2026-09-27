"""Conditional joint-parameter uncertainty, separate from seed variation.

The organ mutation-rate parameters scale time in the original Model II
survival curves because Xcrit*exp(mu*t) depends on mu*t. Hazard ratios
are sampled on the log scale from published point estimates and 95% CIs.
"""
import math
from pathlib import Path
import numpy as np,pandas as pd
from scipy.special import expit,logit
from run_replicates import Hb,Hh,LAMBDA_BG,GROWTH,DRIVERS,WEIGHTS,K_V,VAF_REF

OUT=Path('analysis_outputs');OUT.mkdir(exist_ok=True)
N=50000;STEP=.5;T=np.arange(0,1000+STEP,STEP);MID=(T[:-1]+T[1:])/2
BASE_T=np.arange(0,1000.25,.25)
MU_B=(2.44,1.53,3.95);MU_H=(2.61,1.41,4.62)
HR_D=(1.15,1.00,1.31);HR_O=(1.52,1.33,1.75)
def sigma(triple):return (math.log(triple[2])-math.log(triple[1]))/(2*1.96)

def population(seed=20260819):
 rng=np.random.default_rng(seed);carrier=rng.random(N)<.30
 drivers=np.full(N,'No CHIP',dtype=object);drivers[carrier]=rng.choice(DRIVERS,size=carrier.sum(),p=WEIGHTS)
 onset=np.full(N,np.nan);a=-4.359484788519582;b=.051171374892018846
 lo=expit(a+b*20);hi=expit(a+b*80)
 onset[carrier]=(logit(lo+rng.random(carrier.sum())*(hi-lo))-a)/b
 growth=np.zeros(N)
 for g in DRIVERS:
  mu,sd=GROWTH[g];ix=np.where(drivers==g)[0];v=rng.normal(mu,sd,len(ix));bad=(v<.002)|(v>.5)
  while bad.any():v[bad]=rng.normal(mu,sd,bad.sum());bad=(v<.002)|(v>.5)
  growth[ix]=v
 ix=np.where(carrier)[0]
 delta=MID[None,:]-onset[ix,None]
 vp=np.where(delta>=0,50/(1+24*np.exp(-growth[ix,None]*np.maximum(delta,0))),0)
 F=(vp/(K_V+vp)).astype(np.float32)
 return ix,drivers[ix],F,rng.exponential(size=N),rng.exponential(size=N),rng.exponential(size=N)

ids,drivers,F,Ebg,Eb,Eh=population()
bg=Ebg/LAMBDA_BG
def model(fb,fh,hrd,hro,limit_age=None):
 hb=np.interp(T*fb,BASE_T,Hb,right=Hb[-1]);hh=np.interp(T*fh,BASE_T,Hh,right=Hh[-1])
 tb=np.interp(Eb,hb,T,right=T[-1]+STEP)
 th0=np.interp(Eh,hh,T,right=T[-1]+STEP)
 life0=np.minimum.reduce([tb,th0,bg]);t0=float(np.median(life0));p0=float(np.mean(th0<tb))
 th=th0.copy();dH=np.diff(hh)
 betas=np.where(drivers=='DNMT3A',math.log(hrd),math.log(hro))/(VAF_REF/(K_V+VAF_REF))
 for st in range(0,len(ids),400):
  sub=ids[st:st+400];lm=betas[st:st+400,None]*F[st:st+400]
  if limit_age is not None:lm=lm.copy();lm[:,MID>=limit_age]=0
  h=np.cumsum(np.exp(lm)*dH[None,:],axis=1)
  for k,i in enumerate(sub):th[i]=np.interp(Eh[i],np.r_[0,h[k]],T,right=T[-1]+STEP)
 life=np.minimum.reduce([tb,th,bg])
 return dict(null_T50=t0,primary_T50=float(np.median(life)),delta_T50=float(np.median(life)-t0),null_P=p0,primary_P=float(np.mean(th<tb)),delta_P=float(np.mean(th<tb)-p0))

if __name__=='__main__':
 central=model(1,1,HR_D[0],HR_O[0]);central_limited=model(1,1,HR_D[0],HR_O[0],100)
 print('central',central,'age<=100',central_limited,flush=True)
 rng=np.random.default_rng(20260927);rows=[]
 for j in range(200):
  fb=float(np.exp(rng.normal(0,sigma(MU_B))));fh=float(np.exp(rng.normal(0,sigma(MU_H))))
  hrd=float(np.exp(rng.normal(math.log(HR_D[0]),sigma(HR_D))))
  hro=float(np.exp(rng.normal(math.log(HR_O[0]),sigma(HR_O))))
  for age,name in [(None,'persistent'),(100,'stop_at_100')]:
   rows.append(dict(draw=j,scenario=name,mu_brain=MU_B[0]*fb,mu_heart=MU_H[0]*fh,HR_DNMT3A=hrd,HR_nonDNMT3A=hro,**model(fb,fh,hrd,hro,age)))
  if j%20==0:print(j,flush=True)
 pd.DataFrame(rows).to_csv(OUT/'joint_uncertainty_draws.csv',index=False)
 summ=[]
 for name,g in pd.DataFrame(rows).groupby('scenario'):
  for col in ['null_T50','primary_T50','delta_T50','null_P','primary_P','delta_P']:
   v=g[col].to_numpy()
   summ.append(dict(scenario=name,quantity=col,central=(central if name=='persistent' else central_limited)[col],median=np.median(v),q025=np.quantile(v,.025),q975=np.quantile(v,.975),min=v.min(),max=v.max(),Pr_positive=np.mean(v>0)))
 pd.DataFrame(summ).to_csv(OUT/'joint_uncertainty_summary.csv',index=False)
