import math,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.special import expit,logit
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_replicates import N,t,tm,dHh,Hh,Hb,LAMBDA_BG,DT,DRIVERS,GROWTH,PRIMARY_HR,K_V,VAF_REF,invert

OUT=Path('analysis_outputs')
ARIC=np.array([540,258,97,18],float);ARIC/=ARIC.sum()
LOTHIAN=np.array([39,18,3,8],float);LOTHIAN/=LOTHIAN.sum()
P_SECOND=74/457 # ARIC baseline carriers with >1 clone, 383/457 had one.

def draw_clone(rng,carrier,weights):
 d=np.full(N,'No CHIP',dtype=object);d[carrier]=rng.choice(DRIVERS,size=carrier.sum(),p=weights)
 onset=np.full(N,np.nan);a=-4.359484788519582;b=.051171374892018846
 lo=expit(a+b*20);hi=expit(a+b*80)
 onset[carrier]=(logit(lo+rng.random(carrier.sum())*(hi-lo))-a)/b
 growth=np.zeros(N)
 for g in DRIVERS:
  mu,sd=GROWTH[g];ix=np.where(d==g)[0];v=rng.normal(mu,sd,len(ix));bad=(v<.002)|(v>.5)
  while bad.any():v[bad]=rng.normal(mu,sd,bad.sum());bad=(v<.002)|(v>.5)
  growth[ix]=v
 return d,onset,growth

def log_mult(driver,onset,growth,ix):
 out=np.zeros((len(ix),len(tm)),dtype=np.float32)
 for g in DRIVERS:
  loc=np.where(driver[ix]==g)[0]
  if not len(loc):continue
  ids=ix[loc]
  delta=tm[None,:]-onset[ids,None]
  vp=np.where(delta>=0,50/(1+24*np.exp(-growth[ids,None]*np.maximum(delta,0))),0)
  beta=math.log(PRIMARY_HR[g])/(VAF_REF/(K_V+VAF_REF))
  out[loc]=beta*vp/(K_V+vp)
 return out

def run(seed,weights,label):
 rng=np.random.default_rng(seed);carrier=rng.random(N)<.30
 first,onset1,growth1=draw_clone(rng,carrier,weights)
 # Draw all second-clone variables before failure thresholds, shared in max/product scenarios.
 second_mask=carrier & (rng.random(N)<P_SECOND)
 second,onset2,growth2=draw_clone(rng,second_mask,weights)
 Ebg=rng.exponential(size=N);Eb=rng.exponential(size=N);Eh=rng.exponential(size=N)
 tb=invert(Eb,Hb);th0=invert(Eh,Hh);bg=Ebg/LAMBDA_BG
 null_life=np.minimum.reduce([tb,th0,bg]);pnull=np.mean(th0<tb)
 th={key:th0.copy() for key in ['single','max','product']}
 ids=np.where(carrier)[0]
 for st in range(0,len(ids),350):
  ix=ids[st:st+350]
  l1=log_mult(first,onset1,growth1,ix);l2=log_mult(second,onset2,growth2,ix)
  for name,lm in [('single',l1),('max',np.maximum(l1,l2)),('product',l1+l2)]:
   h=np.cumsum(np.exp(lm)*dHh[None,:],axis=1)
   for k,ii in enumerate(ix):th[name][ii]=np.interp(Eh[ii],np.r_[0,h[k]],t,right=t[-1]+DT)
 rows=[]
 for name,ht in th.items():
  life=np.minimum.reduce([tb,ht,bg]);p=np.mean(ht<tb)
  rows.append(dict(seed=seed,composition=label,multiple_rule=name,second_clone_carrier_fraction=second_mask.sum()/carrier.sum(),T50=np.median(life),delta_T50=np.median(life)-np.median(null_life),P_heart_before_brain=p,delta_P=p-pnull))
 return rows

if __name__=='__main__':
 rows=[]
 for label,w in [('Lothian',LOTHIAN),('ARIC_four_gene',ARIC)]:
  for seed in range(20260819,20260829):
   rows+=run(seed,w,label);print(label,seed,flush=True)
 df=pd.DataFrame(rows);df.to_csv(OUT/'composition_multiclone_raw.csv',index=False)
 agg=df.groupby(['composition','multiple_rule']).agg({c:['mean','min','max'] for c in ['T50','delta_T50','P_heart_before_brain','delta_P','second_clone_carrier_fraction']})
 agg.columns=[a+'_'+b for a,b in agg.columns]
 agg.reset_index().to_csv(OUT/'composition_multiclone_summary.csv',index=False)
