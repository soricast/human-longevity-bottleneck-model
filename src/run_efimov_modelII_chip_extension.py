#!/usr/bin/env python3
"""
Independent reconstruction of Efimov et al. Model II (brain/heart)
plus the frozen CHIP/VAF cardiovascular extension.

Important reproducibility note
------------------------------
The public repository Config currently contains:
    lambda_bg = 0.0016133681587490935
but the published non-aging median of 1759 years implies:
    lambda_bg = ln(2)/1759 = 0.000394057521637...
The latter is required to reproduce the published ~194 y brain,
~208 y heart, and ~156 y independent brain-heart combined medians.

This script therefore treats the published article as authoritative for
lambda_bg and retains all other Model II parameters/logic from the public
Efimov repository.
"""
import numpy as np
import pandas as pd
import math
from scipy.stats import lognorm
from scipy.special import expit, logit

SEED = 20260819
N = 100_000
N_EFIMOV_MC = 100_000
DT = 0.25
TMAX = 1000.0
LAMBDA_BG = math.log(2)/1759.0

PARAMS = {
    "brain":{
        "SNV":(17.489,16.112,18.865),
        "indels":(6.926,5.484,8.367),
        "p":(12.34e-5,3.87e-5),
        "K":(3.5e9,7e8),
        "xcrit":0.60,
    },
    "heart":{
        "SNV":(36.369,19.519,53.218),
        "indels":(14.403,6.644,23.604),
        "p":(5.67e-5,3.87e-5),
        "K":(3.2e9,7.5e8),
        "xcrit":0.55,
    }
}

DRIVERS = np.array(["DNMT3A","TET2","ASXL1","JAK2"], dtype=object)
WEIGHTS = np.array([39,18,3,8],dtype=float)
WEIGHTS /= WEIGHTS.sum()
GROWTH = {
    "DNMT3A":(0.05,0.05),
    "TET2":(0.10,0.05),
    "ASXL1":(0.10,0.05),
    "JAK2":(0.10,0.10),
}
PRIMARY_HR = {"DNMT3A":1.15,"TET2":1.52,"ASXL1":1.52,"JAK2":1.52}
K_V = 1.0
VAF_REF = 10.0

def lognormal_params(mean,se):
    sigma=np.sqrt(np.log(1+(se/mean)**2))
    mu=np.log(mean)-0.5*sigma**2
    return mu,sigma

def efimov_model_ii_curve(organ,t,seed=123,n_mc=N_EFIMOV_MC):
    c=PARAMS[organ]
    snv_se=(c["SNV"][2]-c["SNV"][1])/(2*1.96)
    indel_se=(c["indels"][2]-c["indels"][1])/(2*1.96)
    mu_mean=c["SNV"][0]*c["p"][0]+c["indels"][0]*c["p"][1]
    mu_se=np.sqrt((snv_se*c["p"][0])**2+(indel_se*c["p"][1])**2)

    mu_ln,mu_sigma=lognormal_params(mu_mean,mu_se)
    K_ln,K_sigma=lognormal_params(*c["K"])
    rng=np.random.default_rng(seed)
    mu=rng.lognormal(mu_ln,mu_sigma,n_mc)

    xcrit=c["xcrit"]*c["K"][0]
    S=np.zeros(len(t))
    for start in range(0,n_mc,2500):
        m=mu[start:start+2500]
        thresholds=xcrit*np.exp(np.outer(m,t))
        S += np.sum(
            1-lognorm.cdf(thresholds,s=K_sigma,scale=np.exp(K_ln)),
            axis=0
        )
    return S/n_mc

def crossing(t,S,q=.5):
    i=np.where(S<=q)[0][0]
    return np.interp(q,[S[i],S[i-1]],[t[i],t[i-1]])

def main():
    t=np.arange(0,TMAX+DT,DT)
    Sb=efimov_model_ii_curve("brain",t)
    Sh=efimov_model_ii_curve("heart",t)
    Sbg=np.exp(-LAMBDA_BG*t)

    print("Published-backbone reproduction:")
    print("Brain T50:",crossing(t,Sb*Sbg))
    print("Heart T50:",crossing(t,Sh*Sbg))
    print("Independent brain-heart T50:",crossing(t,Sb*Sh*Sbg))

    # Synthetic CHIP population
    rng=np.random.default_rng(SEED)
    carrier=rng.random(N)<0.30
    driver=np.full(N,"No CHIP",dtype=object)
    driver[carrier]=rng.choice(DRIVERS,size=carrier.sum(),p=WEIGHTS)

    # Conditional age at crossing VAF>=2% between ages 20 and 80
    a_inc=-4.359484788519582
    b_inc=0.051171374892018846
    L0=expit(a_inc+b_inc*20)
    L1=expit(a_inc+b_inc*80)
    u=rng.random(carrier.sum())
    Lt=L0+u*(L1-L0)
    onset=np.full(N,np.nan)
    onset[carrier]=(logit(Lt)-a_inc)/b_inc

    growth=np.zeros(N)
    for g in DRIVERS:
        mu,sd=GROWTH[g]
        idx=np.where(driver==g)[0]
        vals=rng.normal(mu,sd,len(idx))
        bad=(vals<0.002)|(vals>0.50)
        while bad.any():
            vals[bad]=rng.normal(mu,sd,bad.sum())
            bad=(vals<0.002)|(vals>0.50)
        growth[idx]=vals

    Ebg=rng.exponential(size=N)
    Eb=rng.exponential(size=N)
    Eh=rng.exponential(size=N)

    Hb=-np.log(np.clip(Sb,1e-300,1))
    Hh=-np.log(np.clip(Sh,1e-300,1))
    dHh=np.diff(Hh)
    tm=(t[:-1]+t[1:])/2

    Tb=np.interp(Eb,Hb,t)
    Th=np.interp(Eh,Hh,t)
    Tbg=Ebg/LAMBDA_BG

    ratio=(0.50-0.02)/0.02
    for g in DRIVERS:
        ids=np.where(driver==g)[0]
        beta=math.log(PRIMARY_HR[g])/(VAF_REF/(K_V+VAF_REF))
        for start in range(0,len(ids),600):
            idx=ids[start:start+600]
            delta=tm[None,:]-onset[idx,None]
            rr=growth[idx,None]
            vaf=np.zeros(delta.shape,dtype=np.float32)
            mask=delta>=0
            vf=0.50/(1+ratio*np.exp(-rr*np.maximum(delta,0)))
            vaf[mask]=vf[mask]
            vp=100*vaf
            f=vp/(K_V+vp)
            multiplier=np.exp(beta*f)
            H=np.cumsum(multiplier*dHh[None,:],axis=1)
            crossed=H>=Eh[idx,None]
            has=crossed.any(axis=1)
            j=np.argmax(crossed,axis=1)
            tt=np.full(len(idx),t[-1]+DT)
            q=np.where(has)[0]
            ji=j[q]
            prev=np.where(ji==0,0,H[q,np.maximum(ji-1,0)])
            curr=H[q,ji]
            tt0=t[ji]; tt1=t[ji+1]
            tt[q]=tt0+(Eh[idx[q]]-prev)*(tt1-tt0)/(curr-prev)
            Th[idx]=tt

    life=np.minimum.reduce([Tbg,Tb,Th])
    rows=[]
    for g in ["No CHIP"]+DRIVERS.tolist()+["Global"]:
        idx=np.ones(N,dtype=bool) if g=="Global" else (driver==g)
        rows.append({
            "group":g,
            "n":int(idx.sum()),
            "T50_years":float(np.median(life[idx])),
            "P_heart_first":float(np.mean(Th[idx]<Tb[idx])),
            "P_brain_first":float(np.mean(Tb[idx]<Th[idx])),
        })
    df=pd.DataFrame(rows)
    df.to_csv("Efimov_ModelII_CHIP_primary_results.csv",index=False)
    print(df.to_string(index=False))

if __name__=="__main__":
    main()
