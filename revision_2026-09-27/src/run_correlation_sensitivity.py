import math
import numpy as np,pandas as pd
from run_joint_uncertainty import model,MU_B,MU_H,HR_D,HR_O,sigma
from pathlib import Path
rng=np.random.default_rng(20260928);rows=[]
for rho in [-.5,.5]:
 for j in range(60):
  z1,z2=rng.normal(size=2)
  fb=np.exp(sigma(MU_B)*z1)
  fh=np.exp(sigma(MU_H)*(rho*z1+np.sqrt(1-rho*rho)*z2))
  hrd=np.exp(rng.normal(math.log(HR_D[0]),sigma(HR_D)))
  hro=np.exp(rng.normal(math.log(HR_O[0]),sigma(HR_O)))
  rows.append(dict(mu_log_correlation=rho,draw=j,**model(fb,fh,hrd,hro)))
 if j%20==0:print(rho,j,flush=True)
o=Path('analysis_outputs');pd.DataFrame(rows).to_csv(o/'correlation_sensitivity_raw.csv',index=False)
a=pd.DataFrame(rows).groupby('mu_log_correlation').agg({c:['median',lambda x:np.quantile(x,.025),lambda x:np.quantile(x,.975)] for c in ['delta_P','delta_T50','primary_P','primary_T50']})
a.columns=[p+'_'+q.replace('<lambda_0>','q025').replace('<lambda_1>','q975') for p,q in a.columns]
a.reset_index().to_csv(o/'correlation_sensitivity_summary.csv',index=False)
