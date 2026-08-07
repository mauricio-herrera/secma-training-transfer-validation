"""Reproduce the observed-data statistics and posterior-predictive design analysis.

Run from repository root:
    python code/reproduce_statistics.py

Inputs are de-identified participant/trial metrics in data/. No raw personal identifiers
or synthetic participants are used in the observed-data estimates.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import multivariate_t
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

SEED = 20260807
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "reproduced_results"
OUT.mkdir(exist_ok=True)
METRICS = ["active_time_s", "total_path_m", "p95_angspeed_rad_s"]

phase = pd.read_csv(DATA / "subject_phase_metrics_v2.csv")
trial = pd.read_csv(DATA / "trial_metrics_anonymized.csv")
wide = phase.pivot(index=["participant", "group"], columns="phase", values=METRICS)
wide.columns = [f"{m}_{ph}" for m, ph in wide.columns]
wide = wide.reset_index()
wide["SECMA"] = (wide.group == "SECMA").astype(int)

# Baseline balance
rows = []
for m in METRICS:
    b = wide.loc[wide.group == "BOX", f"{m}_PRE"]
    s = wide.loc[wide.group == "SECMA", f"{m}_PRE"]
    sp = np.sqrt(((len(b)-1)*b.var(ddof=1)+(len(s)-1)*s.var(ddof=1))/(len(b)+len(s)-2))
    tt = stats.ttest_ind(s, b, equal_var=False)
    rows.append([m,b.mean(),b.std(ddof=1),s.mean(),s.std(ddof=1),(s.mean()-b.mean())/sp,tt.statistic,tt.pvalue])
pd.DataFrame(rows, columns=["metric","BOX_mean","BOX_sd","SECMA_mean","SECMA_sd","SMD_SECMA_minus_BOX","Welch_t","Welch_p"]).to_csv(OUT/"baseline_balance.csv", index=False)

# Paired PRE-POST
rows = []
for g in ["BOX", "SECMA"]:
    w = wide[wide.group == g]
    for m in METRICS:
        pre = w[f"{m}_PRE"].to_numpy(); post = w[f"{m}_POST"].to_numpy(); d = post-pre; n=len(d)
        md=d.mean(); se=d.std(ddof=1)/np.sqrt(n); ci=stats.t.interval(.95,n-1,loc=md,scale=se)
        tt=stats.ttest_rel(post,pre); J=1-3/(4*(n-1)-1); gz=J*md/d.std(ddof=1)
        try: wp=stats.wilcoxon(post,pre).pvalue
        except Exception: wp=np.nan
        rows.append([g,m,n,pre.mean(),post.mean(),md,ci[0],ci[1],tt.statistic,tt.pvalue,wp,gz])
pd.DataFrame(rows, columns=["group","metric","n","PRE_mean","POST_mean","change_POST_minus_PRE","CI95_low","CI95_high","t","p_paired_t","p_wilcoxon","Hedges_gz"]).to_csv(OUT/"paired_results.csv", index=False)

# Absolute-agreement ICC(A,1)
def icc_a1(d, m, ph):
    x=d[(d.phase==ph)&(d.rep<=2)].pivot(index="participant",columns="rep",values=m).dropna()
    Y=x[[1,2]].to_numpy(); n,k=Y.shape; grand=Y.mean(); rm=Y.mean(1); cm=Y.mean(0)
    MSR=k*np.sum((rm-grand)**2)/(n-1); MSC=n*np.sum((cm-grand)**2)/(k-1)
    MSE=np.sum((Y-rm[:,None]-cm[None,:]+grand)**2)/((n-1)*(k-1))
    val=(MSR-MSE)/(MSR+(k-1)*MSE+k*(MSC-MSE)/n)
    return n,val,stats.pearsonr(Y[:,0],Y[:,1]).statistic,np.mean(np.abs(Y[:,0]-Y[:,1]))
rows=[]
for ph in ["PRE","POST"]:
    for m in METRICS:
        rows.append([ph,m,*icc_a1(trial,m,ph)])
pd.DataFrame(rows,columns=["phase","metric","n_pairs","ICC_A1","Pearson_r","mean_abs_difference"]).to_csv(OUT/"reliability_icc.csv",index=False)

# HC3 ANCOVA and TOST
anc=[]; fits={}
for m in METRICS:
    X=sm.add_constant(wide[[f"{m}_PRE","SECMA"]]); y=wide[f"{m}_POST"]
    mod=sm.OLS(y,X).fit(); rb=mod.get_robustcov_results(cov_type="HC3",use_t=True)
    idx=list(X.columns).index("SECMA"); est,se,df=rb.params[idx],rb.bse[idx],rb.df_resid
    ci95=(est-stats.t.ppf(.975,df)*se,est+stats.t.ppf(.975,df)*se)
    ci90=(est-stats.t.ppf(.95,df)*se,est+stats.t.ppf(.95,df)*se)
    anc.append([m,est,se,df,rb.pvalues[idx],*ci95,*ci90,mod.rsquared]); fits[m]=mod
anc=pd.DataFrame(anc,columns=["metric","SECMA_minus_BOX","HC3_SE","df","p","CI95_low","CI95_high","CI90_low","CI90_high","R2"])
anc["p_Holm"]=multipletests(anc.p,method="holm")[1]
anc.to_csv(OUT/"ancova_primary.csv",index=False)

rows=[]; sens=[]
for m in METRICS:
    r=anc.loc[anc.metric==m].iloc[0]; est,se,df=r.SECMA_minus_BOX,r.HC3_SE,r.df
    sdpre=wide[f"{m}_PRE"].std(ddof=1)
    for k in [.2,.3,.4,.5]:
        delta=k*sdpre; p_lower=1-stats.t.cdf((est+delta)/se,df); p_upper=stats.t.cdf((est-delta)/se,df)
        rec=[m,k,delta,est,r.CI90_low,r.CI90_high,p_lower,p_upper,(p_lower<.05 and p_upper<.05)]
        sens.append(rec)
        if k==.2: rows.append(rec)
pd.DataFrame(rows,columns=["metric","margin_SD_PRE","delta_raw","estimate","CI90_low","CI90_high","p_lower","p_upper","TOST_equivalent"]).to_csv(OUT/"tost_primary.csv",index=False)
pd.DataFrame(sens,columns=["metric","margin_SD_PRE","delta_raw","estimate","CI90_low","CI90_high","p_lower","p_upper","TOST_equivalent"]).to_csv(OUT/"tost_margin_sensitivity.csv",index=False)

# Reference-posterior ROPE
rows=[]
for m in METRICS:
    mod=fits[m]; est=mod.params["SECMA"]; se=mod.bse["SECMA"]; df=mod.df_resid
    delta=.2*wide[f"{m}_PRE"].std(ddof=1)
    prob=stats.t.cdf((delta-est)/se,df)-stats.t.cdf((-delta-est)/se,df)
    rows.append([m,est,se,df,delta,prob])
pd.DataFrame(rows,columns=["metric","estimate","classical_SE","df","delta_0.2SD_PRE","posterior_P_within_ROPE"]).to_csv(OUT/"rope_reference_posterior.csv",index=False)

# Posterior-predictive operating characteristics
order=[]
for m in METRICS: order += [f"{m}_PRE",f"{m}_POST"]
Xraw=np.log(wide[order].to_numpy()); mu_global=Xraw.mean(0); sd_global=Xraw.std(0,ddof=1); Z=(Xraw-mu_global)/sd_global
d=Z.shape[1]; kappa0=.2; nu0=d+3; mu0=np.zeros(d); Psi0=(nu0-d-1)*np.eye(d)
def pred_par(Zg):
    n=len(Zg); xb=Zg.mean(0); S=(Zg-xb).T@(Zg-xb); kn=kappa0+n; nun=nu0+n; mun=(kappa0*mu0+n*xb)/kn
    q=(xb-mu0).reshape(-1,1); P=Psi0+S+(kappa0*n/kn)*(q@q.T); df=nun-d+1; shape=P*(kn+1)/(kn*df)
    return mun,shape,df
pred={g:pred_par(Z[(wide.group==g).to_numpy()]) for g in ["BOX","SECMA"]}
def sample_group(g,n,rng):
    mu,shape,df=pred[g]; z=multivariate_t.rvs(loc=mu,shape=shape,df=df,size=n,random_state=rng)
    if z.ndim==1: z=z[None,:]
    return pd.DataFrame(np.exp(z*sd_global+mu_global),columns=order)
def ancova_tost(pre,post,grp,delta):
    n=len(post); X=np.c_[np.ones(n),pre,grp]; inv=np.linalg.inv(X.T@X); beta=inv@X.T@post; resid=post-X@beta; df=n-3
    se=np.sqrt((resid@resid/df)*inv[2,2]); est=beta[2]
    pl=1-stats.t.cdf((est+delta)/se,df); pu=stats.t.cdf((est-delta)/se,df)
    return pl<.05 and pu<.05
obs_sd={m:wide[f"{m}_PRE"].std(ddof=1) for m in METRICS}
rng=np.random.default_rng(SEED+1); rows=[]
for N in [18,30,40,60,80,120,200,400]:
    counts={(m,k):0 for m in METRICS for k in [.2,.3,.4,.5]}
    R=300
    for _ in range(R):
        nb=N//2; ns=N-nb; b=sample_group("BOX",nb,rng); ss=sample_group("SECMA",ns,rng); b["grp"]=0; ss["grp"]=1
        sim=pd.concat([b,ss],ignore_index=True)
        for m in METRICS:
            for k in [.2,.3,.4,.5]:
                counts[(m,k)] += ancova_tost(sim[f"{m}_PRE"].to_numpy(),sim[f"{m}_POST"].to_numpy(),sim.grp.to_numpy(),k*obs_sd[m])
    for m in METRICS:
        for k in [.2,.3,.4,.5]: rows.append([N,m,k,k*obs_sd[m],counts[(m,k)]/R,R])
pd.DataFrame(rows,columns=["N_total","metric","margin_SD_PRE","delta_raw","P_TOST_equivalence","MC_replicates"]).to_csv(OUT/"posterior_predictive_power.csv",index=False)
print(f"Results written to {OUT}")
