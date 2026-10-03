from pathlib import Path
import pandas as pd,numpy as np,json,scipy,scipy.stats as st,statsmodels.api as sm,statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from lifelines import CoxPHFitter,KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test,logrank_test
import lifelines,statsmodels,warnings
B=Path(__file__).resolve().parent;O=B/'results'
G=['EGFR','CASP3','ERBB2','MMP9','PIK3CA','PIK3CG','PIK3CD','PIK3CB','MDM2','VDR','MMP2','MMP7','BRAF','EDNRB']
rows=[];adj=[]
for ds in ['gse30784','gse9844']:
 x=pd.read_csv(O/f'{ds}_target_expression.csv',index_col=0);meta=pd.read_csv(O/f'{ds}_samples.csv',index_col=0);meta=meta[meta.group.isin(['Tumour','Normal'])]
 for g in G:
  a=x.loc[g,meta.index[meta.group=='Tumour']].to_numpy(float);b=x.loc[g,meta.index[meta.group=='Normal']].to_numpy(float)
  t=st.ttest_ind(a,b,equal_var=False);diff=a.mean()-b.mean();se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b));v1=a.var(ddof=1)/len(a);v2=b.var(ddof=1)/len(b);df=(v1+v2)**2/(v1*v1/(len(a)-1)+v2*v2/(len(b)-1));ci=st.t.ppf(.975,df)*se
  rows.append(dict(dataset=ds,gene=g,n_tumour=len(a),n_normal=len(b),tumour_mean=a.mean(),normal_mean=b.mean(),log2FC=diff,ci_low=diff-ci,ci_high=diff+ci,p=t.pvalue))
  if ds=='gse30784':
   d=meta.copy();d['expression']=x.loc[g,d.index];d['tumour']=(d.group=='Tumour').astype(int)
   fit=smf.ols('expression ~ tumour + C(age) + C(sex)',data=d).fit(cov_type='HC3');ci2=fit.conf_int().loc['tumour'];adj.append(dict(gene=g,n=int(fit.nobs),log2FC=fit.params['tumour'],ci_low=ci2.iloc[0],ci_high=ci2.iloc[1],p=fit.pvalues['tumour']))
r=pd.DataFrame(rows);r['q']=r.groupby('dataset')['p'].transform(lambda p:multipletests(p,method='fdr_bh')[1]);r.to_csv(O/'expression_results.csv',index=False)
a=pd.DataFrame(adj);a['q']=multipletests(a.p,method='fdr_bh')[1];a.to_csv(O/'expression_adjusted_sensitivity.csv',index=False)
co=pd.read_csv(O/'tcga_analysis_data.csv');sur=[];ph=[];coeff=[];km=[];log=[]
# Every gene uses its SD in the full eligible survival cohort, also in adjusted/subsite analyses.
z=(co[G]-co[G].mean())/co[G].std(ddof=1);pd.DataFrame({'mean_log2RSEM':co[G].mean(),'sd_log2RSEM':co[G].std(ddof=1)}).to_csv(O/'tcga_expression_standardisation.csv')
for model in ['unadjusted','adjusted','specific_site_adjusted']:
 for g in G:
  d=co[['time','event','age','male','advanced','site']].copy();d['expression_z']=z[g]
  if model=='specific_site_adjusted':d=d[d.site!='Oral Cavity']
  cols=['time','event','expression_z']+(['age','male','advanced'] if model!='unadjusted' else [])
  d=d[cols].dropna();assert d.event.sum()>10*(len(cols)-2)
  with warnings.catch_warnings(record=True) as ww:
   warnings.simplefilter('always');c=CoxPHFitter();c.fit(d,duration_col='time',event_col='event');pt=proportional_hazard_test(c,d,time_transform='rank').summary
  for w in ww:log.append({'gene':g,'model':model,'warning':str(w.message)})
  s=c.summary.loc['expression_z'];sur.append(dict(model=model,gene=g,n=len(d),events=int(d.event.sum()),HR=s['exp(coef)'],ci_low=s['exp(coef) lower 95%'],ci_high=s['exp(coef) upper 95%'],p=s.p,ph_p=pt.loc['expression_z','p'],concordance=c.concordance_index_))
  for cov in pt.index:ph.append(dict(model=model,gene=g,covariate=cov,p=pt.loc[cov,'p']))
  for cov,row in c.summary.iterrows():coeff.append(dict(model=model,gene=g,covariate=cov,coefficient=row['coef'],HR=row['exp(coef)'],p=row.p))
sr=pd.DataFrame(sur);sr['q']=sr.groupby('model')['p'].transform(lambda p:multipletests(p,method='fdr_bh')[1]);sr.to_csv(O/'survival_results.csv',index=False);pd.DataFrame(ph).to_csv(O/'proportional_hazards_diagnostics.csv',index=False);pd.DataFrame(coeff).to_csv(O/'cox_all_coefficients.csv',index=False);json.dump(log,open(O/'model_warnings.json','w'),indent=2)
for g in G:
 high=co[g]>co[g].median();t=logrank_test(co.loc[high,'time'],co.loc[~high,'time'],co.loc[high,'event'],co.loc[~high,'event']);km.append(dict(gene=g,median_log2RSEM=co[g].median(),n_high=int(high.sum()),n_low=int((~high).sum()),p=t.p_value))
k=pd.DataFrame(km);k['q']=multipletests(k.p,method='fdr_bh')[1];k.to_csv(O/'median_split_logrank_secondary.csv',index=False)
# Integrated results retain all targets, regardless of significance.
e1=r[r.dataset=='gse30784'].set_index('gene');e2=r[r.dataset=='gse9844'].set_index('gene');ss=sr[sr.model=='adjusted'].set_index('gene')
integ=pd.DataFrame(index=G)
for tag,e in [('discovery',e1),('replication',e2)]:
 for col in ['log2FC','q']:integ[tag+'_'+col]=e[col]
integ['replicated_expression']=(integ.discovery_log2FC*integ.replication_log2FC>0)&(integ.discovery_q<.05)&(integ.replication_q<.05)
for c in ['HR','ci_low','ci_high','p','q','ph_p']:integ['adjusted_OS_'+c]=ss[c]
integ.index.name='gene';integ.to_csv(O/'integrated_target_evidence.csv')
json.dump({'pandas':pd.__version__,'numpy':np.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__,'lifelines':lifelines.__version__},open(O/'software_versions.json','w'),indent=2)
print(integ.to_string());print('Warnings',log);print('PH',pd.DataFrame(ph).query('model=="adjusted" and p<0.05').to_string(index=False))
