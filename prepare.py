from pathlib import Path
import pandas as pd,numpy as np,gzip,csv,io,json,re
B=Path(__file__).resolve().parent;O=B/'results';O.mkdir(exist_ok=True)
GENES=['EGFR','CASP3','ERBB2','MMP9','PIK3CA','PIK3CG','PIK3CD','PIK3CB','MDM2','VDR','MMP2','MMP7','BRAF','EDNRB']
def tab_gz(p,tag):
 with gzip.open(p,'rt') as f:
  meta=[]
  for l in f:
   if l.startswith(tag):break
   meta.append(next(csv.reader([l],delimiter='\t')))
  df=pd.read_csv(f,sep='\t',comment='!',dtype={0:str})
 return meta,df
_,ann=tab_gz(B/'data/gpl570.gz','!platform_table_begin');ann=ann.dropna(subset=['Gene symbol']);ann['Gene symbol']=ann['Gene symbol'].str.strip();ann=ann[ann['Gene symbol'].isin(GENES)]
assert not ann['ID'].duplicated().any();ann[['ID','Gene symbol','Gene ID']].to_csv(O/'probe_mapping.csv',index=False)
qc={}
for key in ['gse30784','gse9844']:
 meta,x=tab_gz(B/'data'/f'{key}.gz','!series_matrix_table_begin');x=x.set_index('ID_REF');x.index=x.index.astype(str)
 q=np.nanquantile(x.to_numpy(float),[0,.01,.5,.99,1]);print(key,x.shape,q)
 assert q[-1]<30 and q[0]>=0,'Inspect scale before transform'
 md={}
 for a in meta:
  if not a:continue
  if a[0] in ['!Sample_geo_accession','!Sample_title','!Sample_source_name_ch1']:md[a[0].replace('!Sample_','')]=a[1:]
  if a[0]=='!Sample_characteristics_ch1':
   if key=='gse30784':name=a[1].split(':')[0].lower();md[name]=[v.split(':',1)[1].strip() for v in a[1:]]
   else:md['characteristics']=a[1:]
 sm=pd.DataFrame(md).set_index('geo_accession');assert sm.index.tolist()==x.columns.tolist()
 if key=='gse30784':sm['group']=sm['status'].map({'cancer':'Tumour','control':'Normal','dysplasia':'Dysplasia'})
 else:sm['group']=np.where(sm['title'].str.contains('normal control'),'Normal','Tumour')
 assert sm.group.notna().all();sm.to_csv(O/f'{key}_samples.csv')
 mapped=x.join(ann.set_index('ID')['Gene symbol'],how='inner');gx=mapped.groupby('Gene symbol').median().reindex(GENES)
 assert gx.notna().all().all();gx.to_csv(O/f'{key}_target_expression.csv')
 # Basic processed-data QC without outcome-dependent sample removal
 qc[key]={'matrix_shape':list(x.shape),'quantiles':q.tolist(),'groups':sm.group.value_counts().to_dict(),'probe_counts':mapped.groupby('Gene symbol').size().to_dict(),'sample_median_range':[float(x.median().min()),float(x.median().max())]}
 quant=pd.DataFrame({'sample':x.columns,'q25':x.quantile(.25),'median':x.median(),'q75':x.quantile(.75)});quant.to_csv(O/f'{key}_sample_qc.csv',index=False)
 del x
p=pd.read_csv(B/'data/clinical_patient.tsv',sep='\t',comment='#');s=pd.read_csv(B/'data/clinical_sample.tsv',sep='\t',comment='#');xc=pd.read_csv(B/'data/xena_clinical.tsv',sep='\t')
xc['PATIENT_ID']=xc['sampleID'].str[:12]
assert xc.groupby('PATIENT_ID')['anatomic_neoplasm_subdivision'].nunique().max()==1
sites=xc.drop_duplicates('PATIENT_ID').set_index('PATIENT_ID')['anatomic_neoplasm_subdivision'];p['site']=p.PATIENT_ID.map(sites)
print('SAMPLE types',s.SAMPLE_TYPE.value_counts().to_dict());print('site',p.site.value_counts(dropna=False).to_dict())
expr=pd.read_csv(B/'data/tcga_expression.tsv',sep='\t');print('TCGA',expr.shape,expr.iloc[:2,:5].to_string(index=False))
expr=expr[expr.Hugo_Symbol.isin(GENES)];assert expr.Hugo_Symbol.nunique()==14;assert not expr.Hugo_Symbol.duplicated().any()
expr=expr.set_index('Hugo_Symbol').drop(columns='Entrez_Gene_Id').T;expr.index.name='SAMPLE_ID';expr=expr.apply(pd.to_numeric,errors='coerce')
cs=s.merge(p,on='PATIENT_ID',validate='many_to_one').merge(expr,left_on='SAMPLE_ID',right_index=True,how='left',validate='one_to_one')
include={'Oral Tongue','Oral Cavity','Floor of mouth','Buccal Mucosa','Alveolar Ridge','Hard Palate'}
cs['oral_site']=cs.site.isin(include);cs['time']=pd.to_numeric(cs.OS_MONTHS,errors='coerce');cs['event']=cs.OS_STATUS.map({'0:LIVING':0,'1:DECEASED':1});cs['age']=pd.to_numeric(cs.AGE,errors='coerce');cs['male']=cs.SEX.map({'Male':1,'Female':0});cs['advanced']=cs.AJCC_PATHOLOGIC_TUMOR_STAGE.map({'STAGE I':0,'STAGE II':0,'STAGE III':1,'STAGE IVA':1,'STAGE IVB':1,'STAGE IVC':1})
cs['exclusion']='included'
cs.loc[~cs.oral_site,'exclusion']='Non-oral site or missing site'
cs.loc[cs.oral_site & (cs.SAMPLE_TYPE!='Primary'),'exclusion']='Not primary tumour'
cs.loc[(cs.exclusion=='included')&cs[GENES].isna().any(axis=1),'exclusion']='Missing target expression'
cs.loc[(cs.exclusion=='included')&(~(cs.time>0)|cs.event.isna()),'exclusion']='Missing or nonpositive survival endpoint'
cs[['PATIENT_ID','SAMPLE_ID','site','SAMPLE_TYPE','oral_site','exclusion','ICD_O_3_SITE']].to_csv(O/'tcga_inclusion_audit.csv',index=False)
co=cs[cs.exclusion=='included'].copy();assert not co.PATIENT_ID.duplicated().any();assert (co[GENES]>=0).all().all()
for g in GENES:co[g]=np.log2(co[g]+1)
co[['PATIENT_ID','SAMPLE_ID','site','ICD_O_3_SITE','time','event','age','male','advanced','AJCC_PATHOLOGIC_TUMOR_STAGE']+GENES].to_csv(O/'tcga_analysis_data.csv',index=False)
qc['tcga']={'source_patients':len(p),'clinical_samples':len(s),'expression_samples':len(expr),'site_counts':p.site.value_counts(dropna=False).to_dict(),'flow':cs.exclusion.value_counts().to_dict(),'unadjusted_n':len(co),'events':int(co.event.sum()),'adjusted_n':len(co.dropna(subset=['age','male','advanced'])),'adjusted_events':int(co.dropna(subset=['age','male','advanced']).event.sum()),'missing_covariates':co[['age','male','advanced']].isna().sum().to_dict(),'retained_sites':co.site.value_counts().to_dict()}
json.dump(qc,open(O/'cohort_qc.json','w'),indent=2);print(json.dumps(qc,indent=2))
