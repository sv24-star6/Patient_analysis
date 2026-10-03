from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/neem-matplotlib')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

B=Path(__file__).resolve().parent
R=B/'results'; F=B/'figure10_survival'; F.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
e=pd.read_csv(R/'expression_results.csv')
colours=['#4C78A8','#D87532']
def save(fig,name):
    for ext in ['png','pdf']: fig.savefig(F/f'{name}.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)

genes=e[e.dataset=='gse30784'].gene.tolist(); assert len(genes)==14
r=pd.read_csv(R/'survival_results.csv').query("model == 'adjusted'").set_index('gene').loc[genes]
assert (r.n==286).all() and (r.events==129).all()
fig,(ax,txt)=plt.subplots(1,2,figsize=(10,7),gridspec_kw={'width_ratios':[1.8,1]},sharey=True)
y=np.arange(14)
ax.errorbar(r.HR,y,xerr=[r.HR-r.ci_low,r.ci_high-r.HR],fmt='o',color='#386B79',capsize=3)
ax.axvline(1,color='grey',linestyle='--',linewidth=1);ax.set_xscale('log')
ax.set_xlim(.65,1.7);ax.set_xticks([.7,1,1.5],['0.7','1.0','1.5']);ax.set_yticks(y,genes);ax.set_ylim(13.7,-1.5)
ax.set_xlabel('Adjusted hazard ratio per 1 SD\n(95% confidence interval)');ax.grid(axis='x',alpha=.15)
txt.axis('off');txt.text(0,-.8,'HR (95% CI)',weight='bold');txt.text(.84,-.8,'q',weight='bold')
for i,(_,v) in enumerate(r.iterrows()):
    txt.text(0,i,f'{v.HR:.2f} ({v.ci_low:.2f}–{v.ci_high:.2f})',va='center');txt.text(.84,i,f'{v.q:.3f}',va='center')
fig.suptitle('Figure 10. Overall survival in oral-site TCGA cases',fontsize=13,y=.98)
fig.text(.12,.025,'Age-, sex- and stage-adjusted models; n = 286, deaths = 129.\nNo target met BH-adjusted q < 0.05. Confidence intervals are pointwise, not multiplicity-adjusted.',fontsize=9)
fig.tight_layout(rect=[0,.085,1,.94]);save(fig,'Figure_10_adjusted_survival_forest')
print(f'Created Figure 10 in PNG and PDF: {F}')
