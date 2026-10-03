from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/neem-matplotlib')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

B=Path(__file__).resolve().parent
R=B/'results'; F=B/'figure9_expression'; F.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
e=pd.read_csv(R/'expression_results.csv')
colours=['#4C78A8','#D87532']
def save(fig,name):
    for ext in ['png','pdf']: fig.savefig(F/f'{name}.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)

fig,axes=plt.subplots(2,3,figsize=(11,8.2))
for row,ds in enumerate(['gse30784','gse9844']):
    x=pd.read_csv(R/f'{ds}_target_expression.csv',index_col=0)
    m=pd.read_csv(R/f'{ds}_samples.csv',index_col=0)
    rng=np.random.default_rng(20261003)
    for ax,g in zip(axes[row],['EGFR','MMP9','MMP7']):
        vals=[x.loc[g,m.index[m.group==k]].to_numpy(float) for k in ['Normal','Tumour']]
        bp=ax.boxplot(vals,positions=[0,1],widths=.48,patch_artist=True,showfliers=False,medianprops={'color':'black'})
        for i,(v,c) in enumerate(zip(vals,colours)):
            bp['boxes'][i].set_facecolor(c);bp['boxes'][i].set_alpha(.3)
            ax.scatter(rng.normal(i,.055,len(v)),v,s=12,alpha=.48,color=c,edgecolors='none')
        r=e[(e.dataset==ds)&(e.gene==g)].iloc[0]
        ax.set_title(f'{g}\nDifference = {r.log2FC:.2f}; q = {r.q:.3g}',fontsize=10)
        ax.set_xticks([0,1],[f'Normal\n(n={len(vals[0])})',f'Tumour\n(n={len(vals[1])})'])
        ax.set_ylabel('Normalised log₂ expression');ax.grid(axis='y',alpha=.18)
    axes[row,0].text(-.28,1.24,f'{chr(65+row)}  {ds.upper()}',transform=axes[row,0].transAxes,fontsize=12,weight='bold')
fig.suptitle('Figure 9. EGFR, MMP9 and MMP7 expression in oral cancer',fontsize=14,y=.98)
fig.subplots_adjust(top=.86,bottom=.09,left=.09,right=.98,hspace=.75,wspace=.4)
save(fig,'Figure_9_expression_both_cohorts')
print(F.resolve())
