#!/usr/bin/env python3
"""Plot retained integration alternatives without hiding failed configurations."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
fig,axes=plt.subplots(1,3,figsize=(13,4.9));colors=['#64748b','#d97706','#0f766e','#be123c','#2563eb']
for ax,cohort,title,get in [(axes[0],'kang','Kang NK recall',lambda m:m['cellTypeRecall']['NK cells']),(axes[1],'ding','Ding pDC recall',lambda m:m['cellTypeRecall']['Plasmacytoid dendritic cell']),(axes[2],'ding','Ding within-stratum T-receptor',lambda m:m['programs']['T-receptor']['withinStratumSpearman'])]:
 d=json.loads((a.root/cohort/'checks.json').read_text());raw=json.loads((a.root/('scanorama-'+cohort)/'checks.json').read_text());scaled=json.loads((a.root/('scanorama-scaled-'+cohort)/'checks.json').read_text())
 groups=[('Original PCA',[d['baseline']]),('Rigid fixed',[x['metrics'] for x in d['runs'] if x['mode']=='fixed']),('Rigid adaptive',[x['metrics'] for x in d['runs'] if x['mode']=='adaptive']),('Scanorama raw',[raw['metrics']]),('Scanorama scaled',[scaled['metrics']])]
 for i,(name,ms) in enumerate(groups):
  values=[get(m) for m in ms];ax.scatter(i+np.linspace(-.07,.07,len(values)),values,color=colors[i],s=32,zorder=3);ax.plot([i-.17,i+.17],[np.mean(values)]*2,color=colors[i],linewidth=2)
 ax.axhline(get(d['baseline'])-.05,color='#991b1b',linestyle='--',linewidth=1,label='Unchanged preservation margin');ax.set_title(title,fontsize=11);ax.set_xticks(range(len(groups)),[x[0] for x in groups],rotation=35,ha='right');ax.grid(axis='y',alpha=.2);ax.spines[['top','right']].set_visible(False)
axes[0].legend(fontsize=8,loc='lower left');fig.suptitle('Integration development: geometry constraint is insufficient; reference scale matters',fontsize=12)
fig.text(.025,.015,'All original cells retained; rigid seeds 7 / 19 / 41. Scaled Scanorama uses one global median PCA norm, with identical anchors.\nPreviously inspected cohorts; not independent validation. These panels are examples; every type and program gate remains in the report.',fontsize=9,color='#475569')
fig.tight_layout();fig.subplots_adjust(bottom=.30,top=.87,wspace=.28)
for extension in ['png','pdf']:fig.savefig(a.out/('integration-alternatives.'+extension),dpi=180,bbox_inches='tight')
