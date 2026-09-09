#!/usr/bin/env python3
"""Measured full-cohort preservation summary; all three seeds remain visible."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
data={name:json.loads((a.root/('reference-'+name)/'checks.json').read_text()) for name in ['kang','hagai']}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(13,4.7));colors=['#36495c','#177a87','#9874bd'];labels=['Original PCA','Native','Harmony 2']
for ax,key,title in [(axes[0],'sameDonorExcess','Donor mixing within type / condition'),(axes[2],'conditionBalancedAccuracy','Treatment accuracy across donors')]:
 for cohort,offset in [('kang',0),('hagai',4)]:
  d=data[cohort]
  values=[[d['baseline'][key]],[r['metrics'][key] for r in d['native']],[r['metrics'][key] for r in d['references']]]
  for i,v in enumerate(values):
   ax.scatter(np.repeat(offset+i,len(v))+np.linspace(-.09,.09,len(v)),v,color=colors[i],s=35,zorder=3)
   ax.plot([offset+i-.2,offset+i+.2],[np.mean(v)]*2,color=colors[i],linewidth=2)
 ax.set_xticks([0,1,2,4,5,6],labels*2,rotation=32,ha='right');ax.set_title(title,fontsize=11);ax.grid(axis='y',alpha=.18)
 ax.text(1,-.29,'Kang · 24,673 cells',ha='center',transform=ax.get_xaxis_transform());ax.text(5,-.29,'Hagai · 13,863 cells',ha='center',transform=ax.get_xaxis_transform())
axes[0].set_ylabel('Excess same-donor neighbor fraction');axes[0].set_ylim(-.02,.54)
axes[2].set_ylabel('Mean balanced accuracy');axes[2].set_ylim(.9,1.006)
d=data['kang'];baseline=d['baseline']['cellTypeRecall']['NK cells'];ax=axes[1]
for i,v in enumerate([[baseline],[r['metrics']['cellTypeRecall']['NK cells'] for r in d['native']],[r['metrics']['cellTypeRecall']['NK cells'] for r in d['references']]]):
 ax.scatter(np.repeat(i,len(v))+np.linspace(-.09,.09,len(v)),v,s=42,color=colors[i],zorder=3)
ax.axhline(baseline-.05,color='#c63c46',linestyle='--',linewidth=1.4,label='Maximum accepted loss: 5 points');ax.set_ylim(.80,.98);ax.set_xticks(range(3),labels,rotation=32,ha='right');ax.set_ylabel('Donor-balanced NK-cell recall');ax.set_title('Kang NK-cell preservation fails',fontsize=11,color='#ad2631');ax.legend(loc='upper right',fontsize=8);ax.grid(axis='y',alpha=.18)
fig.suptitle('Integration storage is numerically exact; biological preservation depends on the cohort',fontsize=13,y=1.02)
fig.text(.02,-.015,'Points are individual seeds 7 / 19 / 41. Same native PCA inputs for both methods; labels used for evaluation.\nKang has four unavailable rare-type treatment folds. These are engineering checks, not prospective or experimental validation.',fontsize=9,color='#43515c')
fig.tight_layout();fig.subplots_adjust(bottom=.30,wspace=.38);fig.savefig(a.out/'preservation.png',dpi=180,bbox_inches='tight');fig.savefig(a.out/'preservation.pdf',bbox_inches='tight')
