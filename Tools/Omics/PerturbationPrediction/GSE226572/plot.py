#!/usr/bin/env python3
"""Plot every released donor/time outcome; no smoothing or selected time window."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from download import sha,write
p=argparse.ArgumentParser();p.add_argument('--scores',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=True)
rows=json.loads((a.scores/'scores.json').read_text());donors=sorted({r['donor'] for r in rows});fig,axes=plt.subplots(1,3,figsize=(11.5,4.1),sharey=True)
for ax,donor in zip(axes,donors):
 times=sorted({r['hours'] for r in rows if r['donor']==donor});values={(r['hours'],r['method']):r['RMSE'] for r in rows if r['donor']==donor}
 for method,label,color,marker in [('meanResponse','Mean response','#087e8b','o'),('contextRidge','Context ridge','#a65477','s')]:
  ax.plot(times,[100*(1-values[t,method]/values[t,'noChange']) for t in times],marker=marker,color=color,label=label,linewidth=1.6,markersize=4)
 ax.axhline(0,color='#555555',linewidth=.9,linestyle='--');ax.set_title(donor.rsplit(':',1)[1],fontsize=12);ax.set_xscale('log',base=2);ax.set_xlim(.85,42);ax.set_xticks([1,2,4,8,12,24,36],labels=['1','2','4','8','12','24','36']);ax.set_xlabel('IFN-beta exposure (hours; log scale)');ax.grid(axis='y',alpha=.17);ax.spines[['top','right']].set_visible(False)
axes[0].set_ylabel('RMSE reduction versus no change (%)');fig.legend(*axes[-1].get_legend_handles_labels(),frameon=False,fontsize=10,loc='upper center',bbox_to_anchor=(.5,.88),ncol=2);fig.suptitle('Fixed Kang response across all released IFN-beta durations',fontsize=14,y=.97)
fig.text(.5,.025,'Positive values improve on no change. Three query donors; each predictor is unchanged across time.\nWhole-population initial-QC RNA; assay, culture, health and composition also differ.',ha='center',fontsize=9,color='#444444');fig.tight_layout(rect=[0,.12,1,.90]);fig.savefig(a.out/'duration-transfer.svg',metadata={'Date':None});fig.savefig(a.out/'duration-transfer.png',dpi=170);plt.close(fig)
svg=a.out/'duration-transfer.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
write(a.out/'plot.json',dict(sourceScoresSHA256=sha(a.scores/'scores.json'),scriptSHA256=sha(Path(__file__)),donors=3,donorTimes=18,allTimePointsIncluded=True,smoothing=False))
