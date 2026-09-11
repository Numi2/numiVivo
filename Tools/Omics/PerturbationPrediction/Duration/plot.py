#!/usr/bin/env python3
"""Show all held-out donor-times and all distinct point predictors without smoothing."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import sha,write
p=argparse.ArgumentParser();p.add_argument('--scores',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=True)
rows=json.loads((a.scores/'scores.json').read_text());fig,axes=plt.subplots(1,3,figsize=(11.5,4.4),sharey=True)
for ax,donor in zip(axes,sorted({r['donor'] for r in rows})):
 times=sorted({r['hours'] for r in rows if r['donor']==donor});values={(r['hours'],r['method']):r['RMSE'] for r in rows if r['donor']==donor}
 for method,label,color,marker in [('noChange','No change','#858585','.'),('matchedInvariantMean','Matched time-invariant mean','#ce8834','^'),('durationMean','Duration mean / median','#087e8b','o'),('durationContextRidge','Duration context ridge','#a65477','s')]:
  ax.plot(times,[values[h,method] for h in times],color=color,marker=marker,markersize=4,linewidth=1.5,label=label)
 ax.set_title(donor.rsplit(':',1)[1]);ax.set_xscale('log',base=2);ax.set_xlim(.85,42);ax.set_xticks([1,2,4,8,12,24,36],labels=['1','2','4','8','12','24','36']);ax.set_xlabel('IFN-beta exposure (hours; log scale)');ax.grid(axis='y',alpha=.17);ax.spines[['top','right']].set_visible(False)
axes[0].set_ylabel('All-panel treated logRNA RMSE (lower is better)');fig.suptitle('Exposure duration improves held-out donor RNA prediction in this study',fontsize=13,y=.97);fig.legend(*axes[-1].get_legend_handles_labels(),frameon=False,loc='upper center',bbox_to_anchor=(.5,.89),ncol=4,fontsize=9)
fig.text(.5,.025,'Retrospective development: 3 donor holdouts, all 18 outcomes, 12,993 genes. Mean = median with 2 training donors.\nAll points retained; lines connect observed query times. No independent external validation or phenotype claim.',ha='center',fontsize=9,color='#444444');fig.tight_layout(rect=[0,.14,1,.88]);fig.savefig(a.out/'duration-development.svg',metadata={'Date':None});fig.savefig(a.out/'duration-development.png',dpi=170);plt.close(fig)
svg=a.out/'duration-development.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n');write(a.out/'plot.json',dict(sourceScoresSHA256=sha(a.scores/'scores.json'),scriptSHA256=sha(Path(__file__)),donors=3,donorTimes=18,allMethodsShown=True,meanAndMedianExact=True,smoothing=False))
