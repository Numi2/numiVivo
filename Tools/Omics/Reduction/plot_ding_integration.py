#!/usr/bin/env python3
"""Show all declared Ding seeds and unchanged preservation thresholds."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
b=json.loads((a.root/'evaluation-inputs-final/baseline.json').read_text());native=json.loads((a.root/'native-evaluation/checks.json').read_text());reference=json.loads((a.root/'reference-evaluation/checks.json').read_text())
assert native['status']==reference['status']=='all-declared-runs-measured'
groups=[('Original PCA',[b])]
for name,results in [('Native',native),('Harmony',reference)]:
 for mode in ['fixed','adaptive']:groups.append((name+' '+mode,[r['metrics'] for r in results['runs'] if r['mode']==mode]))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(13,4.8));colors=['#475569','#b45309','#087f8c','#9a6d39','#7654a3']
panels=[('Excess same-method neighbors',lambda m:m['sameMethodExcess'],None),('Plasmacytoid dendritic-cell recall',lambda m:m['cellTypeRecall']['Plasmacytoid dendritic cell'],.05),('Within-stratum T-receptor program',lambda m:m['programs']['T-receptor']['withinStratumSpearman'],.05)]
for ax,(title,get,margin) in zip(axes,panels):
 for i,(label,values) in enumerate(groups):
  values=[get(v) for v in values];ax.scatter(i+np.linspace(-.07,.07,len(values)),values,s=35,color=colors[i],zorder=3);ax.plot([i-.18,i+.18],[np.mean(values)]*2,color=colors[i],linewidth=2)
 if margin is not None:ax.axhline(get(b)-margin,color='#b91c1c',linestyle='--',linewidth=1.2,label='Declared preservation threshold');ax.legend(fontsize=8)
 ax.set_title(title,fontsize=11);ax.set_xticks(range(len(groups)),[g[0] for g in groups],rotation=35,ha='right');ax.grid(axis='y',alpha=.2)
axes[2].set_ylabel('Spearman correlation')
fig.suptitle('Independent Ding study: mixing improves while preservation gates fail',fontsize=13)
fig.text(.025,.015,'44,031 input cells; 29,411 source-assigned labels; seeds 7 / 19 / 41. Lower method excess means better mixing.\nHoldouts are source experiments, with unverified donor identity. Count programs are diagnostic proxies; annotation coverage is partial.',fontsize=9,color='#475569')
fig.tight_layout();fig.subplots_adjust(bottom=.31,top=.87,wspace=.3)
for extension in ['png','pdf']:fig.savefig(a.out/('ding-preservation.'+extension),dpi=180,bbox_inches='tight')
