#!/usr/bin/env python3
"""Plot every declared seed, retaining fixed-ridge failures and adaptive outcomes."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
p=argparse.ArgumentParser(description=__doc__)
for name in ['frozen','adaptive','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
fixed=json.loads((a.frozen/'reference-kang/checks.json').read_text());adaptive=json.loads((a.adaptive/'evaluation-kang/kang/checks.json').read_text())
assert adaptive['status']=='full-cohort-adaptive-native-and-reference-measured'
groups=[('Original PCA',[fixed['baseline']]),('Native fixed',[r['metrics'] for r in fixed['native']]),('Native adaptive',[r['metrics'] for r in adaptive['native']]),('Harmony fixed',[r['metrics'] for r in fixed['references']]),('Harmony adaptive',[r['metrics'] for r in adaptive['references']])]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(13,4.8));colors=['#475569','#b45309','#087f8c','#9a6d39','#7654a3']
panels=[('NK-cell recall',lambda m:m['cellTypeRecall']['NK cells'],.05),('Treatment balanced accuracy',lambda m:m['conditionBalancedAccuracy'],.02),('Excess same-donor neighbors',lambda m:m['sameDonorExcess'],None)]
for ax,(title,get,margin) in zip(axes,panels):
 for i,(label,values) in enumerate(groups):
  v=[get(m) for m in values];ax.scatter(i+np.linspace(-.07,.07,len(v)),v,s=35,color=colors[i],zorder=3)
  ax.plot([i-.18,i+.18],[np.mean(v)]*2,color=colors[i],linewidth=2)
 if margin is not None:ax.axhline(get(fixed['baseline'])-margin,color='#b91c1c',linestyle='--',linewidth=1.2,label='Declared preservation threshold');ax.legend(fontsize=8,loc='best')
 ax.set_title(title,fontsize=11);ax.set_xticks(range(len(groups)),[g[0] for g in groups],rotation=35,ha='right');ax.grid(axis='y',alpha=.2)
fig.suptitle('Full Kang cohort: fixed and expected-mass ridge, all three seeds',fontsize=13)
fig.text(.025,.015,'24,673 cells; seeds 7 / 19 / 41; identical PCA inputs and unchanged margins. Lower donor excess means better mixing.\nRetrospective evaluation; unavailable rare-type folds remain unavailable. Source labels are evaluation annotations.',fontsize=9,color='#475569')
fig.tight_layout();fig.subplots_adjust(bottom=.31,top=.87,wspace=.3)
for extension in ['png','pdf']:fig.savefig(a.out/('adaptive-preservation.'+extension),dpi=180,bbox_inches='tight')
