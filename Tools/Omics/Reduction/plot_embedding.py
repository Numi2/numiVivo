#!/usr/bin/env python3
"""Render measured-data native coordinates; labels are not biological validation."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser()
p.add_argument('--kang',type=Path,required=True);p.add_argument('--pbmc',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
fig,axes=plt.subplots(2,2,figsize=(10,8),layout='constrained')
for row,(name,path) in enumerate([('Kang B cells',a.kang),('PBMC3k',a.pbmc)]):
    report=json.loads(path.read_text());embedding=report['embedding']
    if row==0:
        samples={s['id']:s for s in report['processed']['dataset']['samples']}
        groups=[samples[c['sampleID']]['condition'] for c in embedding['cells']]
        legend_title='Declared condition'
    else:
        groups=[str(v) for v in report['clustering']['labels']]
        legend_title='Graph community'
    categories=sorted(set(groups),key=int) if row else sorted(set(groups)); labels=np.asarray(groups)
    for col,key in enumerate(['initialCoordinates','coordinates']):
        xy=np.asarray(embedding[key]);axis=axes[row,col]
        for index,group in enumerate(categories):
            selected=labels==group
            axis.scatter(xy[selected,0],xy[selected,1],s=5,alpha=.65,color=plt.get_cmap('tab20' if row else 'tab10')(index),label=group,linewidths=0)
        axis.set_title(name+' — '+('initial PCA' if col==0 else 'native UMAP'))
        axis.set_xlabel('Coordinate 1');axis.set_ylabel('Coordinate 2')
        if col==1: axis.legend(title=legend_title,fontsize=7,title_fontsize=8,loc='upper left',bbox_to_anchor=(1,1),markerscale=2)
fig.suptitle('Fixed-epoch native embeddings — descriptive coordinates\nDistances and apparent densities are not calibrated biological measurements',fontsize=11)
a.out.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.out,dpi=160);plt.close(fig)
