from pathlib import Path
import argparse,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser(description='Render complete native and reference cohort layouts without inferred labels');p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();b=a.inputs;out=a.out;out.mkdir(parents=True,exist_ok=False)
fig,axes=plt.subplots(2,3,figsize=(13,8));sources={}
for col,name in enumerate(['baron','hagai','norman']):
 p=b/'real-evidence'/name/'embedding/result.json';r=json.loads(p.read_text());q=b/('reference-'+name)/'reference-7.npz';reference=np.load(q)['coordinates'];native=np.asarray(r['coordinates'])
 for row,(kind,values,color) in enumerate([('Native',native,'#176b75'),('umap-learn, seed 7',reference,'#76529e')]):
  ax=axes[row,col];ax.scatter(values[:,0],values[:,1],s=1.1 if name=='norman' else 2,alpha=.35,linewidths=0,c=color,rasterized=True)
  ax.set_title(f'{name.title()} · {len(values):,} cells\n{kind}',fontsize=11);ax.set_xlabel('Coordinate 1');ax.set_ylabel('Coordinate 2');ax.spines[['top','right']].set_visible(False)
 sources[name]=dict(nativeResultSHA256=hashlib.sha256(p.read_bytes()).hexdigest(),referenceCoordinatesSHA256=hashlib.sha256(q.read_bytes()).hexdigest())
fig.suptitle('Complete-cohort 500-epoch UMAP-compatible layouts',fontsize=15)
fig.text(.5,.015,'All cells shown. Axes have arbitrary orientation and scale; no cell types inferred. Norman quality metrics use 2,048 fixed queries.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.045,1,.96]);fig.savefig(out/'coordinates.png',dpi=160);plt.close(fig)
(out/'figure-source.json').write_text(json.dumps(dict(sources=sources,protocolSHA256=hashlib.sha256((b/'protocol.json').read_bytes()).hexdigest(),scope='Descriptive coordinates; no inference from apparent group separation, density or distance'),indent=2)+'\n')
