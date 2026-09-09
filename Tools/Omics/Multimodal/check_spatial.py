#!/usr/bin/env python3
"""Complete public Visium spot/count/coordinate interchange; no spatial biology claim."""
import argparse, copy, hashlib, importlib.metadata, json, platform, re, shutil, subprocess, tarfile, time
from pathlib import Path
import anndata as ad
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mudata
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

SOURCE_SHA='eb883e48d9aa935d7959600153ad7ba8c7f1748302fcf56eef4c6999574bb37d'
SPATIAL_SHA='064f7e62e43a705730c30911911940dbbe583ebc6b743dd2dc52eba5037a2ef5'
BASE='https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--source',type=Path,required=True)
p.add_argument('--spatial',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
assert sha(a.source)==SOURCE_SHA and sha(a.spatial)==SPATIAL_SHA
a.out.mkdir(parents=True,exist_ok=False)
inputs=a.out/'input';inputs.mkdir();shutil.copy2(a.source,inputs/'original.h5')
# Read only known regular members; never extract arbitrary archive paths/links.
expected={'spatial/aligned_fiducials.jpg','spatial/tissue_hires_image.png','spatial/detected_tissue_image.jpg','spatial/scalefactors_json.json','spatial/tissue_lowres_image.png','spatial/tissue_positions_list.csv'}
(inputs/'spatial').mkdir()
with tarfile.open(a.spatial) as t:
    assert {m.name for m in t.getmembers() if m.isfile()}==expected
    for name in sorted(expected):
        m=t.getmember(name);assert m.isfile() and m.size<10_000_000
        with t.extractfile(m) as f:(inputs/name).write_bytes(f.read())
with h5py.File(a.source) as f:
    m=f['matrix'];original=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr()
    ids=m['features/id'].asstr()[:];names=m['features/name'].asstr()[:];barcodes=m['barcodes'].asstr()[:]
    assert set(m['features/genome'].asstr()[:])=={'GRCh38'}
assert original.shape==(4039,33538) and original.nnz==21210753
positions=pd.read_csv(inputs/'spatial/tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','y','x'],index_col=0)
assert positions.index.is_unique and set(positions.index[positions.in_tissue==1])==set(barcodes)
xy=positions.loc[barcodes,['x','y']].to_numpy(dtype=np.float64)
scan=sc.read_visium(inputs,count_file='original.h5')
np.testing.assert_array_equal(scan.obs_names,barcodes)
np.testing.assert_array_equal(scan.var['gene_ids'],ids)
np.testing.assert_array_equal(scan.obsm['spatial'],xy)
assert (scan.X!=original).nnz==0
del scan
sample=dict(id='lymph-node',biologicalReplicateID='donor1',donorID='donor1',condition='source-unspecified',batchID='one-section',organism='NCBITaxon:9606')
frame=dict(id='full-resolution-image',unit='pixel',axes=['x','y'],sourceDescription='Space Ranger 1.0 full-resolution image pixels: x=pxl_col_in_fullres, y=pxl_row_in_fullres; no physical-scale conversion')
plan=dict(schemaVersion=1,id='visium-human-lymph-node-1.0',evidence='measured',sourceDescription='10x Genomics Human Lymph Node, Visium, Space Ranger 1.0.0, full filtered release; CC BY 4.0; one section, no independent replication',
    samples=[sample],sampleColumn='sample',barcodeColumn='barcode',defaultObservationKind='spot',observationKindColumn='observation_kind',
    spatial=dict(path='obsm/spatial',frame=frame),assays=[dict(sourceName='rna',id='rna',kind='rna',featureNamespace='Ensembl:GRCh38-SpaceRanger-1.0-source-release',countUnit='umiCount',genomeAssembly='GRCh38',matrixPath='X',featureNameColumn='name')])
write(a.out/'plan.json',plan)
obs=pd.DataFrame(dict(sample=pd.Categorical([sample['id']]*len(barcodes)),barcode=barcodes,observation_kind=pd.Categorical(['spot']*len(barcodes))),index=pd.Index(barcodes))
var=pd.DataFrame(dict(name=names),index=pd.Index(ids))
mu=mudata.MuData({'rna':ad.AnnData(X=original,obs=obs,var=var)})
mu.obs=obs.copy();mu.obsm['spatial']=xy;mu.write(a.out/'input.h5mu');del mu
commands=[]
def run(label,args,success=True):
    start=time.monotonic();r=subprocess.run(['/usr/bin/time','-l',str(a.binary),*map(str,args)],capture_output=True,text=True)
    (a.out/(label+'.log')).write_text(r.stdout+r.stderr)
    peak=re.search(r'(\d+)\s+maximum resident set size',r.stderr)
    commands.append(dict(label=label,arguments=list(map(str,args)),exitCode=r.returncode,expectedSuccess=success,seconds=time.monotonic()-start,maximumResidentBytes=int(peak[1]) if peak else None));write(a.out/'commands.json',commands)
    assert (r.returncode==0)==success,(label,r.stderr[-2500:])
bundle=a.out/'bundle'
run('import',['multiassay-h5mu-import',a.out/'input.h5mu','--plan',a.out/'plan.json','--output',bundle])
run('verify',['multiassay-verify',bundle])
run('repeat',['multiassay-h5mu-import',a.out/'input.h5mu','--plan',a.out/'plan.json','--output',a.out/'repeat'])
for name in ['original.h5','plan.json','dataset.json','dataset.h5mu','receipt.json']:assert sha(bundle/name)==sha(a.out/'repeat'/name),name
native=json.loads((bundle/'dataset.json').read_text())
assert native['samples']==[sample] and native['spatialFrames']==[frame]
assert len(native['observations'])==4039 and all(v['kind']=='spot' for v in native['observations'])
assert [v['identity']['barcode'] for v in native['observations']]==list(barcodes)
np.testing.assert_array_equal([v['position']['coordinates'] for v in native['observations']],xy)
assay=native['assays'][0];d=assay['matrix']
assert [v['id'] for v in assay['features']]==list(ids)
assert [v['name'] for v in assay['features']]==list(names)
counts=sparse.csr_matrix((np.asarray(d['counts'],dtype=np.uint64),d['featureIndices'],d['rowOffsets']),shape=original.shape)
assert (counts!=original).nnz==0
np.testing.assert_array_equal(np.asarray(counts.sum(axis=0)),np.asarray(original.sum(axis=0)))
np.testing.assert_array_equal(np.asarray(counts.sum(axis=1)),np.asarray(original.sum(axis=1)))
export=mudata.read_h5mu(bundle/'dataset.h5mu')
np.testing.assert_array_equal(export.obsm['spatial'],xy)
assert (export.mod['rna'].X!=original).nnz==0
metadata=json.loads(export.uns['numivivo'])
assert metadata['spatialArrays']==[dict(frameID=frame['id'],path='obsm/spatial')]
# Standalone plot checks orientation against the source image; no cell labels inferred.
scale=json.loads((inputs/'spatial/scalefactors_json.json').read_text())['tissue_lowres_scalef']
image=plt.imread(inputs/'spatial/tissue_lowres_image.png');total=np.asarray(original.sum(axis=1)).ravel()
fig,axes=plt.subplots(1,2,figsize=(10,6),constrained_layout=True)
for ax,coords,title in zip(axes,[xy,export.obsm['spatial']],['Source coordinates','Native exported coordinates']):
    ax.imshow(image);scatter=ax.scatter(coords[:,0]*scale,coords[:,1]*scale,c=np.log1p(total),s=5,cmap='viridis',alpha=.65)
    ax.set_title(title);ax.set_axis_off()
fig.colorbar(scatter,ax=axes,label='log(1 + total RNA UMIs)',shrink=.6)
fig.suptitle('Human lymph node — all 4,039 filtered spots; full-resolution pixels scaled for display')
fig.savefig(a.out/'coordinate-overlay.png',dpi=160);plt.close(fig)
del native,counts,assay,d,export
# The native export must be a directly readable spatial input using the same plan.
run('native-reimport',['multiassay-h5mu-import',bundle/'dataset.h5mu','--plan',a.out/'plan.json','--output',a.out/'reimport'])
for name in ['dataset.json','dataset.h5mu']:assert sha(bundle/name)==sha(a.out/'reimport'/name),name
run('reimport-verify',['multiassay-verify',a.out/'reimport'])
# Structural multi-frame control: 2D/3D units, absent positions, and name collision.
fixture=json.loads((Path(__file__).parent/'evidence/2026-09-09/fixture.json').read_text())
fixture['assays'][0]['id']='spatial'
fixture['spatialFrames'].append(dict(id='volume',unit='pixel',axes=['x','y','z'],sourceDescription='Synthetic 3D frame'))
fixture['observations'][0]['position']=dict(frameID='volume',coordinates=[1.,2.,3.])
write(a.out/'frames.json',fixture)
run('multiple-frames',['multiassay-h5mu-write',a.out/'frames.json','--output',a.out/'frames.h5mu'])
frames=mudata.read_h5mu(a.out/'frames.h5mu');meta=json.loads(frames.uns['numivivo'])
assert meta['spatialArrays']==[dict(frameID='slide',path='obsm/spatial_1'),dict(frameID='volume',path='obsm/spatial_2')]
np.testing.assert_array_equal(np.asarray(frames.obsm['spatial']).ravel(),[True,True,True])
np.testing.assert_array_equal(frames.obsm['spatial_1'],[[np.nan,np.nan],[np.nan,np.nan],[2.5,4.]])
np.testing.assert_array_equal(frames.obsm['spatial_2'],[[1.,2.,3.],[np.nan,np.nan,np.nan],[np.nan,np.nan,np.nan]])
write(a.out/'checks.json',dict(status='passed',sourceSHA256=SOURCE_SHA,spatialArchiveSHA256=SPATIAL_SHA,inputH5MUSHA256=sha(a.out/'input.h5mu'),binarySHA256=sha(a.binary),checkerSHA256=sha(Path(__file__)),
    sourceURL=BASE+'filtered_feature_bc_matrix.h5',spatialURL=BASE+'spatial.tar.gz',spots=4039,features=33538,nonzeros=original.nnz,totalUMIs=int(original.sum()),
    allCountsAndCoordinatesExact=True,scanpyCountAndCoordinateAgreement=True,nativeSpatialReimportBytesExact=True,repeatBundleBytesExact=True,multipleFramesAndNameCollisionExact=True,
    commands=len(commands),maximumResidentBytes=max(v['maximumResidentBytes'] or 0 for v in commands),scope='Count/spot/pixel-coordinate interchange; not segmentation, physical image registration, spatial biology or tissue prediction',
    platform=platform.platform(),packages={n:importlib.metadata.version(n) for n in ['numpy','scipy','anndata','mudata','scanpy','h5py','matplotlib']}))
