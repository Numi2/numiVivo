from pathlib import Path
import json,anndata
r=Path('/Users/n/numivivo-scverse-chain-20260912');p=Path('/Users/n/numivivo-pca-borrowed-20260912/pca-1');a=anndata.read_h5ad(p/'original.h5ad',backed='r');model=json.loads((p/'model.json').read_text());selected=json.loads((r/'selected-features.json').read_text())
assert set(selected)=={model['features'][i]['featureID'] for i in model['selectedFeatureIndices']}
for folder in ['0-cpu','2-metal']:
 g=json.loads((Path('/Users/n/numivivo-metal-scverse-20260912')/folder/'graph.json').read_text())
 assert [x['barcode'] for x in g['cells']]==a.obs_names.tolist()
 assert [x['sampleID'] for x in g['cells']]==a.obs['native_sample'].astype(str).tolist()
(r/'identity.json').write_text(json.dumps({'allCellBarcodesAndSampleIDsInOrder':True,'independentSelectedFeatureMembershipMatchesNative':True,'cells':len(a.obs),'selectedFeatures':len(selected)},indent=2)+'\n')
