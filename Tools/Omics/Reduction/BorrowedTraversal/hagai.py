import pathlib,subprocess,os,time,json,hashlib
r=pathlib.Path('/Users/n/numivivo-pca-borrowed-20260912');env=os.environ.copy();env.update(NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-integration-reference-py/lib/python3.13/site-packages/h5py/.dylibs/libhdf5.320.0.0.dylib',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
results=[]
for i,backend in enumerate(['old','new','new','old','old','new']):
 binary=pathlib.Path('/Users/n/numivivo-native-metal-knn-20260912/build/numivivo-omics') if backend=='old' else r/'build/numivivo-omics'
 output=r/('hagai-pca-'+str(i));cmd=[str(binary),'singlecell-h5ad-pca','/Users/n/numivivo-native-mnn-20260910/hagai/pca/original.h5ad','--plan','/Users/n/numivivo-native-mnn-20260910/hagai/fit.json','--output',str(output)]
 start=time.monotonic()
 with (r/('hagai-run-'+str(i)+'.log')).open('w') as f:p=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
 elapsed=time.monotonic()-start
 matches={name:h(output/name)==h(pathlib.Path('/Users/n/numivivo-native-mnn-20260910/hagai/pca')/name) for name in ['scores.bin','loadings.bin','metadata.json'] } if p.returncode==0 else {}
 results.append(dict(backend=backend,seconds=elapsed,exitCode=p.returncode,binarySHA256=h(binary),matches=matches))
 (r/'hagai-checks.json').write_text(json.dumps(results,indent=2)+'\n')
 assert p.returncode==0 and all(matches.values())
print(json.dumps(results))
