"""Prepare source-bound counts while discarding the query annotation field."""
import argparse, csv, datetime, gzip, hashlib, io, json, shutil, tarfile
from pathlib import Path
import h5py, numpy as np
from scipy import sparse

ROOT = None
RECIPE = Path(__file__).parent
ARCHIVE = None
BINARY = None
MEMBERS = ['GSM2230761_mouse1_umifm_counts.csv.gz', 'GSM2230762_mouse2_umifm_counts.csv.gz']

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v): Path(p).write_text(json.dumps(v,indent=2,sort_keys=True,allow_nan=False)+'\n')
def dataframe(group, columns, index):
    group.attrs['encoding-type']='dataframe'; group.attrs['encoding-version']='0.2.0'
    group.attrs['_index']='_index'; group.attrs['column-order']=np.array(list(columns),dtype=h5py.string_dtype())
    for name,values in {'_index':index,**columns}.items():
        d=group.create_dataset(name,data=np.array(values,dtype=object),dtype=h5py.string_dtype())
        d.attrs['encoding-type']='string-array'; d.attrs['encoding-version']='0.2.0'
def h5ad(path, matrix, features, cells, columns):
    with h5py.File(path,'w') as h:
        h.attrs['encoding-type']='anndata';h.attrs['encoding-version']='0.1.0'
        g=h.create_group('X');g.attrs['encoding-type']='csr_matrix';g.attrs['encoding-version']='0.1.0';g.attrs['shape']=matrix.shape
        for name,a in [('data',matrix.data),('indices',matrix.indices),('indptr',matrix.indptr)]:g.create_dataset(name,data=a,compression='gzip')
        dataframe(h.create_group('obs'),columns,cells);dataframe(h.create_group('var'),{},features)
        for name in ['obsm','varm','obsp','varp','layers','uns']:
            g=h.create_group(name);g.attrs['encoding-type']='dict';g.attrs['encoding-version']='0.1.0'

def main():
    inputs=ROOT/'inputs';inputs.mkdir(exist_ok=False)
    source=ROOT/'source';source.mkdir(exist_ok=False)
    shutil.copy2(RECIPE/'PROTOCOL.md',ROOT/'PROTOCOL.md')
    # Freeze design, exact source and implementation before extracting any labels.
    write(ROOT/'protocol-freeze.json',dict(createdUTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        protocolSHA256=sha(ROOT/'PROTOCOL.md'),sourceArchiveSHA256=sha(ARCHIVE),sourceArchive=str(ARCHIVE),
        members=MEMBERS,implementationSHA256=sha(BINARY),implementation=str(BINARY),
        preparerSHA256=sha(__file__),testLabelsRead=False,fitStarted=False))
    shutil.copy2(BINARY,ROOT/'numivivo-omics')
    records=[];original_features=None
    with tarfile.open(ARCHIVE) as tar:
        for i,member in enumerate(MEMBERS):
            raw=tar.extractfile(member).read();(source/member).write_bytes(raw)
            sample='mouse'+str(i+1);role='training' if i==0 else 'query'
            data=[];indices=[];indptr=[0];cells=[];barcodes=[];labels=[]
            with gzip.GzipFile(fileobj=io.BytesIO(raw)) as binary:
                reader=csv.reader(io.TextIOWrapper(binary));header=next(reader);features=header[3:]
                assert len(features)==14878 and len(set(features))==14878
                if original_features is None: original_features=features
                assert features==original_features
                for row in reader:
                    assert len(row)==len(header)
                    cells.append(sample+':'+row[0]);barcodes.append(row[1])
                    if i==0: labels.append(row[2])
                    values=np.asarray(row[3:],dtype=np.float64)
                    assert np.all(np.isfinite(values)&(values>=0)&(values==np.floor(values)))
                    nz=np.flatnonzero(values);indices.extend(nz.tolist());data.extend(values[nz].astype(np.uint64).tolist());indptr.append(len(data))
            assert len(cells)==[822,1064][i] and len(cells)==len(set(cells))
            matrix=sparse.csr_matrix((np.asarray(data,dtype=np.uint64),np.asarray(indices,dtype=np.int32),np.asarray(indptr,dtype=np.int64)),shape=(len(cells),len(features)))
            columns={'native_sample':[sample]*len(cells),'source_barcode':barcodes}
            if i==0: columns['cell_type']=labels
            h5ad(inputs/(role+'.h5ad'),matrix,features,cells,columns)
            mapping=dict(schemaVersion=1,id='baron2016-mouse-'+role,evidence='measured',sourceDescription='Original GEO '+member+'; all cells and features; query annotation discarded before inference.',countUnit='umiCount',matrixPath='X',sampleColumn='native_sample',mitochondrialFeatureIDs=[],samples=[dict(id=sample,donorID=sample,biologicalReplicateID=sample,batchID='unreported',condition='source pancreatic islet sample; strain/sample confounded',organism='NCBITaxon:10090')])
            if i==0:mapping['groupColumn']='cell_type'
            plan=dict(schemaVersion=1,mapping=mapping,featureNamespace='BARON2016_MOUSE_SOURCE_GENE_SYMBOL')
            if i==0:
                plan.update(reduction=dict(maximumCacheBytes=2000000000,maximumEntryVisits=2000000000,normalizationTarget=10000,pca=dict(components=20,highlyVariableFeatures=2000,maximumBasis=128,meanBins=20,relativeResidualTolerance=1e-6,seed=7)),labelProvenance='Exact source author assigned_cluster from training GSM2230761 only',neighbors=15)
                write(inputs/'fit-knn.json',plan)
                plan['logistic']=dict(penalty=1,gradientTolerance=1e-8,maximumIterations=2000,maximumWork=20000000000)
                write(inputs/'fit-logistic.json',plan)
            else:write(inputs/'query.json',plan)
            records.append(dict(member=member,sha256=hashlib.sha256(raw).hexdigest(),role=role,cells=len(cells),features=len(features),nonzeros=matrix.nnz,totalCounts=int(matrix.sum()),queryLabelsRetained=False))
    write(ROOT/'input-freeze.json',dict(createdUTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),protocolFreezeSHA256=sha(ROOT/'protocol-freeze.json'),sources=records,files={p.name:sha(p) for p in inputs.iterdir()},testLabelsRead=False,fitStarted=False))
    print(json.dumps({'sources':records,'inputFreezeSHA256':sha(ROOT/'input-freeze.json')}))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,required=True);p.add_argument('--archive',type=Path,required=True);p.add_argument('--binary',type=Path,required=True)
    a=p.parse_args();ROOT=a.study.resolve();ARCHIVE=a.archive.resolve();BINARY=a.binary.resolve()
    if ROOT.exists():raise FileExistsError('study must be new; frozen studies are immutable')
    ROOT.mkdir(parents=True);main()
