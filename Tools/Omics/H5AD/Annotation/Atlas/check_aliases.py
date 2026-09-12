from pathlib import Path
import h5py,numpy as np,json,hashlib,subprocess,shutil,argparse
parser=argparse.ArgumentParser();parser.add_argument('--binary',type=Path,required=True);parser.add_argument('--seed',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);a=parser.parse_args();r=a.out;r.mkdir(exist_ok=False);source=r/'alias-source.h5ad';shutil.copy2(a.seed,source)
with h5py.File(source,'r+') as f:
 f['uns/obs_alias']=f['obs'];f['uns/sample_alias']=f['obs/mapped_sample'];f['uns/journal_alias']=f['uns/numivivo_edits'];oldEntries=list(f['uns/numivivo_edits']);oldColumns=list(f['obs'].attrs['column-order'])
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();before=sha(source);plan=r/'alias-plan.json';plan.write_text(json.dumps({'schemaVersion':1,'source':{'bytes':list(bytes.fromhex(before))},'provenance':'Verify shallow group detachment preserves aliases and old history','edits':[{'path':'obs/mapped_sample','mode':'replace','value':{'string':{'shape':[7],'values':['new']*7}}},{'path':'obs/added','mode':'add','value':{'int64':{'shape':[7],'values':list(range(7))}}}]}));dest=r/'alias-output.h5ad';x=subprocess.run([str(a.binary),'annotate',str(source),str(plan),str(dest)],capture_output=True,text=True);(r/'alias.log').write_text(x.stdout+x.stderr);assert x.returncode==0;assert sha(source)==before
with h5py.File(dest) as f:
 np.testing.assert_equal(f['obs/mapped_sample'].asstr()[()],['new']*7);np.testing.assert_equal(f['uns/obs_alias/mapped_sample'].asstr()[()],['sample']*7);np.testing.assert_equal(f['uns/sample_alias'].asstr()[()],['sample']*7)
 assert list(f['uns/obs_alias'].attrs['column-order'])==oldColumns;assert list(f['obs'].attrs['column-order'])==oldColumns+['added']
 assert set(f['uns/journal_alias'])==set(oldEntries);assert len(f['uns/numivivo_edits'])==len(oldEntries)+1
 addr=lambda p:h5py.h5o.get_info(f[p].id).addr
 assert addr('obs')!=addr('uns/obs_alias');assert addr('obs/index')==addr('uns/obs_alias/index');assert addr('obs/mapped_sample')!=addr('uns/sample_alias');assert addr('uns/numivivo_edits')!=addr('uns/journal_alias')
(r/'alias-verification.json').write_text(json.dumps({'status':'passed','sourceSHA256':before,'outputSHA256':sha(dest),'originalGroupAliasesUnchanged':True,'untouchedChildStorageShared':True,'oldJournalAliasUnchanged':True},indent=2)+'\n');print('Alias preservation passed')
