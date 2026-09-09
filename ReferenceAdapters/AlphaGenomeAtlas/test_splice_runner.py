import hashlib, json, os, tempfile, unittest
from pathlib import Path
import run_pvacsplice

def h(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def rec(path): return {'path':str(path),'sha256':h(path)}

class SpliceRunnerTests(unittest.TestCase):
    def make(self):
        td=tempfile.TemporaryDirectory(); root=Path(td.name)
        inputs={}
        for name in run_pvacsplice.EXPECTED_INPUTS:
            p=root/name; p.write_text('data-'+name); inputs[name]=rec(p)
        reg=root/'regtools'; reg.write_text('#!/usr/bin/env python3\nimport pathlib,sys\nout=pathlib.Path(sys.argv[sys.argv.index("-o")+1]); out.write_text("chrom\\tstart\\tend\\tname\\tscore\\tstrand\\tanchor\\ttranscripts\\n1\\t1\\t3\\tJ1\\t4\\t+\\tD\\tENST1\\n")\n'); reg.chmod(0o755)
        pvac=root/'pvacsplice'; pvac.write_text('#!/usr/bin/env python3\nimport pathlib,sys\nargs=sys.argv[1:]; out=next(pathlib.Path(x) for x in args if x.endswith("/upstream")); (out/"MHC_Class_I").mkdir(parents=True); (out/"MHC_Class_I"/"S.MHC_I.all_epitopes.tsv").write_text("header\\n"); (out/"S.transcripts.fa").write_text(">ALT.x\\nACDEFGHIK\\n>WT.x\\nACDEYGHIK\\n")\n'); pvac.chmod(0o755)
        pred=root/'iedb'; pred.mkdir(); (pred/'resource').write_text('predictor')
        inv, _=run_pvacsplice.inventory(pred)
        binding={'parentReportSHA256':'a'*64,'caseFingerprint':'b'*64,'dataClass':'publicReference','assembly':'GRCh38',
                 'referenceSHA256':inputs['referenceFASTA']['sha256'],'hlaAlleles':['HLA-A*02:01'],'tumorRNASampleID':'rna-1'}
        job={'schema':'numivivo.org/splice-job/v1','binding':binding,'inputs':inputs,'tools':{'regtools':rec(reg),'pvacsplice':rec(pvac)},
             'resourceFiles':{},'sampleName':'S','normalSampleName':'N','rnaStrand':'RF','predictors':['NetMHC'],'epitopeLengths':[9],
             'threads':1,'junctionMinimumReads':10,'variantDistance':100,'expressionMinimum':1.0,
             'resourceVersions':{'regtools':'1.0.0','pvacsplice':'7.1.2','NetMHC':'4.0'},'sourceCitation':'synthetic public-reference runner test'}
        return td,root,job,pred

    def test_runner_writes_importable_bundle_contract(self):
        td,root,job,pred=self.make(); self.addCleanup(td.cleanup)
        out=root/'out'; result=run_pvacsplice.execute(run_pvacsplice.encoded(job),out,30,pred)
        self.assertEqual(result['status'],'external-splice-evidence-written')
        completion=json.loads((out/'complete.json').read_text())
        self.assertEqual(set(completion['files']),{'splice-manifest.json','splice-report.tsv','regtools.tsv','transcripts.fa','splice-execution.json'})
        for name,digest in completion['files'].items(): self.assertEqual(h(out/name),digest)
        manifest=json.loads((out/'splice-manifest.json').read_text()); receipt=json.loads((out/'splice-execution.json').read_text())
        self.assertEqual(manifest['executionReceiptSHA256'],h(out/'splice-execution.json'))
        self.assertEqual([s['tool'] for s in receipt['steps']],['regtools','pvacsplice'])
        self.assertEqual(receipt['outputSHA256']['splice-report.tsv'],h(out/'splice-report.tsv'))

    def test_input_digest_mismatch_rejected_before_execution(self):
        td,root,job,pred=self.make(); self.addCleanup(td.cleanup)
        job['inputs']['rnaBAM']['sha256']='0'*64
        with self.assertRaises(ValueError): run_pvacsplice.execute(run_pvacsplice.encoded(job),root/'out',30,pred)

    def test_non_public_reference_rejected(self):
        td,root,job,pred=self.make(); self.addCleanup(td.cleanup)
        job['binding']['dataClass']='synthetic'
        with self.assertRaises(ValueError): run_pvacsplice.validate_job(job)

if __name__=='__main__': unittest.main()
