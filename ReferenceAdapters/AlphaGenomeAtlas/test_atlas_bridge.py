import hashlib, json, tempfile, unittest
from pathlib import Path
import atlas_bridge

class AtlasBridgeTests(unittest.TestCase):
    def fixture(self):
        td = tempfile.TemporaryDirectory(); root = Path(td.name)
        fasta = root/'ref.fa'; fasta.write_bytes(b'>chr1\nACGT\n')
        candidate='a'*64
        request={
            'schema':'numivivo.org/atlas-request/v1',
            'binding':{
                'parentReportSHA256':'b'*64,'caseFingerprint':'c'*64,'dataClass':'publicReference',
                'assembly':'GRCh38','referenceSHA256':hashlib.sha256(fasta.read_bytes()).hexdigest(),
                'hlaAlleles':['HLA-A*02:01'],'tumorRNASampleID':'rna-1'},
            'requestedScorers':['demo'], 'ontologyTerms':['CL:0000000'],
            'variants':[{'key':'GRCh38|chr1|2|C|T','chromosome':'chr1','position1':2,'reference':'C','alternate':'T','candidateIDs':[candidate]}],
            'excludedCandidates':{}
        }
        terms=b'terms'; usage={'purpose':'noncommercialResearch','termsURI':atlas_bridge.TERMS_URI,
            'termsSHA256':hashlib.sha256(terms).hexdigest(),'acknowledgedBy':'tester',
            'externalSharingApproved':True,'snapshotStorageApproved':True,'modelTrainingAllowed':False}
        return td, fasta, request, terms, usage

    def test_fetch_preserves_explicit_status_and_tables(self):
        td,fasta,request,terms,usage=self.fixture(); self.addCleanup(td.cleanup)
        def call(payload, timeout):
            if payload['operation']=='catalog':
                return {'clientVersion':'0.0.0','clientCommit':atlas_bridge.SDK_COMMIT,'scorers':{'demo':{'isSigned':False,'tracks':[]}}}
            q=payload['query']
            return {'clientVersion':'0.0.0','clientCommit':atlas_bridge.SDK_COMMIT,'tables':[
                {'scorer':'demo','isSigned':False,'observations':[{'variantKey':q['key'],'metadata':{'gene_id':'ENSG1'}}],
                 'tracks':[{'ontology_curie':'CL:0000000'}],'rawScores':[[0.5]],'quantiles':[[0.8]]}]}
        data=atlas_bridge.encoded(request)
        result=atlas_bridge.fetch(data,fasta,usage,terms,10,call=call)
        self.assertEqual(result['outcomes'][0]['status'],'available')
        self.assertEqual(result['tables'][0]['rawScores'],[[0.5]])
        self.assertIsNone(result['serviceVersion'])
        self.assertFalse(result['trainingPermitted'])

    def test_reference_mismatch_is_failure_not_zero(self):
        td,fasta,request,terms,usage=self.fixture(); self.addCleanup(td.cleanup)
        request['variants'][0]['reference']='A'; request['variants'][0]['key']='GRCh38|chr1|2|A|T'
        calls=[]
        def call(payload, timeout):
            calls.append(payload['operation'])
            return {'clientVersion':'0.0.0','clientCommit':atlas_bridge.SDK_COMMIT,'scorers':{'demo':{}}}
        result=atlas_bridge.fetch(atlas_bridge.encoded(request),fasta,usage,terms,10,call=call)
        self.assertEqual(result['outcomes'][0],{'variantKey':'GRCh38|chr1|2|A|T','status':'failed','diagnostic':'referenceMismatch'})
        self.assertEqual(result['tables'],[])
        self.assertEqual(calls,['catalog'])

    def test_wrong_reference_bytes_reject_entire_batch(self):
        td,fasta,request,terms,usage=self.fixture(); self.addCleanup(td.cleanup)
        request['binding']['referenceSHA256']='0'*64
        with self.assertRaises(ValueError):
            atlas_bridge.fetch(atlas_bridge.encoded(request),fasta,usage,terms,10,call=lambda *_: {})

    def test_usage_requires_matching_terms_and_no_training(self):
        td,fasta,request,terms,usage=self.fixture(); self.addCleanup(td.cleanup)
        usage['modelTrainingAllowed']=True
        with self.assertRaises(ValueError): atlas_bridge.validate_usage(usage,terms)
        usage['modelTrainingAllowed']=False; usage['termsSHA256']='0'*64
        with self.assertRaises(ValueError): atlas_bridge.validate_usage(usage,terms)

    def test_request_rejects_duplicate_candidate_mapping(self):
        td,fasta,request,terms,usage=self.fixture(); self.addCleanup(td.cleanup)
        request['variants'].append(dict(request['variants'][0]))
        with self.assertRaises(ValueError): atlas_bridge.validate_request(request)

if __name__=='__main__': unittest.main()
