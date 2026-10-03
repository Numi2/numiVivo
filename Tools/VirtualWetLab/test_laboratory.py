import copy,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from experiment_campaign import expand,verify_barrier,register_and_predict
from qualification import admit_zero_shot
from molecular_adapter import MolecularPerturbationAdapter

class LaboratoryTests(unittest.TestCase):
    def test_matrix_budget_and_determinism(self):
        self.assertEqual(expand({'dose':[0,1],'time':[1,2]}),[{'dose':0,'time':1},{'dose':0,'time':2},{'dose':1,'time':1},{'dose':1,'time':2}])
        with self.assertRaises(ValueError):expand({'dose':list(range(33))})
    def test_incomplete_campaign_cannot_reveal(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(FileNotFoundError):verify_barrier(Path(d))
    def test_arc_rejects_training_like_contexts(self):
        plan={'benchmarkID':'arc-vcc-2026','inputState':'unperturbed','interventionKind':'CRISPRi','trainingProvenanceComplete':True,'trainingPerturbedContexts':['A'],'queryContexts':['B'],'queryUnits':['query'],'trainingUnits':['train'],'heldOutObservationsUsed':False,'modelSHA256':'a'*64,'inputSHA256':'b'*64,'scorerRevision':'c'*40,'queryTargets':['X'],'trainingTargets':['X']}
        self.assertFalse(admit_zero_shot(plan)['zeroShotDemonstrated'])
        for key,value in [('queryContexts',['A']),('trainingProvenanceComplete',False),('queryUnits',['train']),('heldOutObservationsUsed',True),('interventionKind','knockout'),('modelSHA256','unpinned'),('scorerRevision','main'),('queryUnits',[])]:
            with self.assertRaises(ValueError):admit_zero_shot({**plan,key:value})
    def test_compile_preserves_unsupported_boundaries(self):
        a={'id':'assay','specimen':'mouse3','timepoints':['endpoint'],'regions':[{'id':'r1'}],'intervention':{'target':'Clu'},'sourceCitation':'source','files':{'source/report.json':'digest'}}
        selection={'specimen':'mouse3','target':'Clu','regionIDs':['r1'],'timepoint':'endpoint'}
        with patch('molecular_adapter.load',return_value=(a,Path('/tmp'))):
            plan=MolecularPerturbationAdapter().compile('config',selection)
            self.assertEqual([t['state'] for t in plan['transitions']],['HYPOTHESIS','MODEL INFERENCE','UNAVAILABLE','UNAVAILABLE'])
            for s in [{**selection,'dose':1},{**selection,'regionIDs':['missing']},{**selection,'timepoint':'48h'}]:
                with self.assertRaises(ValueError):MolecularPerturbationAdapter().compile('config',s)
    def test_campaign_verifier_rejects_arm_verdict_and_matrix_drift(self):
        from experiment_campaign import verify
        from wetlab import sha,write
        class Owner:
            def summary(self,run):return {'registration':{'plan':{'specimenID':'mouse3'}}}
            def verify(self,*args):return {'status':'verified'}
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);identifier='a'*32;run=root/'arms'/identifier;run.mkdir(parents=True)
            write(run/'seal.json',{})
            write(run/'comparison.json',{'verdict':'criterion-contradicted'})
            registration={'plans':[{'specimenID':'mouse3'}],'request':{'matrix':{'target':['Clu']},'simulationReplicates':1}}
            write(root/'registration.json',registration)
            seal={'registrationSHA256':sha(root/'registration.json'),'arms':[{'run':identifier,'selection':{'target':'Clu'},'replicate':0,'sealSHA256':sha(run/'seal.json')}]}
            write(root/'prediction-seal.json',seal)
            comparison={'predictionSealSHA256':sha(root/'prediction-seal.json'),'arms':[{'run':identifier,'verdict':'criterion-contradicted','comparisonSHA256':sha(run/'comparison.json')}],'criterionPassed':False}
            write(root/'comparison.json',comparison)
            with patch('laboratory.adapter_for_run',return_value=Owner()):
                self.assertEqual(verify(root,{})['status'],'verified')
                comparison['arms'][0]['verdict']='criterion-met-on-development-specimen'
                (root/'comparison.json').write_text(json.dumps(comparison))
                with self.assertRaises(ValueError):verify(root,{})
                seal['arms'][0]['selection']['target']='Srf'
                (root/'prediction-seal.json').write_text(json.dumps(seal))
                with self.assertRaises(ValueError):verify_barrier(root)
    def test_failed_campaign_retains_record_and_blocks_barrier(self):
        class Owner:
            def compile(self,c,s):return {'specimenID':'same-mouse'}
            def predict(self,*args):raise ValueError('native rejection')
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);config=root/'config.json';config.write_text('{}')
            with patch('laboratory.adapter_for_config',return_value=Owner()):
                with self.assertRaises(ValueError):register_and_predict(config,{},root,{'hypothesis':'test','matrix':{'target':['A','B']},'simulationReplicates':1,'successCriterion':'adapter-primary-criterion-in-every-arm'})
            campaign=next(root.glob('campaign-*'));self.assertTrue((campaign/'failure.json').exists());self.assertFalse((campaign/'prediction-seal.json').exists())
if __name__=='__main__':unittest.main()
