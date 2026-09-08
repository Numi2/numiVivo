"""Input and independent observation controls; no GPU or MD result is mocked as evidence."""
import copy
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from ensemble_campaign import digest, read, write, validate_policy, rigid_triangle_degrees_of_freedom, verify_series, campaign_grid, temperature_seeds, verify_runtime

TOOLS=Path(__file__).parent


def fixture():
    policy=read(TOOLS/"ensemble_policy.json");dt=policy["timeStepsPS"][0]
    def vector(x=0,y=0,z=0):return dict(x=x,y=y,z=z)
    xyz=[vector(),vector(.125),vector(0,.125)];zero=[vector() for _ in xyz]
    cell=dict(a=vector(1),b=vector(0,1),c=vector(0,0,1))
    system=dict(particles=[dict(index=i,role="atom",massDa=1) for i in range(3)],constraints=[
        dict(a=0,b=1,distanceNM=.125),dict(a=0,b=2,distanceNM=.125),dict(a=1,b=2,distanceNM=math.sqrt(2)*.125)])
    config=dict(ensemble="nvt",thermostat="langevinMiddle",barostat="none",positionPrecision="compensated",
        targetTemperatureK=300,frictionPerPS=policy["frictionPerPS"],randomSeed=policy["seeds"][0],timeStepPS=dt,
        maximumConstraintIterations=policy["maximumConstraintIterations"],constraintTolerance=policy["constraintTolerance"],cutoffNM=.4)
    geometry=dict(positionsNM=xyz,periodicCell=cell)
    reference=dict(identifier="exact",geometry=geometry,forcesKJPerMolNM=zero,energyKJPerMol=0)
    request=dict(identifier="water-orthogonal-smooth",system=system,configuration=config,references=[reference],
        dynamicsSteps=round(policy["durationPS"]/dt),dynamicsObserveEvery=round(policy["observationIntervalPS"]/dt),
        dynamicsPreparation="projectConstraints",limits=dict(energyAbsolutePerParticleKJPerMol=.002,
        forceNormalizedRMS=.001,forceNormalizedMaximum=.01,forceNormalizationFloor=1))
    checkpoint=dict(positionsNM=xyz,positionHighNM=xyz,positionCorrectionsNM=zero,
        velocitiesNMPerPS=[vector(1) for _ in xyz],positionPrecision="compensated",periodicCell=cell,
        systemFingerprint="system",configurationFingerprint="configuration",numericalContract="test-only",acceptedStep=0,timePS=0)
    samples=[dict(stepIndex=s,timePS=s*dt,systemFingerprint="system",configurationFingerprint="configuration",
        degreesOfFreedom=6,kineticEnergyKJPerMol=1.5,potentialEnergyKJPerMol=0,totalEnergyKJPerMol=1.5,
        temperatureK=3/(6*policy["boltzmannKJPerMolK"])) for s in range(0,request["dynamicsSteps"]+1,request["dynamicsObserveEvery"])]
    final=copy.deepcopy(checkpoint);final.update(acceptedStep=samples[-1]["stepIndex"],timePS=samples[-1]["timePS"])
    report=dict(outcome="passed",numericalContract="test-only",evaluations=[dict(evaluatedGeometry=geometry,
        physicalParticleForcesKJPerMolNM=zero,energyKJPerMol=0)],comparisons=[dict(identifier="exact",outcome="passed",
        energyErrorPerParticleKJPerMol=0,forceNormalizedRMS=0,forceNormalizedMaximum=0)],
        dynamics=dict(requestedSteps=request["dynamicsSteps"],committedSteps=request["dynamicsSteps"],rejected=False,
        preparedCheckpoint=checkpoint,finalCheckpoint=final,observations=samples,start=samples[0],end=samples[-1]))
    return request,report,policy


class EnsembleCampaignTests(unittest.TestCase):
    def test_missing_or_changed_shader_bundle_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);binary=root/"numivivo";binary.write_text("test fixture, not executable")
            bundle=root/"NumiVivo_NumiVivoShaders.bundle";bundle.mkdir();shader=bundle/"fixture.metal";shader.write_text("test shader")
            manifest=dict(files=[dict(path=str(p.relative_to(root)),sha256=digest(p)) for p in (binary,shader)])
            verify_runtime(binary,manifest)
            with self.assertRaises(ValueError):verify_runtime(binary,dict(files=manifest["files"][:1]))
            shader.write_text("changed")
            with self.assertRaises(ValueError):verify_runtime(binary,manifest)

    def test_temperature_random_streams_are_separate_and_time_steps_are_paired(self):
        policy=read(TOOLS/"ensemble_policy.json");validate_policy(policy)
        low,high=(set(temperature_seeds(policy,t)) for t in policy["temperaturesK"])
        self.assertFalse(low&high);self.assertEqual(len(campaign_grid(policy)),12)
        for temperature in policy["temperaturesK"]:
            sets=[{s for d,t,s in campaign_grid(policy) if d==dt and t==temperature} for dt in policy["timeStepsPS"]]
            self.assertEqual(*sets)
        policy["higherTemperatureSeedOffset"]=policy["seeds"][1]-policy["seeds"][0]
        with self.assertRaises(ValueError):validate_policy(policy)

    def test_prescribed_schedule_and_independent_degrees(self):
        request,report,policy=fixture();validate_policy(policy)
        result=verify_series(request,report,policy,300,policy["seeds"][0],.001)
        self.assertEqual(result["degreesOfFreedom"],6)
        self.assertEqual(len(result["kinetic"]),800)
        self.assertEqual(result["kinetic"],[1.5]*800)

    def test_unqualified_constraint_graphs_refused(self):
        system=fixture()[0]["system"]
        variants=[]
        s=copy.deepcopy(system);s["constraints"].pop();variants.append(s)
        s=copy.deepcopy(system);s["constraints"].append(s["constraints"][0]);variants.append(s)
        s=copy.deepcopy(system);s["particles"][0]["massDa"]=0;variants.append(s)
        s=copy.deepcopy(system);s["constraints"][2]["distanceNM"]=.25;variants.append(s)
        for s in variants:
            with self.subTest(system=s),self.assertRaises(ValueError):rigid_triangle_degrees_of_freedom(s)

    def test_series_corruptions_rejected(self):
        mutations=[lambda r:r["dynamics"]["observations"][200].update(degreesOfFreedom=9),
            lambda r:r["dynamics"]["observations"][200].update(timePS=math.nan),
            lambda r:r["dynamics"]["observations"][200].update(stepIndex=1),
            lambda r:r["dynamics"]["observations"].pop(200),
            lambda r:r["dynamics"]["observations"][200].update(totalEnergyKJPerMol=2),
            lambda r:r["dynamics"]["finalCheckpoint"].update(configurationFingerprint="other"),
            lambda r:r["dynamics"]["finalCheckpoint"].update(numericalContract="other"),
            lambda r:r["dynamics"].update(committedSteps=1),
            lambda r:r["dynamics"]["finalCheckpoint"]["velocitiesNMPerPS"][0].update(x=10)]
        for index,mutation in enumerate(mutations):
            request,report,policy=fixture();mutation(report)
            with self.subTest(mutation=index),self.assertRaises((ValueError,AssertionError)):
                verify_series(request,report,policy,300,policy["seeds"][0],.001)

    def test_prescribed_controls_cannot_silently_change(self):
        for field,value in [("barostat","monteCarlo"),("frictionPerPS",1),("randomSeed",42),("constraintTolerance",.001)]:
            request,report,policy=fixture();request["configuration"][field]=value
            with self.subTest(field=field),self.assertRaises(ValueError):verify_series(request,report,policy,300,policy["seeds"][0],.001)

    def test_hamiltonian_identity_catches_changed_physics(self):
        request,report,policy=fixture()
        a=verify_series(request,report,policy,300,policy["seeds"][0],.001)
        request["configuration"]["cutoffNM"]*=.9
        b=verify_series(request,report,policy,300,policy["seeds"][0],.001)
        self.assertNotEqual(a["hamiltonianIdentity"],b["hamiltonianIdentity"])

    def test_invalid_or_unachievable_policies_refused(self):
        for field,value in [("seeds",[1,1,2]),("cases",["../water"]),("durationPS",1000),
                            ("bootstrapBlockPS",.073),("minimumBlocksPerReplica",1000),
                            ("bootstrapSeed",True),("frictionPerPS",math.nan),("equilibrationPS",10.025)]:
            p=read(TOOLS/"ensemble_policy.json");p[field]=value
            with self.subTest(field=field),self.assertRaises(ValueError):validate_policy(p)

    def test_derivation_retains_failures_physics_and_source(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/"parent";source.mkdir();request=fixture()[0]
            cases=[]
            for name in ("water-orthogonal-smooth","excluded"):
                directory=source/name;directory.mkdir();write(directory/"request.json",request)
                (directory/"reference-system.xml").write_text("<independent-test-fixture/>\n")
                cases.append(dict(identifier=name,status="prepared",requestSHA256=digest(directory/"request.json"),
                                  serializedSystemSHA256=digest(directory/"reference-system.xml")))
            failure=dict(identifier="original-water",status="preparation-failed",error="retained fixture failure")
            cases.append(failure);write(source/"manifest.json",dict(schema="numivivo.org/md-reference-campaign/v1",cases=cases))
            before={str(p.relative_to(source)):digest(p) for p in source.rglob("*") if p.is_file()}
            command=[sys.executable,str(TOOLS/"derive_campaign.py"),"--references",str(source),"--out",str(root/"derived"),
                     "--cases","water-orthogonal-smooth","--ensemble","nvt","--temperature-k","305","--friction-per-ps","10",
                     "--seed","31415","--constraint-iterations","64","--steps","100000","--observe-every","100"]
            result=subprocess.run(command,capture_output=True,text=True);self.assertEqual(result.returncode,0,result.stderr)
            manifest=read(root/"derived/manifest.json");self.assertEqual(manifest["cases"][-1],failure)
            self.assertEqual([c["identifier"] for c in manifest["excludedPreparedCases"]],["excluded"])
            derived=read(root/"derived/water-orthogonal-smooth/request.json")
            for key in ("system","references","limits"):self.assertEqual(request[key],derived[key])
            self.assertEqual(derived["configuration"]["randomSeed"],31415)
            self.assertEqual(derived["configuration"]["targetTemperatureK"],305)
            self.assertEqual(before,{str(p.relative_to(source)):digest(p) for p in source.rglob("*") if p.is_file()})
            self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0,"must not overwrite evidence")
            command[command.index(str(root/"derived"))]=str(root/"invalid")
            command+= ["--seed","-1"]
            self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0)
            self.assertFalse((root/"invalid").exists())


if __name__=="__main__":unittest.main()
