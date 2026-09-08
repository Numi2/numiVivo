"""Shared input/observation validation for the independent NVT campaign."""
import json
import hashlib
import itertools
import platform
from importlib import metadata
import math
from pathlib import Path
from run_campaign import digest, verify_comparisons, verify_checkpoint, write


def read(path):return json.loads(Path(path).read_text())


def analysis_environment():
    return dict(python=platform.python_version(),numpy=metadata.version("numpy"),scipy=metadata.version("scipy"))


def verify_runtime(binary,manifest):
    bundle=binary.parent/"NumiVivo_NumiVivoShaders.bundle"
    shaders={str(p.relative_to(binary.parent)) for p in bundle.rglob("*") if p.is_file()}
    paths=[r["path"] for r in manifest["files"]]
    if not shaders or len(set(paths))!=len(paths) or set(paths)!=(shaders|{binary.name}):
        raise ValueError("runtime manifest must cover exactly the executable and complete shader bundle")
    for row in manifest["files"]:
        if digest(binary.parent/row["path"])!=row["sha256"]:raise ValueError("altered executable/shader bundle")


def temperature_seeds(policy,temperature):
    index=policy["temperaturesK"].index(temperature)
    return [seed+index*policy["higherTemperatureSeedOffset"] for seed in policy["seeds"]]


def campaign_grid(policy):
    return [(dt,temperature,seed) for dt,temperature in itertools.product(policy["timeStepsPS"],policy["temperaturesK"])
            for seed in temperature_seeds(policy,temperature)]


def validate_policy(policy):
    if policy.get("schema")!="numivivo.org/md-ensemble-policy/v1":raise ValueError("ensemble policy schema")
    for key in ("durationPS","equilibrationPS","observationIntervalPS","frictionPerPS","constraintTolerance",
                "bootstrapBlockPS","maximumStatisticalZ","minimumEffectiveSamplesPerReplica","minimumTotalEffectiveSamples",
                "maximumBlockAutocorrelation","maximumSplitRHat","minimumDistributionBootstrapP",
                "maximumRelativeMeanTemperatureError","maximumRelativeWidthTemperatureError",
                "maximumRelativeSlopeStandardError","minimumTemperatureOverlap","boltzmannKJPerMolK"):
        if type(policy[key]) not in (int,float) or not math.isfinite(policy[key]) or policy[key]<=0:raise ValueError("invalid policy field "+key)
    for key,lower,upper in [("maximumConstraintIterations",1,16384),("bootstrapReplicates",200,10000),
                            ("minimumBlocksPerReplica",2,1000),("bootstrapSeed",0,2**64-1),
                            ("higherTemperatureSeedOffset",1,2**64-1)]:
        if type(policy[key]) is not int or not lower<=policy[key]<=upper:raise ValueError("invalid policy field "+key)
    if not 0<policy["equilibrationPS"]<policy["durationPS"]:raise ValueError("equilibration must leave production observations")
    if not 1<=len(policy["cases"])<=7 or len(set(policy["cases"]))!=len(policy["cases"]):raise ValueError("case identities")
    if any(not isinstance(c,str) or not c or not all(v.isalnum() or v in "-_" for v in c) for c in policy["cases"]):raise ValueError("unsafe case identifier")
    if not 3<=len(policy["seeds"])<=8 or len(set(policy["seeds"]))!=len(policy["seeds"]):raise ValueError("independent seed identities")
    if any(type(s) is not int or not 0<=s<2**64 for s in policy["seeds"]):raise ValueError("seed must fit UInt64")
    temperatures=policy["temperaturesK"];steps=policy["timeStepsPS"]
    if len(temperatures)!=2 or not 0<temperatures[0]<temperatures[1]<=10000:raise ValueError("two ordered temperatures required")
    all_seeds=[s for t in temperatures for s in temperature_seeds(policy,t)]
    if len(set(all_seeds))!=len(all_seeds) or max(all_seeds)>=2**64:raise ValueError("temperatures require disjoint UInt64 random seeds")
    if len(steps)!=2 or not 0<steps[1]<steps[0]<=0.01:raise ValueError("two ordered time steps required")
    for dt in steps:
        n=round(policy["durationPS"]/dt);interval=round(policy["observationIntervalPS"]/dt)
        if not 0<n<=100000 or interval<=0 or n//interval>1000 or n%interval:raise ValueError("bounded native series capacity")
        for value in (policy["durationPS"],policy["observationIntervalPS"],policy["equilibrationPS"]):
            if not math.isclose(round(value/dt)*dt,value,abs_tol=1e-12):raise ValueError("nonintegral time schedule")
    production=round((policy["durationPS"]-policy["equilibrationPS"])/policy["observationIntervalPS"])
    block=round(policy["bootstrapBlockPS"]/policy["observationIntervalPS"])
    if block<=0 or production%block or not math.isclose(block*policy["observationIntervalPS"],policy["bootstrapBlockPS"],abs_tol=1e-12):raise ValueError("partial production blocks are not allowed")
    if production<40 or production//block<policy["minimumBlocksPerReplica"]:raise ValueError("planned production cannot meet the sample/block minimum")
    if not math.isclose(round(policy["equilibrationPS"]/policy["observationIntervalPS"])*policy["observationIntervalPS"],policy["equilibrationPS"],abs_tol=1e-12):raise ValueError("equilibration must end at an observation boundary")
    for key in ("maximumBlockAutocorrelation","minimumDistributionBootstrapP","maximumRelativeMeanTemperatureError",
                "maximumRelativeWidthTemperatureError","maximumRelativeSlopeStandardError","minimumTemperatureOverlap"):
        if policy[key]>=1:raise ValueError("fractional policy field must be less than one: "+key)


def rigid_triangle_degrees_of_freedom(system):
    """Count independent constraints for disjoint noncollinear rigid triangles.

    This deliberately refuses other graphs instead of assuming every listed
    constraint is independent. Positive masses exclude fixed/virtual particles.
    """
    particles=system["particles"];n=len(particles);adj=[set() for _ in range(n)];lengths={}
    if not n or any(p["index"]!=i or p["role"]!="atom" or not math.isfinite(p["massDa"]) or p["massDa"]<=0 for i,p in enumerate(particles)):
        raise ValueError("rigid-triangle gate requires indexed physical particles with positive masses")
    for c in system["constraints"]:
        a,b=c["a"],c["b"];d=c["distanceNM"];pair=tuple(sorted((a,b)))
        if type(a) is not int or type(b) is not int or not 0<=a<n or not 0<=b<n or a==b or pair in lengths or not math.isfinite(d) or d<=0:raise ValueError("invalid or repeated constraint")
        lengths[pair]=d;adj[a].add(b);adj[b].add(a)
    unseen=set(range(n));molecules=0
    while unseen:
        a=min(unseen);group={a}|adj[a]
        if len(group)!=3 or any(adj[i]!=group-{i} for i in group):raise ValueError("unqualified constraint graph; expected separate rigid triangles")
        a,b,c=sorted(group);distances=sorted([lengths[(a,b)],lengths[(a,c)],lengths[(b,c)]])
        if distances[0]+distances[1]<=distances[2]+1e-10:raise ValueError("degenerate constraint triangle")
        unseen-=group;molecules+=1
    if 3*molecules!=n:raise ValueError("incomplete rigid-molecule partition")
    return 6*molecules


def verify_series(request,report,policy,temperature,seed,dt):
    config=request["configuration"];dof=rigid_triangle_degrees_of_freedom(request["system"])
    expected={"ensemble":"nvt","thermostat":"langevinMiddle","barostat":"none","positionPrecision":"compensated",
              "targetTemperatureK":temperature,"frictionPerPS":policy["frictionPerPS"],"randomSeed":seed,
              "timeStepPS":dt,"maximumConstraintIterations":policy["maximumConstraintIterations"],
              "constraintTolerance":policy["constraintTolerance"]}
    if any(config.get(k)!=v for k,v in expected.items()):raise ValueError("native request differs from prescribed execution controls")
    if request.get("dynamicsPreparation")!="projectConstraints":raise ValueError("missing prescribed initial constraint projection")
    steps=round(policy["durationPS"]/dt);every=round(policy["observationIntervalPS"]/dt)
    if request["dynamicsSteps"]!=steps or request["dynamicsObserveEvery"]!=every:raise ValueError("request schedule mismatch")
    if report["outcome"]!="passed" or report.get("executionError") or not verify_comparisons(request,report):raise ValueError("native/static execution did not pass")
    d=report["dynamics"];samples=d["observations"];first=d["preparedCheckpoint"];last=d["finalCheckpoint"]
    if d["committedSteps"]!=steps or d["requestedSteps"]!=steps or d.get("rejected"):raise ValueError("incomplete native dynamics")
    if len(samples)!=steps//every+1 or samples[0]!=d["start"] or samples[-1]!=d["end"]:raise ValueError("observation boundaries or count mismatch")
    if last["acceptedStep"]-first["acceptedStep"]!=steps:raise ValueError("checkpoint clock mismatch")
    if first["acceptedStep"]!=0 or first["timePS"]!=0:raise ValueError("benchmark did not start from its prescribed initial state")
    for key in ("systemFingerprint","configurationFingerprint"):
        if first[key]!=last[key]:raise ValueError("checkpoint state identity changed")
    if last["timePS"]!=samples[-1]["timePS"]:raise ValueError("final checkpoint/observation time mismatch")
    for checkpoint in (first,last):
        if checkpoint["numericalContract"]!=report["numericalContract"]:raise ValueError("numerical contract mismatch")
        if not verify_checkpoint(request,checkpoint)["passed"]:raise ValueError("independent checkpoint constraints failed")
        if checkpoint["periodicCell"]!=request["references"][0]["geometry"]["periodicCell"]:raise ValueError("NVT cell changed")
    kinetic=[];potential=[];equilibration_steps=round(policy["equilibrationPS"]/dt)
    for i,s in enumerate(samples):
        step=i*every;expected_time=first["timePS"]+step*dt
        if not math.isfinite(s["timePS"]) or s["stepIndex"]!=first["acceptedStep"]+step or abs(s["timePS"]-expected_time)>max(1e-12,step*math.ulp(expected_time)):
            raise ValueError("observation physical clock mismatch")
        if s["degreesOfFreedom"]!=dof:raise ValueError("native degree count disagrees with independent rigid-triangle count")
        if s["configurationFingerprint"]!=first["configurationFingerprint"] or s["systemFingerprint"]!=first["systemFingerprint"]:
            raise ValueError("observation state identity mismatch")
        k=s["kineticEnergyKJPerMol"];u=s["potentialEnergyKJPerMol"];total=s["totalEnergyKJPerMol"]
        if not all(math.isfinite(v) for v in (k,u,total,s["temperatureK"])) or k<0:raise ValueError("invalid observable")
        if not math.isclose(total,k+u,rel_tol=1e-12,abs_tol=1e-9):raise ValueError("energy accounting mismatch")
        if not math.isclose(s["temperatureK"],2*k/(dof*policy["boltzmannKJPerMolK"]),rel_tol=1e-12,abs_tol=1e-9):raise ValueError("temperature conversion mismatch")
        if step>equilibration_steps:kinetic.append(k);potential.append(u)
    invariant=dict(system=request["system"],references=request["references"],limits=request["limits"],
                   configuration={k:v for k,v in config.items() if k not in expected})
    model_hash=hashlib.sha256(json.dumps(invariant,sort_keys=True,allow_nan=False).encode()).hexdigest()
    return dict(degreesOfFreedom=dof,kinetic=kinetic,potential=potential,hamiltonianIdentity=model_hash,prepared={k:first[k] for k in
        ("positionsNM","positionHighNM","positionCorrectionsNM","velocitiesNMPerPS","acceptedStep","timePS")})
