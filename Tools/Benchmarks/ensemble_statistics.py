"""Independent, bounded ensemble diagnostics for retained MD observations.

No native residual, fitted temperature, or iid interpretation of a correlated
trajectory supplies the reference distribution. Failed and insufficient evidence
are distinct. These necessary tests do not establish general equilibration.
"""
import math
import numpy as np
from scipy.special import expit, gammainc


def traces_array(traces):
    values=np.asarray(traces,dtype=np.float64)
    if values.ndim!=2 or values.shape[0]<3 or values.shape[1]<40 or not np.isfinite(values).all():
        raise ValueError("at least three equal-length finite replica traces are required")
    return values


def effective_samples(values):
    """Geyer-style initial positive, monotone pairs; never exceeds raw count."""
    x=np.asarray(values,dtype=np.float64);n=len(x);x=x-x.mean()
    variance=float(np.dot(x,x))
    if not variance>0:return 0.0
    size=1 << (2*n-1).bit_length()
    spectrum=np.fft.rfft(x,n=size)
    acf=np.fft.irfft(spectrum*np.conjugate(spectrum),n=size)[:n]/variance
    total=0.0;previous=math.inf
    for lag in range(1,n//2-1,2):
        pair=float(acf[lag]+acf[lag+1])
        if pair<=0:break
        pair=min(previous,pair);total+=pair;previous=pair
    return float(n/max(1.0,1+2*total))


def split_rhat(values):
    n=values.shape[1]//2
    split=np.concatenate((values[:,:n],values[:,-n:]),axis=0)
    within=float(np.var(split,axis=1,ddof=1).mean())
    if not within>0:return None
    between=n*float(np.var(split.mean(axis=1),ddof=1))
    return float(max(1.0,math.sqrt(((n-1)*within+between)/(n*within))))


def block_values(values,policy):
    size=round(policy["bootstrapBlockPS"]/policy["observationIntervalPS"])
    if size<1 or not math.isclose(size*policy["observationIntervalPS"],policy["bootstrapBlockPS"],abs_tol=1e-12):
        raise ValueError("bootstrap blocks must use an integer observation count")
    if values.shape[1]%size:raise ValueError("bootstrap must retain every production observation")
    return values.reshape(values.shape[0],-1,size)


def quality(traces,policy):
    values=traces_array(traces);blocks=block_values(values,policy).mean(axis=2)
    ess=[effective_samples(row) for row in values];rhat=split_rhat(values)
    centered=blocks-blocks.mean(axis=1,keepdims=True)
    denominator=float(np.square(centered).sum())
    lag=float((centered[:,:-1]*centered[:,1:]).sum()/denominator) if denominator>0 else None
    half=blocks.shape[1]//2
    first=blocks[:,:half].ravel();last=blocks[:,-half:].ravel()
    se=math.sqrt(float(np.var(first,ddof=1)/len(first)+np.var(last,ddof=1)/len(last)))
    drift=float(abs(first.mean()-last.mean())/se) if se>0 else None
    reasons=[]
    if min(ess)<policy["minimumEffectiveSamplesPerReplica"]:reasons.append("insufficient effective samples in a replica")
    if sum(ess)<policy["minimumTotalEffectiveSamples"]:reasons.append("insufficient total effective samples")
    if blocks.shape[1]<policy["minimumBlocksPerReplica"]:reasons.append("insufficient independent blocks")
    if lag is None or abs(lag)>policy["maximumBlockAutocorrelation"]:reasons.append("unresolved correlation between bootstrap blocks")
    if rhat is None or rhat>policy["maximumSplitRHat"]:reasons.append("replica/half-trace disagreement")
    if drift is None or drift>policy["maximumStatisticalZ"]:reasons.append("unresolved first-half/last-half drift")
    return dict(effectiveSamplesPerReplica=ess,totalEffectiveSamples=sum(ess),splitRHat=rhat,
                blocksPerReplica=int(blocks.shape[1]),blockLagOneCorrelation=lag,halfTraceDriftZ=drift,
                inconclusiveReasons=reasons)


def bootstrap_weights(blocks,policy,offset):
    count=policy["bootstrapReplicates"]
    if not isinstance(count,int) or not 200<=count<=10000:raise ValueError("bounded bootstrap replicate count")
    rng=np.random.default_rng(policy["bootstrapSeed"]+offset)
    replicas,per_replica,_=blocks.shape
    # Each replica contributes exactly its original number of full blocks.
    return np.concatenate([rng.multinomial(per_replica,np.full(per_replica,1/per_replica),size=count)
                           for _ in range(replicas)],axis=1)


def kinetic_assessment(traces,degrees_of_freedom,temperature,policy,seed_offset=0):
    values=traces_array(traces)
    if not isinstance(degrees_of_freedom,int) or degrees_of_freedom<=0 or temperature<=0 or np.any(values<0):
        raise ValueError("positive prescribed degrees of freedom, temperature and kinetic energies required")
    k=policy["boltzmannKJPerMolK"];shape=degrees_of_freedom/2;scale=k*temperature
    blocks=block_values(values,policy);weights=bootstrap_weights(blocks,policy,seed_offset)
    flat=values.ravel();n=len(flat);matrix=blocks.reshape(-1,blocks.shape[-1])
    mean=float(flat.mean());variance=float(np.var(flat,ddof=1))
    bootstrap_mean=weights @ matrix.sum(axis=1)/n
    bootstrap_variance=np.maximum(0,(weights @ np.square(matrix).sum(axis=1)-n*np.square(bootstrap_mean))/(n-1))
    tmean=mean/(shape*k);twidth=math.sqrt(variance/shape)/k
    mean_se=float(np.std(bootstrap_mean/(shape*k),ddof=1))
    width_se=float(np.std(np.sqrt(bootstrap_variance/shape)/k,ddof=1))
    # Kolmogorov distance to the prescribed gamma CDF. The null fluctuation is
    # a centered empirical-process block bootstrap, including tied observations.
    order=np.argsort(flat,kind="stable");sorted_values=flat[order]
    unique,first,counts=np.unique(sorted_values,return_index=True,return_counts=True)
    cumulative=np.cumsum(counts);cdf=gammainc(shape,unique/scale)
    distance=float(max(np.max(np.abs(cumulative/n-cdf)),np.max(np.abs((cumulative-counts)/n-cdf))))
    sample_blocks=np.repeat(np.arange(matrix.shape[0]),matrix.shape[1])[order]
    fluctuations=[]
    for row in weights:
        weighted=np.add.reduceat(row[sample_blocks],first)
        fluctuations.append(float(np.max(np.abs(np.cumsum(weighted)-cumulative))/n))
    distribution_p=(1+sum(d>=distance for d in fluctuations))/(1+len(fluctuations))
    mean_z=abs(tmean-temperature)/mean_se if mean_se>0 else None
    width_z=abs(twidth-temperature)/width_se if width_se>0 else None
    diagnostics={"kineticEnergy":quality(values,policy),"squaredKineticDeviation":quality(np.square(values-shape*scale),policy)}
    inconclusive=[key+": "+reason for key,q in diagnostics.items() for reason in q["inconclusiveReasons"]]
    z=policy["maximumStatisticalZ"]
    if mean_se<=0 or z*mean_se/temperature>policy["maximumRelativeMeanTemperatureError"]:inconclusive.append("mean-temperature uncertainty is too large or degenerate")
    if width_se<=0 or z*width_se/temperature>policy["maximumRelativeWidthTemperatureError"]:inconclusive.append("width-temperature uncertainty is too large or degenerate")
    failures=[]
    if abs(tmean/temperature-1)>policy["maximumRelativeMeanTemperatureError"] or (mean_z is not None and mean_z>z):failures.append("kinetic mean misses the prescribed temperature")
    if abs(twidth/temperature-1)>policy["maximumRelativeWidthTemperatureError"] or (width_z is not None and width_z>z):failures.append("kinetic width misses the prescribed gamma distribution")
    if distribution_p<policy["minimumDistributionBootstrapP"]:failures.append("kinetic distribution shape rejected by block-bootstrap gamma CDF check")
    return dict(outcome="inconclusive" if inconclusive else ("failed" if failures else "passed"),
        degreesOfFreedom=degrees_of_freedom,targetTemperatureK=temperature,samples=int(n),replicas=int(values.shape[0]),
        meanKineticEnergyKJPerMol=mean,standardDeviationKineticEnergyKJPerMol=math.sqrt(variance),
        meanTemperatureK=tmean,meanTemperatureStandardErrorK=mean_se,meanTemperatureZ=mean_z,
        widthTemperatureK=twidth,widthTemperatureStandardErrorK=width_se,widthTemperatureZ=width_z,
        gammaCDFDistance=distance,gammaCDFBootstrapP=distribution_p,diagnostics=diagnostics,
        failureReasons=failures,inconclusiveReasons=inconclusive)


def _logistic(x,y,weights,initial=None):
    theta=np.zeros(2) if initial is None else np.array(initial,copy=True)
    for _ in range(60):
        linear=theta[0]+theta[1]*x;probability=expit(linear)
        error=weights*(probability-y);gradient=np.array([error.sum(),np.dot(error,x)])
        if np.max(np.abs(gradient))<=1e-9*weights.sum():return theta
        curvature=weights*probability*(1-probability)
        hessian=np.array([[curvature.sum(),np.dot(curvature,x)],[np.dot(curvature,x),np.dot(curvature,x*x)]])
        if not np.linalg.det(hessian)>1e-12:raise ValueError("unidentifiable temperature-overlap fit")
        step=np.linalg.solve(hessian,gradient)
        cost=float(np.dot(weights,np.logaddexp(0,linear)-y*linear));fraction=1.0
        # Near the optimum, an objective decrease smaller than its FP64 ULP
        # cannot satisfy a useful line search. The Newton decrement supplies
        # a curvature-aware numerical termination test, not a scientific limit.
        if float(np.dot(gradient,step))<=64*np.finfo(np.float64).eps*max(1,abs(cost)):
            return theta-step
        for _ in range(30):
            trial=theta-fraction*step;z=trial[0]+trial[1]*x
            if float(np.dot(weights,np.logaddexp(0,z)-y*z))<=cost-1e-4*fraction*np.dot(gradient,step):
                theta=trial;break
            fraction*=0.5
        else:raise ValueError("temperature-overlap fit line search did not converge")
    raise ValueError("temperature-overlap fit exhausted its iteration bound")


def slope_assessment(low_traces,high_traces,low_temperature,high_temperature,policy,seed_offset=0):
    low=traces_array(low_traces);high=traces_array(high_traces)
    if low.shape!=high.shape or not 0<low_temperature<high_temperature:raise ValueError("matched replicas at two ordered temperatures required")
    low_blocks=block_values(low,policy);high_blocks=block_values(high,policy)
    energy=np.concatenate((low.ravel(),high.ravel()));scale=float(np.std(energy))
    if not scale>0:return dict(outcome="inconclusive",failureReasons=[],inconclusiveReasons=["constant potential energy cannot identify an ensemble slope"])
    x=(energy-energy.mean())/scale;y=np.concatenate((np.zeros(low.size),np.ones(high.size)))
    try:estimate=_logistic(x,y,np.ones_like(x))
    except ValueError as error:return dict(outcome="inconclusive",failureReasons=[],inconclusiveReasons=[str(error)])
    slope=float(estimate[1]/scale)
    probability=expit(estimate[0]+estimate[1]*x)
    overlap=float(np.mean(4*probability*(1-probability)))
    weights_low=bootstrap_weights(low_blocks,policy,seed_offset)
    weights_high=bootstrap_weights(high_blocks,policy,seed_offset+1)
    bootstrap=[];bootstrap_failures=0
    for a,b in zip(weights_low,weights_high):
        weights=np.concatenate((np.repeat(a,low_blocks.shape[-1]),np.repeat(b,high_blocks.shape[-1])))
        try:bootstrap.append(float(_logistic(x,y,weights,estimate)[1]/scale))
        except ValueError:bootstrap_failures+=1
    expected=(1/low_temperature-1/high_temperature)/policy["boltzmannKJPerMolK"]
    se=float(np.std(bootstrap,ddof=1)) if len(bootstrap)>1 else None
    z=abs(slope-expected)/se if se is not None and se>0 else None
    diagnostics={"lowerTemperature":quality(low,policy),"higherTemperature":quality(high,policy)}
    inconclusive=[key+": "+reason for key,q in diagnostics.items() for reason in q["inconclusiveReasons"]]
    if overlap<policy["minimumTemperatureOverlap"]:inconclusive.append("insufficient potential-energy distribution overlap")
    if bootstrap_failures:inconclusive.append("one or more prescribed bootstrap fits failed; none were replaced")
    if se is None or se<=0 or se/expected>policy["maximumRelativeSlopeStandardError"]:inconclusive.append("slope uncertainty is too large or degenerate")
    failures=[]
    if z is not None and z>policy["maximumStatisticalZ"]:failures.append("configurational slope misses the prescribed temperature difference")
    if not failures and se is not None and slope<=policy["maximumStatisticalZ"]*se:inconclusive.append("temperature difference is not statistically resolved")
    return dict(outcome="inconclusive" if inconclusive else ("failed" if failures else "passed"),
        expectedSlopeMolPerKJ=expected,fittedSlopeMolPerKJ=slope,slopeStandardErrorMolPerKJ=se,slopeZ=z,
        estimatedTemperatureDifferenceK=slope*policy["boltzmannKJPerMolK"]*low_temperature*high_temperature,
        overlapScore=overlap,bootstrapFitFailures=bootstrap_failures,diagnostics=diagnostics,
        failureReasons=failures,inconclusiveReasons=inconclusive)
