"""Executable reference for every formula family in the manuscript.

The formula registry links these functions and production functions to numbered
equations and pseudocode. This module makes definitions independently testable.
"""
import numpy as np
from scipy.special import ndtr, ndtri
from scipy.stats import t

TARGET = np.array([[-50., 50.], [-10., 10.], [137.5, 162.5]])
PHYSICAL = np.array([[-500., 500.], [-500., 500.], [120., 270.]])

def domain_volume(bounds=TARGET):
    """F01: product of positive side lengths."""
    bounds=np.asarray(bounds)
    assert np.all(bounds[:,1]>bounds[:,0])
    return float(np.prod(bounds[:,1]-bounds[:,0]))

def normal_pdf(x):
    """F02: standard normal density phi; ndtr is its CDF Phi."""
    return np.exp(-np.asarray(x)**2/2)/np.sqrt(2*np.pi)

def truncation_normalizer(mu,sd,lo,hi):
    """F02: mass between physical bounds before conditioning."""
    a=(lo-mu)/sd;b=(hi-mu)/sd
    # Use the survival-tail identity when both bounds are above the mean.
    # This avoids subtracting two floating-point values rounded to one.
    return np.where(a>0,ndtr(-a)-ndtr(-b),ndtr(b)-ndtr(a))

def truncated_density(x,mu,sd,lo,hi):
    """F02: normalized bounded parent-Gaussian marginal."""
    x=np.asarray(x);z=truncation_normalizer(mu,sd,lo,hi)
    return normal_pdf((x-mu)/sd)/(sd*z)*((x>=lo)&(x<=hi))

def truncated_mean(mu,sd,lo,hi):
    """F03: conditional actual mean and mean additive error = mean - mu."""
    z=truncation_normalizer(mu,sd,lo,hi)
    return mu+sd*(normal_pdf((lo-mu)/sd)-normal_pdf((hi-mu)/sd))/z

def joint_density(x,mu,sd,bounds=PHYSICAL):
    """F04: product of coordinate marginals, under stated independence."""
    return np.prod(truncated_density(np.asarray(x),np.asarray(mu),np.asarray(sd),bounds[:,0],bounds[:,1]),axis=-1)

def spherical_hit(actual,target,radius):
    """F05/F13: one target, complete array, exact Euclidean sensing event."""
    return np.any(np.sum((actual-target[...,None,:])**2,axis=-1)<=radius**2,axis=-1)

def pointwise_union(p):
    """F06: product across independent sensors at the SAME target."""
    return 1-np.prod(1-np.asarray(p),axis=-1)

def spatial_average(field,weights):
    """F07: normalized target quadrature (last dimension = targets)."""
    weights=np.asarray(weights)
    assert np.isclose(weights.sum(),1.)
    return np.asarray(field)@weights

def covariance_audit(p,weights):
    """F08: exact two-sensor difference between aggregation orders."""
    assert p.shape[1]==2
    means=weights@p
    correct=weights@pointwise_union(p)
    wrong=pointwise_union(means)
    cov=weights@(p[:,0]*p[:,1])-means[0]*means[1]
    return {'correct':correct,'wrong':wrong,'covariance':cov,'difference':wrong-correct}

def horizontal_noncentrality(q,s,scale):
    """F09: squared horizontal separation divided by parent variance."""
    return np.sum((np.asarray(q)[:2]-np.asarray(s)[:2])**2)/scale**2

def depth_limits(qz,radius,lo=120.,hi=270.):
    """F10: sphere/physical-depth intersection and nonnegative width."""
    lower=max(lo,qz-radius);upper=min(hi,qz+radius)
    return lower,max(lower,upper)

def sample_actual(uniform,centres,sd=np.array([120.,120.,40.]),bounds=PHYSICAL):
    """F14: inverse-CDF sampling, broadcasting over trials and sensors."""
    a=(bounds[:,0]-centres)/sd;b=(bounds[:,1]-centres)/sd
    ordinary=ndtri(ndtr(a)+uniform*(ndtr(b)-ndtr(a)))
    upper_tail=-ndtri(ndtr(-b)+(1-uniform)*(ndtr(-a)-ndtr(-b)))
    return centres+sd*np.where(a>0,upper_tail,ordinary)

def batch_interval(estimates,alpha=.05):
    """F15: interval over independent MC batches, not quadrature nodes."""
    estimates=np.asarray(estimates);n=len(estimates)
    se=estimates.std(ddof=1)/np.sqrt(n)
    half=t.ppf(1-alpha/2,n-1)*se
    return float(estimates.mean()),float(estimates.mean()-half),float(estimates.mean()+half)

def coverage_gains(fixed,optimized):
    """F16: pp gain, relative coverage gain, relative uncovered reduction."""
    return {'absolute_percentage_points':100*(optimized-fixed),
            'relative_coverage_percent':100*(optimized/fixed-1),
            'relative_uncovered_percent':100*(optimized-fixed)/(1-fixed)}

def pso_update(position,velocity,personal,global_best,r1,r2):
    """F17: exact update used by the benchmark."""
    velocity=.7298*velocity+1.49618*r1*(personal-position)+1.49618*r2*(global_best-position)
    velocity=np.clip(velocity,-.2,.2)
    proposed=position+velocity
    outside=(proposed<0)|(proposed>1)
    return np.clip(proposed,0,1),np.where(outside,0,velocity)

def afsa_schedule(evaluations,budget):
    """F18: shrinking visual radius and movement scale."""
    fraction=evaluations/budget
    return .5*(1-fraction)+.05,.1*(1-fraction)+.01

def prescribed_bias(step):
    """F19: synthetic 13-step illustrative forcing, never ocean data."""
    return np.array([80*np.sin(step*np.pi/6),40*np.sin(step*np.pi/3)])

def movement_distance(positions):
    """F20: sum distance over time and nodes; output metres, not joules."""
    return float(np.linalg.norm(np.diff(positions,axis=0),axis=-1).sum())

def energy_penalized_objective(coverage,energy_joules,weight_per_joule):
    """F21: proposed interface ONLY; no calibrated energy values supplied."""
    if energy_joules is None:raise ValueError('A calibrated energy input is required')
    return coverage-weight_per_joule*energy_joules

def current_displacement(velocity_m_s,horizon_s):
    """F22: frozen-current kinematic replay, not measured sensor displacement."""
    return np.asarray(velocity_m_s)*horizon_s

def empirical_current_mean(conditional_coverages,weights=None):
    """F23: average AFTER computing conditional pointwise union and volume."""
    x=np.asarray(conditional_coverages)
    w=np.ones(len(x))/len(x) if weights is None else np.asarray(weights)
    assert np.isclose(w.sum(),1.)
    return float(w@x)

def block_bootstrap_interval(daily_difference,block=14,repetitions=2000,seed=7100):
    """F24: circular moving-block percentile interval for a daily replay mean."""
    x=np.asarray(daily_difference);n=len(x);rng=np.random.default_rng(seed)
    starts=rng.integers(n,size=(repetitions,int(np.ceil(n/block))))
    idx=((starts[...,None]+np.arange(block))%n).reshape(repetitions,-1)[:,:n]
    resampled=x[idx].mean(axis=1)
    return tuple(np.quantile(resampled,[.025,.975]))
