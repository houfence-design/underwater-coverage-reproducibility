"""Mathematical contract checks; fixtures are not ocean observations."""
from pathlib import Path
import json, numpy as np
from scipy.integrate import quad as integrate
from scipy.stats import truncnorm
from formula_reference import *
from run_experiments import objective, layout, decode, quad
from run_ocean_experiments import conditional_coverage,training_states

ROOT=Path(__file__).resolve().parents[1];checks=[]
def check(name,ok,detail):
    assert ok,(name,detail);checks.append({'check':name,'pass':bool(ok),'detail':detail})

check('Monitoring volume units',domain_volume()==50000.,'100 x 20 x 25 m^3')
for mu,sd,lo,hi in [(120.,40.,120.,270.),(150.,40.,120.,270.),(-1100.,120.,-500.,500.),(1100.,120.,-500.,500.)]:
    mass=integrate(lambda x:float(truncated_density(x,mu,sd,lo,hi)),lo,hi,epsabs=1e-10)[0]
    mean=integrate(lambda x:x*float(truncated_density(x,mu,sd,lo,hi)),lo,hi,epsabs=1e-8)[0]
    a=(lo-mu)/sd;b=(hi-mu)/sd
    expected=truncnorm.mean(a,b,loc=mu,scale=sd)
    check('Normalized density and truncated mean '+str(mu),abs(mass-1)<1e-9 and abs(mean-expected)<1e-7,
          {'integrated_mass':mass,'mean':mean,'reference_mean':float(expected)})

p=np.array([[.9,.8],[.1,.2]]);audit=covariance_audit(p,np.array([.5,.5]))
check('Aggregation order covariance identity',abs(audit['difference']-audit['covariance'])<1e-14 and abs(audit['difference'])>.1,
      {k:float(v) for k,v in audit.items()})
check('Union matches enumerated independent events',abs(pointwise_union([.3,.4])-(.3*.6+.7*.4+.3*.4))<1e-14,'Pointwise independent two-sensor event')

rng=np.random.default_rng(9920);u=rng.random((20000,3))
for centres in [np.array([0.,0.,120.]),np.array([-1100.,1100.,150.])]:
    a=sample_actual(u,centres)
    check('Inverse CDF physical support '+str(centres.tolist()),np.isfinite(a).all() and (a>=PHYSICAL[:,0]).all() and (a<=PHYSICAL[:,1]).all(),
          {'sample_mean':a.mean(axis=0).tolist()})

q,w=quad((9,5,11));s=layout(40,120)
zero=conditional_coverage(s,np.zeros((3,2)),order=32,orders=(9,5,11))
expected=objective(s,order=32,orders=(9,5,11))
check('Zero-current evaluator agrees with base integral',np.max(abs(zero-expected))<1e-12,{'maximum_error':float(np.max(abs(zero-expected)))})
forcing=np.array([[80.,20.],[-80.,-20.]])
daily=conditional_coverage(s,forcing,order=32,orders=(9,5,11))
check('Shared-current conditional symmetry',abs(daily[0]-daily[1])<1e-12,daily.tolist())
check('Velocity displacement units',np.array_equal(current_displacement([.1,-.2],600),[60.,-120.]),'m/s times seconds yields metres')
# A shared latent event shows why unconditional sensor products are invalid.
conditional=np.array([[.9,.9],[.1,.1]])
correct=pointwise_union(conditional).mean();wrong=pointwise_union(conditional.mean(axis=0))
check('Shared forcing requires conditional union before averaging',abs(correct-.59)<1e-14 and abs(wrong-.75)<1e-14,
      {'correct':float(correct),'incorrect_unconditional_product':float(wrong)})
centres,weights=training_states(rng.normal(size=(365,2)),16)
check('Training-state frequency normalization',np.isclose(weights.sum(),1) and centres.shape==(16,2),'365 declared fixture states; no field data generated')
g=coverage_gains(.02591994996,.04343887880)
check('Percent units audit',abs(g['absolute_percentage_points']-1.751892884)<1e-8 and 67<g['relative_coverage_percent']<68,g)
positions=np.array([[[0.,0.,0.]],[[3.,4.,0.]]])
check('Travel proxy units',movement_distance(positions)==5.,'Metres; not an energy measurement')
try:energy_penalized_objective(.5,None,.1)
except ValueError:check('Uncalibrated energy is rejected',True,'Interface cannot manufacture a joule-valued result')
else:raise AssertionError('Uncalibrated energy should be rejected')
ci=block_bootstrap_interval(np.ones(365)*.02)
check('Block bootstrap constant-series identity',np.allclose(ci,[.02,.02]),[float(x) for x in ci])
(ROOT/'data'/'formula_verification.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2))
