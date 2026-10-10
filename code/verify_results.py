"""Independent event-integral and stored-data consistency checks."""
from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.special import ndtr
from numpy.polynomial.legendre import leggauss
from run_experiments import C,pointwise,decode,quad,objective,layout
R=Path(__file__).resolve().parents[1];D=R/'data'
checks=[]
def verify(name,condition,detail):
    checks.append({'check':name,'pass':bool(condition),'detail':detail})
    assert condition,(name,detail)
q,w=quad((11,7,15));verify('Target weights',abs(w.sum()-1)<1e-12,float(w.sum()))
rng=np.random.default_rng(1289)
for u in rng.random((100,4)):
    s=decode(u);d=np.linalg.norm(s[:,None,:]-s[None,:,:],axis=2)+np.eye(9)*1e6
    assert d.min()>=40-1e-9 and s[:,:2].min()>=-500 and s[:,:2].max()<=500 and s[:,2].min()>=120 and s[:,2].max()<=270
verify('100 random feasible grids',True,'bounds and nominal spacing >=40 m')
# Independent disk quadrature with exact conditional-depth CDF.
def disk(q,s,order=160):
    g,w=leggauss(order);x=C.R*g;hy=np.sqrt(C.R**2-x*x);y=hy[:,None]*g
    hz=np.sqrt(np.maximum(C.R**2-x[:,None]**2-y*y,0))
    xyx=q[0]+x[:,None];xyy=q[1]+y
    norm=np.prod(ndtr((np.array([500,500,270])-s)/[120,120,40])-ndtr((np.array([-500,-500,120])-s)/[120,120,40]))
    fxy=np.exp(-.5*((xyx-s[0])/120)**2-.5*((xyy-s[1])/120)**2)/(2*np.pi*120**2)
    lo=np.maximum(q[2]-hz,120);hi=np.minimum(q[2]+hz,270)
    mass=np.maximum(ndtr((hi-s[2])/40)-ndtr((lo-s[2])/40),0)
    return np.sum(fxy*mass*C.R*w[:,None]*hy[:,None]*w[None,:])/norm
errs=[]
for pos,s in [(np.array([0.,0.,150]),np.array([0.,0.,150.])),(np.array([-50.,10.,137.5]),np.array([40.,40.,120.])),(np.array([50.,-10.,162.5]),np.array([100.,100.,150.]))]:
    a=pointwise(pos[None,:],s[None,:],order=80)[0,0];b=disk(pos,s);errs.append(abs(a-b))
verify('Independent spherical integral',max(errs)<2e-7,errs)
mc=pd.read_csv(D/'mc_batches.csv')
for rep in range(30):
    with np.load(D/'mc_raw'/f'run_{rep:02d}.npz') as z:
        for name in z.files:
            count=int(np.unpackbits(z[name])[:50000].sum());row=mc[(mc.rep==rep)&(mc.layout==name)].iloc[0]
            assert count==row.hits and abs(row.coverage-count/50000)<1e-12
verify('Raw Monte Carlo counts',True,'150 arrays of 50000 Bernoulli outcomes agree with CSV')
run=pd.read_csv(D/'optimizer_runs.csv');verify('Optimization seeds and budgets',len(run)==80 and (run.evaluations<=600).all(),'20 runs each for four methods; common maximum 600')
ss=json.loads((D/'summary.json').read_text());verify('Coverage monotonicity across sensor addition',ss['Fixed']['quadrature']>=ss['Single']['quadrature'],'fixed array contains reference single sensor')
maps=pd.read_csv(D/'coverage_maps.csv');verify('Probability range',maps.coverage.between(0,1).all(),'all archived spatial probabilities within [0,1]')
cv=pd.read_csv(D/'quadrature_convergence.csv');changes={}
for name in cv.layout.unique():
    a=cv[cv.layout==name].coverage.to_numpy();changes[name]=abs(a[-1]-a[-2])
verify('Fine-to-finest refinement',max(changes.values())<1e-6,changes)
(D/'verification.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
