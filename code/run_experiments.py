"""Reproducible synthetic coverage experiments; no field observations.
Run: python code/run_experiments.py. Outputs CSV/JSON and packed IID MC outcomes.
The spherical integral uses an exact horizontal Gaussian disk CDF and GL depth
quadrature. All target spheres lie within the horizontal physical box.
"""
from pathlib import Path
from dataclasses import dataclass, asdict, replace
import csv,json,time,sys,platform
import numpy as np
import scipy
from scipy.special import ndtr,ndtri
from scipy.stats import ncx2,t
from scipy.optimize import minimize
from numpy.polynomial.legendre import leggauss
from formula_reference import pso_update, afsa_schedule

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'data'; OUT.mkdir(exist_ok=True)
@dataclass(frozen=True)
class Config:
    R:float=20.; sh:float=120.; sz:float=40.
    xmin:float=-500.; xmax:float=500.; zmin:float=120.; zmax:float=270.
    dmin:float=40.; dmax:float=500.
C=Config()
def save(name,rows):
    if rows:
        with (OUT/name).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader();w.writerows(rows)
def quad(orders):
    axes=[];ws=[]
    for n,(a,b) in zip(orders,[(-50,50),(-10,10),(137.5,162.5)]):
        g,w=leggauss(n);axes.append((a+b)/2+(b-a)*g/2);ws.append(w/2)
    pts=np.column_stack([a.ravel() for a in np.meshgrid(*axes,indexing='ij')])
    weights=np.prod(np.array(np.meshgrid(*ws,indexing='ij')),axis=0).ravel()
    return pts,weights
def decode(u,c=C):
    ax,ay,d,z=np.asarray(u)*[2,2,c.dmax-c.dmin,c.zmax-c.zmin]+[-1,-1,c.dmin,c.zmin]
    cx=ax*(500-d);cy=ay*(500-d)
    return np.array([[cx+i*d,cy+j*d,z] for i in [-1,0,1] for j in [-1,0,1]])
def layout(d=100,z=150,cx=0,cy=0):
    return np.array([[cx+i*d,cy+j*d,z] for i in [-1,0,1] for j in [-1,0,1]],float)
def pointwise(q,s,c=C,order=24):
    q=np.asarray(q);s=np.asarray(s)
    assert np.all(q[:,:2]-c.R>=c.xmin) and np.all(q[:,:2]+c.R<=c.xmax)
    nc=((q[:,None,:2]-s[None,:,:2])**2).sum(axis=2)/c.sh**2
    normxy=np.prod(ndtr((c.xmax-s[:,:2])/c.sh)-ndtr((c.xmin-s[:,:2])/c.sh),axis=1)
    if c.sz==0:
        r2=c.R**2-(q[:,None,2]-s[None,:,2])**2
        return ncx2.cdf(np.maximum(r2,0)/c.sh**2,2,nc)/normxy
    g,w=leggauss(order)
    lo=np.maximum(q[:,2]-c.R,c.zmin);hi=np.maximum(lo,np.minimum(q[:,2]+c.R,c.zmax))
    z=lo[:,None]+(g+1)[None,:]*(hi-lo)[:,None]/2
    r2=np.maximum(c.R**2-(q[:,None,2]-z)**2,0)/c.sh**2
    nz=ndtr((c.zmax-s[:,2])/c.sz)-ndtr((c.zmin-s[:,2])/c.sz)
    f=np.exp(-.5*((z[:,:,None]-s[None,None,:,2])/c.sz)**2)/(c.sz*np.sqrt(2*np.pi)*nz)
    out=np.sum(ncx2.cdf(r2[:,:,None],2,nc[:,None,:])*f*w[None,:,None],axis=1)*(hi-lo)[:,None]/2/normxy
    return np.clip(out,0,1)
def objective(s,c=C,order=24,orders=(7,5,7)):
    q,w=quad(orders);p=pointwise(q,s,c,order)
    return float(w@(1-np.prod(1-p,axis=1)))
class BudgetEnd(Exception):pass
def optimize(method,seed,c=C,budget=600,evaluator=None,refined_evaluator=None):
    rng=np.random.default_rng(seed);pop=rng.random((20,4));initial=pop[0].copy()
    q,w=quad((5,3,5));trace=[];best=[-1.,None];t0=time.perf_counter()
    def ev(u):
        if len(trace)>=budget:raise BudgetEnd
        if evaluator is None:
            p=pointwise(q,decode(u,c),c,16);v=float(w@(1-np.prod(1-p,axis=1)))
        else:v=float(evaluator(decode(u,c)))
        if v>best[0]:best[:]=[v,np.array(u).copy()]
        trace.append(best[0]);return v
    status='budget'
    try:
        if method=='Powell':
            res=minimize(lambda u:-ev(u),initial,method='Powell',bounds=[(0,1)]*4,
                options={'maxfev':budget,'maxiter':100,'ftol':1e-9,'xtol':1e-5})
            status=str(res.message)
        elif method=='PSO':
            v=np.array([ev(u) for u in pop]);personal=pop.copy();pv=v.copy();vel=np.zeros_like(pop)
            while True:
                gb=personal[np.argmax(pv)].copy()
                pop,vel=pso_update(pop,vel,personal,gb,rng.random(pop.shape),rng.random(pop.shape))
                for i in range(20):
                    v[i]=ev(pop[i])
                    if v[i]>pv[i]:pv[i]=v[i];personal[i]=pop[i].copy()
        elif method=='AFSA':
            # Explicit bounded AFSA variant: following, swarming, prey trials,
            # crowding fraction, shrinking visual/step, and global best bulletin.
            v=np.array([ev(u) for u in pop])
            while True:
                for i in range(20):
                    visual,step=afsa_schedule(len(trace),budget)
                    dist=np.linalg.norm(pop-pop[i],axis=1);nb=np.flatnonzero((dist>0)&(dist<visual))
                    candidates=[]
                    if 0<len(nb)<=.6*20:
                        j=nb[np.argmax(v[nb])]
                        if v[j]>v[i]:candidates.append((v[j],pop[j].copy()))
                        cen=pop[nb].mean(axis=0);vc=ev(cen)
                        if vc>v[i]:candidates.append((vc,cen))
                    for _ in range(3):
                        delta=rng.normal(size=4);delta/=np.linalg.norm(delta)
                        trial=np.clip(pop[i]+visual*rng.random()**.25*delta,0,1);vt=ev(trial)
                        if vt>v[i]:candidates.append((vt,trial));break
                    if candidates:
                        target=max(candidates,key=lambda x:x[0])[1];delta=target-pop[i];norm=np.linalg.norm(delta)
                        nxt=np.clip(pop[i]+min(step*rng.random(),norm)*delta/max(norm,1e-15),0,1)
                    else:nxt=np.clip(pop[i]+step*rng.uniform(-1,1,4),0,1)
                    pop[i]=nxt;v[i]=ev(nxt)
        else:
            while True:ev(rng.random(4))
    except BudgetEnd:pass
    s=decode(best[1],c)
    row={'method':method,'seed':seed,'evaluations':len(trace),'seconds':time.perf_counter()-t0,
         'training':best[0],'test':objective(s,c,48,(11,7,15)) if refined_evaluator is None else float(refined_evaluator(s)),
         'cx':float(s[:,0].mean()),'cy':float(s[:,1].mean()),'pitch':float(s[3,0]-s[0,0]),'z':float(s[0,2]),'status':status}
    row.update({f'initial_u{i}':float(x) for i,x in enumerate(initial)})
    return row,[{'method':method,'seed':seed,'evaluation':i+1,'best':v} for i,v in enumerate(trace)]

def main():
    start=time.perf_counter()
    (OUT/'parameters.json').write_text(json.dumps({'physical':asdict(C),'domain':[[-50,50],[-10,10],[137.5,162.5]],'n':9,'budget':600,'optimizer_seeds':list(range(2100,2120)),'MC_repetitions':30,'MC_trials_per_repetition':50000,'provenance':'new synthetic computational experiments; not field data'},indent=2))
    runs=[];traces=[]
    for method in ['Powell','PSO','AFSA','Random']:
        for seed in range(2100,2120):
            row,tr=optimize(method,seed);runs.append(row);traces+=tr
            save('optimizer_runs.csv',runs);save('optimizer_trace.csv',traces)
            print(method,seed,round(row['test'],8),row['evaluations'],flush=True)
    best=max((r for r in runs if r['method']=='Powell'),key=lambda r:r['training'])
    opt=layout(best['pitch'],best['z'],best['cx'],best['cy'])
    layouts={'Single':np.array([[0.,0.,150.]]),'Fixed':layout(),'Pitch-only':layout(40,150),'Depth-only':layout(100,120),'Optimized':opt}
    conv=[]
    for name,s in layouts.items():
        for order,orders in [(16,(5,3,5)),(24,(7,5,7)),(48,(11,7,15)),(80,(15,9,25))]:
            conv.append({'layout':name,'sphere_order':order,'target_orders':str(orders),'coverage':objective(s,order=order,orders=orders)})
    save('quadrature_convergence.csv',conv)
    # Paired MC uses independent observations and runs; common random numbers
    # couple layouts only, never sensors within the same layout.
    raw=OUT/'mc_raw';raw.mkdir(exist_ok=True);mc=[]
    for rep in range(30):
        rng=np.random.default_rng(31000+rep);n=50000
        q=rng.uniform([-50,-10,137.5],[50,10,162.5],(n,3));u=rng.random((n,9,3));bits={}
        for name,s in layouts.items():
            sd=np.array([120,120,40]);lo=np.array([-500,-500,120]);hi=np.array([500,500,270])
            a=ndtr((lo-s)/sd);b=ndtr((hi-s)/sd);actual=s+sd*ndtri(a+u[:,:len(s)]*(b-a))
            hit=(((actual-q[:,None,:])**2).sum(axis=2)<=400).any(axis=1)
            mc.append({'layout':name,'rep':rep,'seed':31000+rep,'trials':n,'hits':int(hit.sum()),'coverage':float(hit.mean())})
            bits[name]=np.packbits(hit)
        np.savez_compressed(raw/f'run_{rep:02d}.npz',**bits)
    save('mc_batches.csv',mc)
    summary={}
    for name in layouts:
        x=np.array([r['coverage'] for r in mc if r['layout']==name]);half=t.ppf(.975,29)*x.std(ddof=1)/np.sqrt(30)
        summary[name]={'quadrature':next(r['coverage'] for r in conv if r['layout']==name and r['sphere_order']==80),'mc_mean':x.mean(),'mc_sd':x.std(ddof=1),'ci_low':x.mean()-half,'ci_high':x.mean()+half}
    diff=np.array([r['coverage'] for r in mc if r['layout']=='Optimized'])-np.array([r['coverage'] for r in mc if r['layout']=='Fixed'])
    h=t.ppf(.975,29)*diff.std(ddof=1)/np.sqrt(30);summary['paired_gain']={'mean':diff.mean(),'low':diff.mean()-h,'high':diff.mean()+h}
    sens=[]
    for sh in [30.,60.,120.,180.]:
        for sz in [0.,15.,40.,60.]:
            for name in ['Fixed','Optimized']:
                sens.append({'layout':name,'sigma_h':sh,'sigma_z':sz,'coverage':objective(layouts[name],replace(C,sh=sh,sz=sz),48,(11,7,25))})
    save('uncertainty_sensitivity.csv',sens)
    sweep=[]
    for variable,values in [('pitch',np.arange(40,201,10)),('depth',np.arange(120,271,5)),('radius',np.arange(10,41,5)),('x',np.arange(-200,201,10)),('y',np.arange(-200,201,10))]:
        for value in values:
            cfg=replace(C,R=float(value)) if variable=='radius' else C
            for name,s in [('Fixed',layout()),('Optimized',opt)]:
                ss=s.copy()
                if variable=='pitch':ss=layout(value,ss[0,2])
                if variable=='depth':ss[:,2]=value
                if variable=='x':ss[:,0]+=value
                if variable=='y':ss[:,1]+=value
                sweep.append({'variable':variable,'value':float(value),'layout':name,'coverage':objective(ss,cfg,32,(9,5,11))})
    save('parameter_sweeps.csv',sweep)
    rng=np.random.default_rng(43210)
    # Random design points for descriptive scatter plots, not hypothesis tests.
    design=[]
    for cx,cy,z in rng.uniform([-160,-160,120],[160,160,270],(600,3)):
        design.append({'cx':cx,'cy':cy,'z':z,'coverage':objective(layout(40,z,cx,cy),order=16,orders=(5,3,5))})
    save('random_designs.csv',design)
    # Spatial maps preserve field analysis. A common colour scale is used.
    q=np.array([[x,y,z] for z in [137.5,150.,162.5] for y in np.linspace(-10,10,25) for x in np.linspace(-50,50,61)])
    maps=[]
    for name,s in layouts.items():
        p=pointwise(q,s,order=48);cov=1-np.prod(1-p,axis=1)
        for pos,v in zip(q,cov):maps.append({'layout':name,'x':pos[0],'y':pos[1],'z':pos[2],'coverage':v})
    save('coverage_maps.csv',maps)
    # Correct aggregation versus invalid aggregation of spatial means.
    agg=[];q,w=quad((11,7,15))
    for name,s in layouts.items():
        p=pointwise(q,s,order=48);means=w@p
        agg.append({'layout':name,'correct':float(w@(1-np.prod(1-p,axis=1))),'wrong_product_of_means':float(1-np.prod(1-means)),'union_bound_sum':float(means.sum())})
    save('aggregation_audit.csv',agg)
    # Ideal known horizontal bias scenario: not CFD/LSTM or real-time control.
    dyn=[];previous=np.zeros(2);movement=0.
    for k in range(13):
        bias=np.array([80*np.sin(k*np.pi/6),40*np.sin(k*np.pi/3)])
        nominal=opt.copy();nominal[:,:2]-=bias
        movement+=9*np.linalg.norm(-bias-previous);previous=-bias
        actual_mean=opt.copy();actual_mean[:,:2]+=bias
        dyn.append({'step':k,'bias_x':bias[0],'bias_y':bias[1],'static_coverage':objective(actual_mean,order=48,orders=(11,7,15)),
                    'ideal_compensated_coverage':objective(nominal+np.r_[bias,0],order=48,orders=(11,7,15)),
                    'cumulative_nominal_movement_m':movement})
    save('dynamic_bias_scenario.csv',dyn)
    save('optimized_coordinates.csv',[{'sensor':i+1,'x':s[0],'y':s[1],'z':s[2]} for i,s in enumerate(opt)])
    summary['best_Powell']=best;summary['relative_gain']=(summary['Optimized']['quadrature']/summary['Fixed']['quadrature']-1)
    summary['elapsed_seconds']=time.perf_counter()-start
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    (OUT/'environment.json').write_text(json.dumps({'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'cpu':platform.processor()},indent=2))
    print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
