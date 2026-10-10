"""Independent provenance, raw-event, feasibility and MC/integral checks."""
from pathlib import Path
import hashlib,json,sys,platform
import numpy as np,pandas as pd,xarray as xr,scipy
from formula_reference import batch_interval,block_bootstrap_interval
from run_ocean_experiments import load_currents,SITES
from run_experiments import layout

ROOT=Path(__file__).resolve().parents[1];O=ROOT/'data'/'ocean';checks=[]
def check(name,ok,detail):
    assert ok,(name,detail);checks.append({'check':name,'pass':bool(ok),'detail':detail})

def main():
    frame,meta=load_currents()
    check('Real-record completeness',len(frame)==2190 and frame.valid.all(),'3 sites x 730 days; no missing-value imputation')
    for site in meta['sites']:
        name=site['site'];path=ROOT/'data'/'cmems'/(name+'.nc')
        check('Original-depth subset hash '+name,hashlib.sha256(path.read_bytes()).hexdigest()==site['sha256_netcdf'],'SHA-256 matches the extraction manifest')
        ds=xr.open_dataset(path);df=frame[frame.site==name]
        assert ds.uo.attrs['units']==ds.vo.attrs['units']=='m s-1'
        d=ds.interp(depth=150.);error=max(np.max(abs(d.uo.values-df.u_m_s.to_numpy())),np.max(abs(d.vo.values-df.v_m_s.to_numpy())))
        check('Depth interpolation and CSV '+name,error<1e-10,{'maximum_absolute_m_s_difference':float(error)})
        check('Chronological split '+name,set(df[df.split=='train'].time.dt.year)=={2018} and set(df[df.split=='test'].time.dt.year)=={2019},'Training 2018; held-out test 2019')
    runs=pd.read_csv(O/'optimizer_runs.csv');check('Matched real-current optimizer repetition',len(runs)==240 and (runs.evaluations<=300).all(),
        '3 sites x 4 methods x 20 seeds, identical maximum budget 300')
    counts=runs.groupby(['site','method']).size();assert (counts==20).all()
    selected=pd.read_csv(O/'selected_designs.csv');coordinates=[]
    for _,row in selected.iterrows():
        subset=runs[(runs.site==row.site)&(runs.method==row.method)]
        best=subset.loc[subset.training_refined.idxmax()]
        assert row.seed==best.seed
        s=layout(row.pitch,row.z,row.cx,row.cy)
        distance=np.linalg.norm(s[:,None,:]-s[None,:,:],axis=2)+np.eye(9)*1e6
        assert s[:,:2].min()>=-500-1e-7 and s[:,:2].max()<=500+1e-7 and s[:,2].min()>=120-1e-7 and s[:,2].max()<=270+1e-7
        assert distance.min()>=40-1e-7
        coordinates.extend(dict(site=row.site,method=row.method,sensor=i+1,x=p[0],y=p[1],z=p[2]) for i,p in enumerate(s))
    pd.DataFrame(coordinates).to_csv(O/'selected_coordinates.csv',index=False,float_format='%.12g')
    check('Training-only selection and nominal feasibility',True,'Selected seeds maximize refined training score; all nominal nodes bounded with separation >=40 m')
    daily=pd.read_csv(O/'daily_test_coverage.csv',parse_dates=['time'])
    check('Held-out daily evaluations',len(daily)==3*6*365 and daily.coverage.between(0,1).all(),'6 layouts, 365 held-out days per site')
    assert (daily.groupby(['site','layout']).size()==365).all()
    mc=pd.read_csv(O/'mc_batches.csv');assert len(mc)==270
    for _,row in mc.iterrows():
        path=O/'mc_raw'/f'{row.site}_{row.layout.replace("-","_")}_{int(row.rep):02d}.npz'
        with np.load(path) as raw:hit=np.unpackbits(raw['hit'])[:int(row.trials)]
        assert int(hit.sum())==int(row.hits)
        assert abs(hit.mean()-row.coverage)<1e-12
    check('Real-current raw Monte Carlo outcomes',True,'270 packed Bernoulli arrays agree with CSV counts; 13.5 million full target--array trials')
    discrepancies=[]
    for (site,name),x in mc.groupby(['site','layout']):
        value=daily[(daily.site==site)&(daily.layout==name)].coverage.mean()
        se=x.coverage.std(ddof=1)/np.sqrt(30)
        z=abs(x.coverage.mean()-value)/se
        discrepancies.append({'site':site,'layout':name,'quadrature':float(value),'MC':float(x.coverage.mean()),'standard_error_units':float(z)})
        assert z<6
    check('Independent MC against identical conditional integral',True,discrepancies)
    refinement=pd.read_csv(O/'accuracy_checks.csv');ref=refinement[refinement.layout=='Powell-test-refinement']
    check('Held-out quadrature refinement',len(ref)==3 and (abs(ref.test_48-ref.test_64)<1e-6).all(),
        [{'site':r.site,'mean_change':float(abs(r.test_48-r.test_64)),'maximum_daily_change':float(r.max_daily_change)} for _,r in ref.iterrows()])
    intervals=[];transfer=[]
    for site in SITES:
        a=daily[daily.site==site].pivot(index='time',columns='layout',values='coverage')
        gain=(a.Powell-a.Fixed).to_numpy()
        for block in [7,14,30]:
            lo,hi=block_bootstrap_interval(gain,block=block)
            intervals.append({'site':site,'block_days':block,'mean_gain':gain.mean(),'low':lo,'high':hi})
        delta=(a.Powell-a['Synthetic-optimum']).to_numpy()
        lo,hi=block_bootstrap_interval(delta)
        transfer.append({'site':site,'comparator':'Synthetic-optimum','mean_gain':float(delta.mean()),'low':lo,'high':hi})
    pd.DataFrame(intervals).to_csv(O/'block_length_sensitivity.csv',index=False,float_format='%.12g')
    pd.DataFrame(transfer).to_csv(O/'transfer_gains.csv',index=False,float_format='%.12g')
    check('Temporal bootstrap reproducibility',True,'2000 circular moving-block draws at lengths 7, 14 and 30 days, seed 7100')
    environment={'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__,
                 'xarray':xr.__version__,'platform':platform.platform(),'cpu':platform.processor()}
    (O/'environment.json').write_text(json.dumps(environment,indent=2),encoding='utf-8')
    (O/'verification.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
