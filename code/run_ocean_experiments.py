"""CMEMS-forced coverage replay with training-only deployment selection.

Input must be a real subset produced by download_public_cmems.py. The code never
manufactures replacement ocean records. Currents are reanalysis velocities;
the displacement rule and residual position law are declared assumptions.
"""
from pathlib import Path
from dataclasses import replace
import argparse, hashlib, json, time
import numpy as np
import pandas as pd
from scipy.cluster.vq import kmeans2, vq
from scipy.stats import ncx2
from numpy.polynomial.legendre import leggauss
from run_experiments import C, quad, layout, optimize, objective
from formula_reference import (truncation_normalizer, current_displacement,
    pointwise_union, sample_actual, spherical_hit, batch_interval,
    block_bootstrap_interval)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data'/'ocean'; INPUT=ROOT/'data'/'cmems'
SITES=('south_china_sea','kuroshio','north_atlantic')

def conditional_coverage(s,displacements,c=C,order=24,orders=(7,5,7),chunk=16):
    """F25: integrate each current-conditioned union, then return daily means.

    A displacement is shared by every sensor in a realization. Residual errors
    are conditionally independent, but sensors are not unconditionally
    independent after averaging over shared currents.
    """
    q,wq=quad(orders);g,wz=leggauss(order)
    lo=np.maximum(q[:,2]-c.R,c.zmin)
    hi=np.maximum(lo,np.minimum(q[:,2]+c.R,c.zmax))
    z=lo[:,None]+(g+1)*(hi-lo)[:,None]/2
    r2=np.maximum(c.R**2-(q[:,None,2]-z)**2,0)/c.sh**2
    nz=truncation_normalizer(s[:,2],c.sz,c.zmin,c.zmax)
    fz=np.exp(-.5*((z[:,:,None]-s[None,None,:,2])/c.sz)**2)/(c.sz*np.sqrt(2*np.pi)*nz)
    wzfz=fz*wz[None,:,None]*(hi-lo)[:,None,None]/2
    outputs=[]
    for start in range(0,len(displacements),chunk):
        b=np.asarray(displacements[start:start+chunk])
        centre=s[None,:,:2]+b[:,None,:]
        nxy=np.prod(truncation_normalizer(centre,c.sh,c.xmin,c.xmax),axis=2)
        if not np.all(nxy>0):raise FloatingPointError('Unresolved truncated Gaussian tail')
        nc=np.sum((q[None,:,None,:2]-centre[:,None,:,:])**2,axis=3)/c.sh**2
        disk=ncx2.cdf(r2[None,:,:,None],2,nc[:,:,None,:])
        p=np.sum(disk*wzfz[None,:,:,:],axis=2)/nxy[:,None,:]
        if not np.isfinite(p).all():raise FloatingPointError('Nonfinite sphere probability')
        if p.max()>1+1e-8 or p.min()<-1e-12:raise ArithmeticError('Invalid sphere probability')
        outputs.extend(pointwise_union(np.clip(p,0,1))@wq)
    return np.asarray(outputs)

def training_states(displacements,k=16,seed=6200):
    """F26: deterministic training-only representative states and frequencies."""
    centres,labels=kmeans2(np.asarray(displacements),k,iter=60,minit='++',
                           seed=np.random.default_rng(seed),missing='raise')
    # Reassign after the last centroid update; frequencies sum to one.
    labels,_=vq(displacements,centres)
    weights=np.bincount(labels,minlength=k)/len(displacements)
    return centres,weights

def load_currents():
    path=INPUT/'currents_150m.csv';manifest_path=INPUT/'provenance.json'
    if not path.exists() or not manifest_path.exists():
        raise FileNotFoundError('Real CMEMS subset is missing. Extract frozen_data.zip or run code/download_public_cmems.py; no account is required.')
    meta=json.loads(manifest_path.read_text(encoding='utf-8'))
    if hashlib.sha256(path.read_bytes()).hexdigest()!=meta['sha256_csv']:
        raise ValueError('CMEMS CSV hash does not match the extraction manifest')
    frame=pd.read_csv(path,parse_dates=['time'])
    if not frame.valid.all() or not np.isfinite(frame[['u_m_s','v_m_s']]).all().all():
        raise ValueError('Missing currents require an explicit exclusion protocol, not imputation')
    for site in SITES:
        sub=frame[frame.site==site]
        assert len(sub)==730 and sub.time.is_unique
        assert (sub.groupby('split').size()==365).all()
        assert set(sub.split)=={'train','test'}
    return frame,meta

def write_csv(name,rows):
    pd.DataFrame(rows).to_csv(OUT/name,index=False,float_format='%.17g')

def run(budget=300,seeds=20,horizon=600.,mc_batches=30,mc_trials=50000,resume=False):
    frame,meta=load_currents();OUT.mkdir(parents=True,exist_ok=True)
    if seeds!=20 or mc_batches!=30 or mc_trials!=50000:
        raise ValueError('Publication run requires 20 optimizer seeds and 30 x 50000 validation trials')
    cfg={'horizon_seconds':horizon,'budget':budget,'optimizer_seeds':list(range(5100,5120)),
         'training_year':2018,'held_out_year':2019,'sites':list(SITES),
         'representative_states':16,'clustering_seed':6200,
         'coarse_rule':[12,3,2,3],'refined_rule':[32,9,5,11],
         'validation_rule':[48,11,7,15],'mc_batches':mc_batches,'mc_trials':mc_trials,
         'block_bootstrap':{'block_days':14,'repetitions':2000},
         'source_csv_sha256':meta['sha256_csv'],
         'interpretation':'Observed-constrained ocean reanalysis forcing; assumed frozen-current displacement and residual Gaussian errors; not a sensor sea trial.'}
    previous_runs=pd.DataFrame();previous_traces=pd.DataFrame()
    if resume and (OUT/'optimizer_runs.csv').exists():
        assert json.loads((OUT/'parameters.json').read_text())==cfg,'Resume parameters differ'
        previous_runs=pd.read_csv(OUT/'optimizer_runs.csv');previous_traces=pd.read_csv(OUT/'optimizer_traces.csv')
    (OUT/'parameters.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    runs=[];traces=[];daily=[];states=[];selected=[];checks=[];mc=[];sens=[];summaries=[]
    synthetic=json.loads((ROOT/'data'/'summary.json').read_text())['best_Powell']
    synthetic_layout=layout(synthetic['pitch'],synthetic['z'],synthetic['cx'],synthetic['cy'])
    t0=time.perf_counter()
    for site_index,site in enumerate(SITES):
        sub=frame[frame.site==site]
        train=sub[sub.split=='train'];test=sub[sub.split=='test']
        bt=current_displacement(train[['u_m_s','v_m_s']].to_numpy(),horizon)
        bv=current_displacement(test[['u_m_s','v_m_s']].to_numpy(),horizon)
        centres,weights=training_states(bt)
        for j,(b,w) in enumerate(zip(centres,weights)):
            states.append({'site':site,'state':j,'bias_x_m':b[0],'bias_y_m':b[1],'weight':w})
        coarse=lambda s:float(weights@conditional_coverage(s,centres,order=12,orders=(3,2,3)))
        refined=lambda s:float(weights@conditional_coverage(s,centres,order=32,orders=(9,5,11)))
        site_runs=[]
        for method in ['Powell','PSO','AFSA','Random']:
            for seed in range(5100,5120):
                cached=previous_runs[(previous_runs.site==site)&(previous_runs.method==method)&(previous_runs.seed==seed)] if len(previous_runs) else pd.DataFrame()
                if len(cached):
                    assert len(cached)==1
                    row=cached.iloc[0].to_dict()
                    tr=previous_traces[(previous_traces.site==site)&(previous_traces.method==method)&(previous_traces.seed==seed)]
                    trace=tr.drop(columns='site').to_dict('records')
                else:
                    row,trace=optimize(method,seed,budget=budget,evaluator=coarse,refined_evaluator=refined)
                # The inherited field called 'test' is a finer TRAINING objective.
                if 'test' in row:row['training_refined']=row.pop('test')
                row['site']=site
                site_runs.append(row);runs.append(row)
                traces.extend(dict(site=site,**r) for r in trace)
                write_csv('optimizer_runs.csv',runs);write_csv('optimizer_traces.csv',traces)
                print(site,method,seed,'train=',round(row['training_refined'],7),flush=True)
        designs={'Fixed':layout(),'Synthetic-optimum':synthetic_layout}
        for method in ['Powell','PSO','AFSA','Random']:
            chosen=max((r for r in site_runs if r['method']==method),key=lambda r:r['training_refined'])
            designs[method]=layout(chosen['pitch'],chosen['z'],chosen['cx'],chosen['cy'])
            selected.append(chosen.copy())
        for name,s in designs.items():
            print(site,name,'full-year quadrature validation',flush=True)
            values=conditional_coverage(s,bv,order=48,orders=(11,7,15))
            for day,value in zip(test.time,values):daily.append({'site':site,'layout':name,'time':day,'coverage':value})
            train_full=conditional_coverage(s,bt,order=48,orders=(11,7,15)).mean()
            c32,w32=training_states(bt,32)
            check32=float(w32@conditional_coverage(s,c32,order=48,orders=(11,7,15)))
            check16=float(weights@conditional_coverage(s,centres,order=48,orders=(11,7,15)))
            checks.append({'site':site,'layout':name,'full_365_train':train_full,'states_16':check16,'states_32':check32})
            if name=='Powell':
                fine=conditional_coverage(s,bv,order=64,orders=(13,9,19))
                checks.append({'site':site,'layout':'Powell-test-refinement','test_48':values.mean(),'test_64':fine.mean(),'max_daily_change':np.max(abs(values-fine))})
            write_csv('daily_test_coverage.csv',daily)
        for name in ['Fixed','Synthetic-optimum','Powell']:
            print(site,name,'independent Monte Carlo batches',flush=True)
            s=designs[name]
            # Common random numbers pair layouts inside each independent batch.
            for rep in range(mc_batches):
                seed=81000+1000*site_index+rep;rng=np.random.default_rng(seed)
                n=mc_trials;indices=rng.integers(len(bv),size=n)
                target=rng.uniform([-50,-10,137.5],[50,10,162.5],(n,3))
                uniform=rng.random((n,9,3));means=np.broadcast_to(s,(n,9,3)).copy()
                means[:,:,:2]+=bv[indices,None,:]
                actual=sample_actual(uniform,means)
                hit=spherical_hit(actual,target,20.)
                raw=OUT/'mc_raw';raw.mkdir(exist_ok=True)
                np.savez_compressed(raw/f'{site}_{name.replace("-","_")}_{rep:02d}.npz',hit=np.packbits(hit))
                mc.append({'site':site,'layout':name,'rep':rep,'seed':seed,'trials':n,'hits':int(hit.sum()),'coverage':hit.mean()})
        site_daily=pd.DataFrame(daily);site_daily=site_daily[site_daily.site==site]
        fixed=site_daily[site_daily.layout=='Fixed'].coverage.to_numpy()
        for name in designs:
            x=site_daily[site_daily.layout==name].coverage.to_numpy();ci=block_bootstrap_interval(x)
            gain_ci=block_bootstrap_interval(x-fixed)
            r={'site':site,'layout':name,'test_coverage':x.mean(),'temporal_ci_low':ci[0],'temporal_ci_high':ci[1],
               'gain_over_fixed':(x-fixed).mean(),'gain_ci_low':gain_ci[0],'gain_ci_high':gain_ci[1]}
            vals=[r['coverage'] for r in mc if r['site']==site and r['layout']==name]
            if vals:r.update(zip(['mc_mean','mc_ci_low','mc_ci_high'],batch_interval(vals)))
            summaries.append(r)
        # Fixed selected layouts; no test-based tuning or reoptimization here.
        for tau in [0.,60.,300.,600.,1200.]:
            for name in ['Fixed','Synthetic-optimum','Powell']:
                d=current_displacement(test[['u_m_s','v_m_s']].to_numpy(),tau)
                value=conditional_coverage(designs[name],d,order=32,orders=(9,5,11)).mean()
                sens.append({'site':site,'layout':name,'horizon_s':tau,'sigma_h_m':120.,'coverage':value})
        for sh in [30.,60.,180.]:
            for name in ['Fixed','Synthetic-optimum','Powell']:
                value=conditional_coverage(designs[name],bv,c=replace(C,sh=sh),order=32,orders=(9,5,11)).mean()
                sens.append({'site':site,'layout':name,'horizon_s':horizon,'sigma_h_m':sh,'coverage':value})
        for filename,rows in [('daily_test_coverage.csv',daily),('training_states.csv',states),('selected_designs.csv',selected),
            ('accuracy_checks.csv',checks),('mc_batches.csv',mc),('sensitivity.csv',sens),('summary.csv',summaries)]:write_csv(filename,rows)
    (OUT/'execution.json').write_text(json.dumps({'elapsed_seconds':time.perf_counter()-t0,'source_hash':meta['sha256_csv']},indent=2))
    print('Completed three real-current replay sites',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--budget',type=int,default=300)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args();run(budget=args.budget,resume=args.resume)
