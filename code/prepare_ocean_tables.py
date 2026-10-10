"""Data-driven LaTeX tables and macros for the real-current replay."""
from pathlib import Path
import json,numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data';O=D/'ocean';G=ROOT/'results'/'tables'
NAMES={'south_china_sea':'South China Sea','kuroshio':'Kuroshio vicinity','north_atlantic':'North Atlantic'}
def table(name,head,rows,fmt):
    text='\\begin{tabular}{'+fmt+'}\n\\toprule\n'+head+r' \\'+ '\n\\midrule\n'+'\n'.join(' & '.join(map(str,r))+r' \\' for r in rows)+'\n\\bottomrule\n\\end{tabular}\n'
    (G/name).write_text(text,encoding='utf-8')
def main():
    G.mkdir(parents=True,exist_ok=True);s=pd.read_csv(O/'summary.csv');r=pd.read_csv(O/'optimizer_runs.csv')
    meta=json.loads((D/'cmems'/'provenance.json').read_text());rows=[];macros=[]
    for i,site in enumerate(meta['sites']):
        a=s[s.site==site['site']];fixed=a[a.layout=='Fixed'].iloc[0];synthetic=a[a.layout=='Synthetic-optimum'].iloc[0];opt=a[a.layout=='Powell'].iloc[0]
        rows.append([NAMES[site['site']],f'{fixed.test_coverage*100:.4f}',f'{synthetic.test_coverage*100:.4f}',
            f'{opt.test_coverage*100:.4f}',f'{opt.gain_over_fixed*100:.4f}',f'[{opt.gain_ci_low*100:.4f}, {opt.gain_ci_high*100:.4f}]'])
        for key,value in [('Fixed',fixed.test_coverage*100),('Synthetic',synthetic.test_coverage*100),('Opt',opt.test_coverage*100),
                          ('Gain',opt.gain_over_fixed*100),('GainLow',opt.gain_ci_low*100),('GainHigh',opt.gain_ci_high*100)]:
            macros.append('\\newcommand{\\Ocean'+['SCS','Kuro','NA'][i]+key+'}{'+f'{value:.4f}'+'}')
    table('ocean_coverage.tex',r'Site & Fixed (\%) & Synthetic optimum (\%) & Current-aware Powell (\%) & Gain (pp) & Temporal gain 95\% CI (pp)',rows,'lrrrrl')
    rows=[]
    for site in meta['sites']:
        rows.append([NAMES[site['site']],f"{site['actual_longitude']:.1f}",f"{site['actual_latitude']:.1f}",
                     '130.666, 155.851','365 / 365'])
    table('ocean_sites.tex',r'Site & Longitude ($^\circ$E) & Latitude ($^\circ$N) & Retained depths (m) & Train / test days',rows,'lrrll')
    rows=[]
    for _,x in s[s.layout.isin(['Fixed','Synthetic-optimum','Powell'])].iterrows():
        rows.append([NAMES[x.site],x.layout,f'{x.test_coverage*100:.4f}',f'{x.mc_mean*100:.4f}',f'[{x.mc_ci_low*100:.4f}, {x.mc_ci_high*100:.4f}]'])
    table('ocean_mc.tex',r'Site & Layout & Integral (\%) & MC mean (\%) & MC 95\% CI (\%)',rows,'llrrl')
    rows=[]
    for (site,method),x in r.groupby(['site','method'],sort=False):
        rows.append([NAMES[site],method,f'{x.training_refined.mean()*100:.4f}',f'{x.training_refined.std()*100:.4f}',
                     f'{x.training_refined.min()*100:.4f}',f'{x.evaluations.mean():.1f}',f'{x.seconds.median():.3f}'])
    table('ocean_algorithms.tex',r'Site & Method & Mean (\%) & SD (pp) & Minimum (\%) & Mean calls & Median seconds',rows,'llrrrrr')
    selected=pd.read_csv(O/'selected_designs.csv');rows=[]
    for _,x in selected.iterrows():rows.append([NAMES[x.site],x.method,int(x.seed),*[f'{x[v]:.3f}' for v in ['cx','cy','pitch','z']]])
    table('ocean_designs.tex',r'Site & Method & Seed & $c_x$ (m) & $c_y$ (m) & $d$ (m) & $z_0$ (m)',rows,'llrrrrr')
    checks=pd.read_csv(O/'accuracy_checks.csv');rows=[]
    for _,x in checks[checks.layout.isin(['Fixed','Synthetic-optimum','Powell'])].iterrows():
        rows.append([NAMES[x.site],x.layout,f'{x.full_365_train*100:.5f}',f'{x.states_16*100:.5f}',f'{x.states_32*100:.5f}'])
    table('ocean_states.tex',r'Site & Layout & Full training year (\%) & 16 states (\%) & 32 states (\%)',rows,'llrrr')
    transfer=pd.read_csv(O/'transfer_gains.csv');rows=[]
    for _,x in transfer.iterrows():
        rows.append([NAMES[x.site],f'{x.mean_gain*100:.4f}',f'[{x.low*100:.4f}, {x.high*100:.4f}]'])
    table('ocean_transfer.tex',r'Site & Powell minus synthetic optimum (pp) & Paired temporal 95\% CI (pp)',rows,'lrl')
    (G/'ocean_numbers.tex').write_text('\n'.join(macros),encoding='utf-8')
    print('Generated 7 real-current tables and numeric macros')
if __name__=='__main__':main()
