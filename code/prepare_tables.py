from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.special import ndtr
from run_experiments import layout,pointwise,quad,C
R=Path(__file__).resolve().parents[1];D=R/'data';P=R/'results';G=P/'tables';G.mkdir(parents=True,exist_ok=True)
s=json.loads((D/'summary.json').read_text());runs=pd.read_csv(D/'optimizer_runs.csv');mc=pd.read_csv(D/'mc_batches.csv')
def table(name,head,rows,fmt):
    text='\\begin{tabular}{'+fmt+'}\n\\toprule\n'+head+r' \\'+ '\n\\midrule\n'+'\n'.join(' & '.join(map(str,r))+r' \\' for r in rows)+'\n\\bottomrule\n\\end{tabular}\n'
    (G/name).write_text(text,encoding='utf-8')
rows=[]
for name in ['Single','Fixed','Pitch-only','Depth-only','Optimized']:
    a=s[name];rows.append([name,'1' if name=='Single' else '9',f'{a["quadrature"]*100:.4f}',f'{a["mc_mean"]*100:.4f}',f'{a["mc_sd"]*100:.4f}',f'[{a["ci_low"]*100:.4f}, {a["ci_high"]*100:.4f}]'])
table('coverage.tex',r'Layout & $N$ & Quadrature (\%) & MC mean (\%) & SD (pp) & 95\% CI (\%)',rows,'lrrrrl')
rows=[]
for method in ['Powell','PSO','AFSA','Random']:
    a=runs[runs.method==method];rows.append([method,f'{a.test.mean()*100:.5f}',f'{a.test.std()*100:.5f}',f'{a.test.min()*100:.5f}',f'{a.evaluations.mean():.1f}',f'{a.seconds.median():.3f}'])
table('algorithms.tex',r'Method & Mean (\%) & SD (pp) & Minimum (\%) & Mean calls & Median time (s)',rows,'lrrrrr')
cv=pd.read_csv(D/'quadrature_convergence.csv');rows=[]
for name in ['Single','Fixed','Optimized']:
    a=cv[cv.layout==name];rows.append([name]+[f'{x:.9f}' for x in a.coverage])
table('quadrature.tex','Layout & Coarse & Intermediate & Fine & Finest',rows,'lrrrr')
opt=pd.read_csv(D/'optimized_coordinates.csv')[['x','y','z']].to_numpy();q,w=quad((15,9,25));p=pointwise(q,opt,order=80);rows=[]
for i,(ss,prob) in enumerate(zip(opt,(w@p))):rows.append([i+1,*[f'{v:.3f}' for v in ss],f'{prob*100:.5f}'])
table('coordinates.tex',r'Node & $x$ (m) & $y$ (m) & $z$ (m) & Single-node spatial mean (\%)',rows,'rrrrr')
rows=[]
for _,r in pd.read_csv(D/'aggregation_audit.csv').iterrows():rows.append([r.layout,f'{r.correct:.8f}',f'{r.wrong_product_of_means:.8f}',f'{r.union_bound_sum:.8f}'])
table('aggregation.tex','Layout & Correct union & Product of spatial means & Sum of spatial means',rows,'lrrr')
rows=[]
for sh in [30,60,120,180]:
    a=pd.read_csv(D/'uncertainty_sensitivity.csv');a=a[(a.sigma_h==sh)&(a.layout=='Fixed')];rows.append([sh]+[f'{v*100:.4f}' for v in a.coverage])
table('horizontal.tex',r'$\sigma_h$ (m) & $\sigma_z=0$ & $\sigma_z=15$ & $\sigma_z=40$ & $\sigma_z=60$',rows,'rrrrr')
gain=(s['Optimized']['quadrature']-s['Fixed']['quadrature'])*100
blind=gain/(1-s['Fixed']['quadrature'])
macros={'FixedCov':f'{s["Fixed"]["quadrature"]*100:.4f}','OptCov':f'{s["Optimized"]["quadrature"]*100:.4f}','AbsGain':f'{gain:.4f}','RelGain':f'{s["relative_gain"]*100:.2f}','BlindGain':f'{blind:.2f}','MCOpt':f'{s["Optimized"]["mc_mean"]*100:.4f}','MCLow':f'{s["Optimized"]["ci_low"]*100:.4f}','MCHigh':f'{s["Optimized"]["ci_high"]*100:.4f}','PairGain':f'{s["paired_gain"]["mean"]*100:.4f}','PairLow':f'{s["paired_gain"]["low"]*100:.4f}','PairHigh':f'{s["paired_gain"]["high"]*100:.4f}'}
(G/'numbers.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+v+'}' for k,v in macros.items()))
print('Numerical tables and macros generated; no bibliography is produced')
