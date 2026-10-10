from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.special import ndtr
from run_experiments import layout,pointwise,quad,C
R=Path(__file__).resolve().parents[1];D=R/'data';P=R/'paper';G=P/'generated';G.mkdir(exist_ok=True)
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
# Primary DOI registry records; replace the erroneous source citations, preserve
# valid source references and their research topics. Full audit is separate.
records=json.loads((D/'reference_records.json').read_text());missing={'10.1145/972627.972631':('Sensor deployment and target localization in distributed sensor networks','Zou, Yi and Chakrabarty, Krishnendu','ACM Transactions on Embedded Computing Systems','3','61--91','2004'),'10.1146/annurev-fluid-010719-060214':('Machine Learning for Fluid Mechanics','Brunton, Steven L. and Noack, Bernd R. and Koumoutsakos, Petros','Annual Review of Fluid Mechanics','52','477--508','2020'),'10.1016/j.coastaleng.2024.104532':('Storm surge modeling in the AI era: Using LSTM-based machine learning for enhancing forecasting accuracy','Giaremis, S. and Nader, N. and Dawson, C. and Kaiser, C. and Nikidis, E. and Kaiser, H.','Coastal Engineering','191','104532','2024')}
keys={'10.1145/972627.972631':'zou2004','10.1109/JOE.2019.2935481':'zhao2020','10.1080/01621459.1949.10483310':'metropolis1949','10.1016/j.engappai.2023.107062':'li2023','10.1016/j.joes.2021.03.001':'murray2022','10.1016/j.oceaneng.2023.115486':'qiu2023','10.1038/s42254-023-00622-y':'vinuesa2023','10.1146/annurev-fluid-010719-060214':'brunton2020','10.1016/j.coastaleng.2024.104532':'giaremis2024','10.1016/j.swevo.2024.101486':'wang2024','10.1109/JIOT.2023.3324392':'zhou2024','10.1126/scirobotics.ade7811':'masmitja2023','10.1038/s42256-020-00272-0':'novati2021','10.1016/j.oceaneng.2025.121980':'dong2025','10.1038/s41467-024-50926-5':'lu2024','10.1016/j.oceaneng.2022.113479':'mauro2023','10.1016/j.adhoc.2005.01.004':'akyildiz2005','10.1145/941350.941367':'huang2003','10.1016/j.jpdc.2004.03.019':'zouuncertainty2004','10.1093/comjnl/7.2.155':'powell1964','10.1109/ICNN.1995.488968':'kennedy1995','10.1109/4235.985692':'clerc2002'}
keys['10.1109/JOE.2019.2916119']='liu2019'
bib=[]
for r in records:
    doi=r['doi']
    if doi not in keys:continue
    if 'error' in r:
        if doi not in missing:continue
        title,author,journal,volume,pages,year=missing[doi]
    else:
        title=r['title'][0];author=' and '.join(x.get('family','')+', '+x.get('given','') for x in r['authors']);journal=r['journal'][0];volume=r.get('volume') or '';pages=r.get('pages') or r.get('article_number') or '';year=str(r['year'].get('date-parts',[[keys[doi][-4:]]])[0][0])
    bib.append('@article{'+keys[doi]+',\n  '+',\n  '.join(k+'={'+str(v)+'}' for k,v in [('title','{'+title+'}'),('author',author),('journal',journal),('volume',volume),('pages',pages),('year',year),('doi',doi)])+'\n}')
bib.append('@article{li2002,author={Li, Xiao Lei and Shao, Zhi Jiang and Qian, Ji Xin},title={An Optimizing Method Based on Autonomous Animats: Fish-swarm Algorithm},journal={Systems Engineering---Theory and Practice},volume={22},number={11},pages={32--38},year={2002},doi={10.12011/1000-6788(2002)11-32}}')
bib.append('@book{papoulis2002,author={Papoulis, Athanasios and Pillai, S. Unnikrishna},title={Probability, Random Variables, and Stochastic Processes},edition={4},publisher={McGraw-Hill},year={2002}}')
# Preserve the separately verified marine-data citations on every regeneration.
bib.append((P/'ocean_references.bib').read_text(encoding='utf-8').strip())
(P/'references.bib').write_text('\n\n'.join(bib),encoding='utf-8')
print('Tables, numerical macros and bibliography generated')
