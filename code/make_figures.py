"""Publication figures; every numerical panel reads the archived experiment data."""
from pathlib import Path
import sys,json,textwrap
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from scipy.special import ndtr,ndtri
from run_experiments import C,layout,pointwise,quad
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data';F=ROOT/'paper'/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.labelsize':9,'axes.titlesize':10,'legend.fontsize':8,'pdf.fonttype':42,'savefig.dpi':320})
colors=['#2166ac','#d95f02','#1b9e77','#7570b3']
def finish(fig,name):
    fig.savefig(F/(name+'.pdf'),bbox_inches='tight');fig.savefig(F/(name+'.png'),bbox_inches='tight');plt.close(fig)
def label(ax,s):ax.set_title(s,loc='left',fontweight='bold')
def main():
    summary=json.loads((D/'summary.json').read_text());maps=pd.read_csv(D/'coverage_maps.csv');sweep=pd.read_csv(D/'parameter_sweeps.csv');sens=pd.read_csv(D/'uncertainty_sensitivity.csv');design=pd.read_csv(D/'random_designs.csv')
    fig,ax=plt.subplots(figsize=(10,6.1));ax.set(xlim=(0,10),ylim=(0,6));ax.axis('off')
    boxes=[(.1,4.15,3.0,1.6,'1  Physical model','Target volume 100 × 20 × 25 m\nActual positions in bounded water box\nGaussian horizontal / truncated depth\nPositive depth points downward'),(3.5,4.15,3.0,1.6,'2  Pointwise probability',r'$p_i(q)=P(\|S_i-q\|\leq R)$'+'\n'+r'$P_{\rm cov}(q)=1-\prod_i(1-p_i(q))$'+'\n'+r'$C=|\Omega|^{-1}\int_\Omega P_{\rm cov}(q)\,dq$'),(6.9,4.15,3.,1.6,'3  Feasible deployment','9-node planar grid\nCentre, pitch and common depth\nNominal separation ≥ 40 m\nBounds enforced by parameterization'),(.1,1.95,3.,1.65,'4  Numerical optimization','Bounded Powell direction search\nPSO / specified AFSA / random search\nSame objective and 600-call allowance\n20 independent initializations'),(3.5,1.95,3.,1.65,'5  Independent validation','Refined deterministic quadrature\n30 × 50,000 IID target-array trials\nPaired comparison of layouts\nSampling CI separate from grid error'),(6.9,1.95,3.,1.65,'6  Analysis and outputs','Corrected maps, coordinates and gains\nHorizontal / vertical sensitivity\nConvergence and method comparison\nAll data, code and random seeds')]
    for i,(x,y,w,h,title,body) in enumerate(boxes):
        if i==3:x=6.9
        if i==5:x=.1
        body=body.replace('Actual positions in bounded water box','Actual positions in physical bounds').replace('Gaussian horizontal / truncated depth','Gaussian errors; depth truncation').replace('Bounds enforced by parameterization','Bounds built into grid variables').replace('Same objective and 600-call allowance','Same objective; 600-call allowance').replace('Refined deterministic quadrature','Refined quadrature').replace('Sampling CI separate from grid error','Sampling CI versus numerical error').replace('Corrected maps, coordinates and gains','Maps, coordinates and coverage gain').replace('Convergence and method comparison','Convergence and algorithm comparison')
        co=colors[i%4];ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.02',facecolor=co+'10',edgecolor=co,lw=1.1));ax.text(x+.12,y+h-.23,title,weight='bold',color=co,fontsize=9.3);ax.text(x+.12,y+h-.53,body,va='top',fontsize=8,linespacing=1.6)
    for a,b in [((3.12,4.9),(3.46,4.9)),((6.52,4.9),(6.86,4.9)),((8.4,4.1),(8.4,3.65)),((6.86,2.8),(6.52,2.8)),((3.46,2.8),(3.12,2.8))]:ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=18,color='#4d4d4d'))
    ax.text(.15,1.25,'RQ1  Consistent probability model    RQ2  Sensitivity to 3D errors    RQ3  Optimization benefit and cost',fontsize=9,weight='bold')
    ax.text(.15,.65,'Extension interface: time-varying error parameters → re-optimization → movement / energy constraints',fontsize=10,color='#5e3c99')
    ax.text(.15,.2,'CFD–LSTM and closed-loop real-time control require independent development and validation.',fontsize=9,color='#555555')
    finish(fig,'fig1_framework')
    fig=plt.figure(figsize=(10,7));axs=[fig.add_subplot(221),fig.add_subplot(222),fig.add_subplot(223,projection='3d'),fig.add_subplot(224,projection='3d')]
    rng=np.random.default_rng(7654);x=np.linspace(-500,500,600);a=ndtr(-500/120);b=ndtr(500/120);sample=120*ndtri(a+(b-a)*rng.random(10000))
    axs[0].hist(sample,bins=40,density=True,color=colors[0],alpha=.45);axs[0].plot(x,np.exp(-.5*(x/120)**2)/(120*np.sqrt(2*np.pi)*(b-a)),c=colors[0]);axs[0].set(xlabel='Horizontal error (m)',ylabel='Density (1/m)');label(axs[0],'a  Horizontal positioning uncertainty')
    z=np.linspace(120,270,500)
    for mu,col in [(150,colors[0]),(120,colors[1])]:
        pdf=np.exp(-.5*((z-mu)/40)**2)/(40*np.sqrt(2*np.pi)*(ndtr((270-mu)/40)-ndtr((120-mu)/40)));axs[1].plot(z,pdf,label=f'Nominal depth {mu} m',c=col)
    axs[1].set(xlabel='Actual absolute depth (m)',ylabel='Density (1/m)');axs[1].legend();label(axs[1],'b  Physical depth truncation')
    mm=maps[(maps.layout=='Optimized')&(maps.z==150)];xx=mm.x.unique();yy=mm.y.unique();X,Y=np.meshgrid(xx,yy);Z=mm.coverage.to_numpy().reshape(len(yy),len(xx))*100
    axs[2].plot_surface(X,Y,Z,cmap='viridis',linewidth=0);axs[2].set(xlabel='x (m)',ylabel='y (m)',zlabel='Coverage (%)');label(axs[2],'c  Coverage field at depth 150 m')
    sc=axs[3].scatter(design.cx,design.cy,design.z,c=design.coverage*100,cmap='viridis',s=5,alpha=.8);axs[3].set(xlabel='Centre x (m)',ylabel='Centre y (m)',zlabel='Nominal depth (m)');axs[3].invert_zaxis();label(axs[3],'d  600 sampled feasible designs');fig.colorbar(sc,ax=axs[3],shrink=.55,pad=.11,label='Coverage (%)');fig.tight_layout();finish(fig,'fig2_distributions_field')
    fig,axs=plt.subplots(2,2,figsize=(10,6.8))
    for ax,name in zip(axs[0],['Fixed','Optimized']):
        tab=sens[sens.layout==name].pivot(index='sigma_z',columns='sigma_h',values='coverage');im=ax.imshow(tab*100,origin='lower',aspect='auto',cmap='viridis');ax.set_xticks(range(4),tab.columns.astype(int));ax.set_yticks(range(4),tab.index.astype(int));ax.set(xlabel='Horizontal SD (m)',ylabel='Vertical SD (m)');label(ax,f'{"a" if name=="Fixed" else "b"}  {name} layout');
        for i in range(4):
            for j in range(4):ax.text(j,i,f'{tab.iloc[i,j]*100:.2f}',ha='center',va='center',color='white' if tab.iloc[i,j]<tab.to_numpy().max()*.6 else 'black',fontsize=8)
        fig.colorbar(im,ax=ax,label='Coverage (%)')
    for ax,var,title in [(axs[1,0],'radius','c  Sensing-radius sensitivity'),(axs[1,1],'depth','d  Nominal-depth sensitivity')]:
        for name,col in [('Fixed',colors[0]),('Optimized',colors[1])]:
            df=sweep[(sweep.variable==var)&(sweep.layout==name)];ax.plot(df.value,df.coverage*100,label=name,c=col)
        ax.set(xlabel=('Radius (m)' if var=='radius' else 'Nominal depth (m)'),ylabel='Coverage (%)');label(ax,title);ax.legend()
    fig.tight_layout();finish(fig,'fig3_sensitivity')
    fig,axs=plt.subplots(2,3,figsize=(10,6))
    for ax,var,letter in zip(axs[0],['x','y','pitch'],'abc'):
        for name,col in [('Fixed',colors[0]),('Optimized',colors[1])]:
            df=sweep[(sweep.variable==var)&(sweep.layout==name)];ax.plot(df.value,df.coverage*100,c=col,label=name)
        ax.set(xlabel={'x':'Centre x offset (m)','y':'Centre y offset (m)','pitch':'Grid pitch (m)'}[var],ylabel='Coverage (%)');label(ax,letter+'  '+{'x':'X sensitivity','y':'Y sensitivity','pitch':'Nominal separation'}[var]);ax.legend()
    for ax,name,letter in zip(axs[1],['Fixed','Pitch-only','Optimized'],'def'):
        mm=maps[(maps.layout==name)&(maps.z==150)];arr=mm.coverage.to_numpy().reshape(25,61)*100;im=ax.imshow(arr,origin='lower',extent=[-50,50,-10,10],vmin=2,vmax=4.6,cmap='viridis',aspect='auto');ax.set(xlabel='Target x (m)',ylabel='Target y (m)');label(ax,letter+'  '+name)
    fig.subplots_adjust(wspace=.38,hspace=.42,right=.88)
    cax=fig.add_axes([.91,.13,.018,.30]);fig.colorbar(im,cax=cax,label='Pointwise coverage (%)');finish(fig,'fig4_constraints_maps')
    mc=pd.read_csv(D/'mc_batches.csv');cv=pd.read_csv(D/'quadrature_convergence.csv')
    fig,axs=plt.subplots(1,2,figsize=(10,3.6));names=['Single','Fixed','Pitch-only','Depth-only','Optimized']
    xs=np.arange(5);ms=np.array([summary[n]['mc_mean'] for n in names])*100;lo=np.array([summary[n]['ci_low'] for n in names])*100;hi=np.array([summary[n]['ci_high'] for n in names])*100
    axs[0].errorbar(xs,ms,yerr=[ms-lo,hi-ms],fmt='o',capsize=4,label='MC mean and 95% CI',color=colors[0]);axs[0].scatter(xs,[summary[n]['quadrature']*100 for n in names],marker='x',s=55,c=colors[1],label='Refined quadrature');axs[0].set_xticks(xs,names,rotation=20);axs[0].set(ylabel='Expected covered volume (%)');axs[0].legend();label(axs[0],'a  Independent same-target validation')
    for name,col in zip(['Fixed','Pitch-only','Optimized'],colors):
        df=cv[cv.layout==name];ref=df.coverage.iloc[-1];axs[1].semilogy(np.arange(1,4),abs(df.coverage.to_numpy()[:3]-ref)+1e-16,'o-',label=name,c=col)
    axs[1].set_xticks([1,2,3],['Coarse','Intermediate','Fine']);axs[1].set(ylabel='Absolute change from finest rule');axs[1].legend();label(axs[1],'b  Numerical refinement');fig.tight_layout();finish(fig,'fig5_validation')
    runs=pd.read_csv(D/'optimizer_runs.csv');trace=pd.read_csv(D/'optimizer_trace.csv');methods=['Powell','PSO','AFSA','Random'];fig,axs=plt.subplots(1,3,figsize=(10,3.6))
    for method,col in zip(methods,colors):
        paths=[]
        for seed in range(2100,2120):
            a=trace[(trace.method==method)&(trace.seed==seed)].best.to_numpy();paths.append(np.pad(a,(0,600-len(a)),mode='edge'))
        paths=np.array(paths)*100;median=np.median(paths,axis=0);axs[0].plot(np.arange(1,601),median,label=method,c=col);axs[0].fill_between(np.arange(1,601),*np.quantile(paths,[.25,.75],axis=0),alpha=.12,color=col)
    axs[0].set(xlabel='Objective evaluations (maximum 600)',ylabel='Best training coverage (%)');axs[0].legend();label(axs[0],'a  Median and interquartile range')
    for ax,col,title in [(axs[1],'test','b  Refined final coverage'),(axs[2],'seconds','c  Measured optimizer runtime')]:
        data=[runs[runs.method==m][col].to_numpy()*(100 if col=='test' else 1) for m in methods];ax.boxplot(data,tick_labels=methods,showfliers=False)
        for i,arr in enumerate(data):ax.scatter(i+1+np.linspace(-.08,.08,len(arr)),arr,s=8,color=colors[i],alpha=.7)
        ax.tick_params(axis='x',rotation=25);ax.set(ylabel='Coverage (%)' if col=='test' else 'Elapsed seconds');label(ax,title)
    fig.tight_layout();finish(fig,'fig6_algorithms')
    # Supplementary panels preserve the original boxplot analysis with correct units.
    fig,axs=plt.subplots(1,3,figsize=(10,3.3))
    for ax,col,letter in zip(axs,['cx','cy','z'],'abc'):
        bins=np.linspace(design[col].min(),design[col].max()+1e-5,5);groups=[design[(design[col]>=a)&(design[col]<b)].coverage*100 for a,b in zip(bins[:-1],bins[1:])];ax.boxplot(groups,tick_labels=[f'{a:.0f}–{b:.0f}' for a,b in zip(bins[:-1],bins[1:])]);ax.tick_params(axis='x',rotation=20);ax.set(xlabel=f'{col} interval (m)',ylabel='Coverage (%)');label(ax,letter+'  Descriptive sampled-design spread')
    fig.tight_layout();finish(fig,'figS1_design_boxplots')
    fig,axs=plt.subplots(1,2,figsize=(10,3.8));dyn=pd.read_csv(D/'dynamic_bias_scenario.csv');axs[0].plot(dyn.step,dyn.static_coverage*100,'o-',label='Static nominal layout');axs[0].plot(dyn.step,dyn.ideal_compensated_coverage*100,'s-',label='Perfect known-bias compensation');axs[0].set(xlabel='Synthetic scenario step',ylabel='Coverage (%)');axs[0].legend();label(axs[0],'a  Idealized time-varying bias')
    axs[1].plot(dyn.step,dyn.cumulative_nominal_movement_m,'o-',c=colors[1]);axs[1].set(xlabel='Synthetic scenario step',ylabel='Cumulative movement of all 9 nodes (m)');label(axs[1],'b  Movement proxy, not measured energy');fig.tight_layout();finish(fig,'figS2_dynamic_movement')
    fig=plt.figure(figsize=(9,3.8))
    for i,(name,s) in enumerate([('Fixed',layout()),('Optimized',pd.read_csv(D/'optimized_coordinates.csv')[['x','y','z']].to_numpy())]):
        ax=fig.add_subplot(1,2,i+1,projection='3d');ax.scatter(*s.T,s=38,c=colors[i]);ax.set(xlim=(-120,120),ylim=(-120,120),zlim=(175,115),xlabel='x (m)',ylabel='y (m)',zlabel='Depth (m)');label(ax,('a  ' if i==0 else 'b  ')+name+' nominal deployment')
        for z in [137.5,162.5]:ax.plot([-50,50,50,-50,-50],[-10,-10,10,10,-10],[z]*5,color='gray',alpha=.6)
    fig.tight_layout();finish(fig,'figS3_deployment')
    print('Generated 9 vector PDF / 320 dpi PNG figure pairs')
if __name__=='__main__':main()
