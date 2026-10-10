"""Figures from downloaded real currents and completed replay output only."""
from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from make_figures import finish,label

ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data';O=D/'ocean';F=ROOT/'results'/'figures'
SITES=['south_china_sea','kuroshio','north_atlantic']
NAMES=['South China Sea','Kuroshio vicinity','North Atlantic']
COLORS={'Fixed':'#777777','Synthetic-optimum':'#d95f02','Powell':'#2166ac','PSO':'#1b9e77','AFSA':'#7570b3','Random':'#b2182b'}

def main():
    plt.rcParams.update({'font.size':10.5,'axes.labelsize':10.5,'axes.titlesize':11,'legend.fontsize':8.5})
    currents=pd.read_csv(D/'cmems'/'currents_150m.csv',parse_dates=['time'])
    daily=pd.read_csv(O/'daily_test_coverage.csv',parse_dates=['time'])
    runs=pd.read_csv(O/'optimizer_runs.csv');trace=pd.read_csv(O/'optimizer_traces.csv')
    summary=pd.read_csv(O/'summary.csv');sens=pd.read_csv(O/'sensitivity.csv')
    fig,axes=plt.subplots(2,3,figsize=(10,6.8))
    for j,(site,name) in enumerate(zip(SITES,NAMES)):
        c=currents[currents.site==site]
        for split,col,alpha in [('train','#8fb7d8',.45),('test','#2166ac',.7)]:
            x=c[c.split==split];axes[0,j].scatter(x.u_m_s,x.v_m_s,s=7,color=col,alpha=alpha,
                    label='2018 training' if split=='train' else '2019 held out')
        axes[0,j].axhline(0,color='#aaaaaa',lw=.5);axes[0,j].axvline(0,color='#aaaaaa',lw=.5)
        axes[0,j].set(xlabel='Eastward velocity (m/s)',ylabel='Northward velocity (m/s)')
        axes[0,j].legend(loc='best',markerscale=1.7);label(axes[0,j],f'{chr(97+j)}  {name}')
        for layout in ['Fixed','Synthetic-optimum','Powell']:
            x=daily[(daily.site==site)&(daily.layout==layout)]
            axes[1,j].plot(x.time,x.coverage*100,color=COLORS[layout],lw=.65,
                          label={'Synthetic-optimum':'Synthetic optimum','Powell':'Current-aware Powell'}.get(layout,layout))
        axes[1,j].set(xlabel='Held-out dates in 2019',ylabel='Daily expected coverage (%)')
        axes[1,j].set_xlim(pd.Timestamp('2019-01-01'),pd.Timestamp('2019-12-31'))
        axes[1,j].xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,4,7,10]))
        axes[1,j].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        label(axes[1,j],f'{chr(100+j)}  {name}')
    handles,labels=axes[1,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=9)
    fig.tight_layout(pad=1.3,rect=(0,.045,1,1));finish(fig,'fig7_real_ocean')
    fig,axes=plt.subplots(3,2,figsize=(10,8.5))
    for j,(site,name) in enumerate(zip(SITES,NAMES)):
        c=currents[currents.site==site]
        axes[j,0].plot(c.time,c.u_m_s,lw=.65,label='Eastward',color='#2166ac')
        axes[j,0].plot(c.time,c.v_m_s,lw=.65,label='Northward',color='#d95f02')
        axes[j,0].axvline(pd.Timestamp('2019-01-01'),color='#333333',ls='--',lw=.8,label='Train/test split')
        axes[j,0].set(xlabel='Date',ylabel='Real reanalysis velocity (m/s)');axes[j,0].legend(ncol=3,fontsize=8)
        axes[j,0].tick_params(axis='x',rotation=20);label(axes[j,0],f'{chr(97+2*j)}  {name}')
        for split,col in [('train','#8fb7d8'),('test','#2166ac')]:
            x=c[c.split==split];speed=np.hypot(x.u_m_s,x.v_m_s)
            axes[j,1].hist(speed,bins=22,density=True,histtype='step',lw=1.3,color=col,label=split)
        axes[j,1].set(xlabel='Current speed at 150 m (m/s)',ylabel='Density (s/m)');axes[j,1].legend()
        label(axes[j,1],f'{chr(98+2*j)}  {name}')
    fig.tight_layout(pad=1.4);finish(fig,'figS5_ocean_forcing')
    fig,axes=plt.subplots(2,3,figsize=(10,6.5))
    for j,(site,name) in enumerate(zip(SITES,NAMES)):
        for layout in ['Fixed','Synthetic-optimum','Powell']:
            a=sens[(sens.site==site)&(sens.layout==layout)&(sens.sigma_h_m==120)].sort_values('horizon_s')
            axes[0,j].plot(a.horizon_s/60,a.coverage*100,'o-',color=COLORS[layout],label=layout)
            a=sens[(sens.site==site)&(sens.layout==layout)&(sens.horizon_s==600)].sort_values('sigma_h_m')
            axes[1,j].plot(a.sigma_h_m,a.coverage*100,'o-',color=COLORS[layout],label=layout)
        axes[0,j].set(xlabel='Assumed drift horizon (min)',ylabel='Held-out mean coverage (%)')
        axes[1,j].set(xlabel='Assumed horizontal residual scale (m)',ylabel='Held-out mean coverage (%)')
        label(axes[0,j],f'{chr(97+j)}  {name}');label(axes[1,j],f'{chr(100+j)}  {name}')
    axes[0,0].legend(fontsize=8);axes[1,0].legend(fontsize=8)
    fig.tight_layout(pad=1.2);finish(fig,'figS6_ocean_sensitivity')
    fig,axes=plt.subplots(3,3,figsize=(10,8.4));methods=['Powell','PSO','AFSA','Random']
    for i,(site,name) in enumerate(zip(SITES,NAMES)):
        for method in methods:
            tr=trace[(trace.site==site)&(trace.method==method)]
            arrays=[]
            for _,x in tr.groupby('seed'):
                v=x.sort_values('evaluation').best.to_numpy();arrays.append(np.pad(v,(0,300-len(v)),constant_values=v[-1]))
            a=np.array(arrays)*100;x=np.arange(1,301);lo,med,hi=np.quantile(a,[.25,.5,.75],axis=0)
            axes[i,0].plot(x,med,color=COLORS[method],label=method);axes[i,0].fill_between(x,lo,hi,color=COLORS[method],alpha=.15)
        a=runs[runs.site==site]
        for j,col in enumerate(['training_refined','seconds'],1):
            vals=[a[a.method==m][col].to_numpy()*(100 if col=='training_refined' else 1) for m in methods]
            axes[i,j].boxplot(vals,tick_labels=methods,showfliers=False)
            rng=np.random.default_rng(9410+i+j)
            for k,(m,v) in enumerate(zip(methods,vals),1):axes[i,j].scatter(k+rng.uniform(-.10,.10,len(v)),v,s=9,color=COLORS[m],alpha=.75)
            axes[i,j].tick_params(axis='x',rotation=20)
        axes[i,0].set(xlabel='Counted training-objective calls',ylabel='Training coverage (%)')
        axes[i,1].set(ylabel='Refined training coverage (%)');axes[i,2].set(ylabel='Optimization runtime (s)')
        for j in range(3):label(axes[i,j],f'{chr(97+3*i+j)}  {name}')
    axes[0,0].legend(fontsize=8,ncol=2);fig.tight_layout(pad=1.25);finish(fig,'figS7_ocean_algorithms')
    print('Generated 4 real-current PDF/PNG figure pairs from archived outputs')

if __name__=='__main__':main()
