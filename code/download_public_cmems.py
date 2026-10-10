"""Extract real GLORYS12V1 currents from Mercator Ocean's public archive.

Provider documentation: mercator-ocean/oceanbench, reference.py,
glorys_reanalysis_1_degree. No Copernicus account is needed for this explicitly
public EDITO archive. This is the provider's 1-degree derivative, not native
1/12-degree data. Original metadata, all request hashes and small native-depth
subsets are preserved; the published CSV can be used completely offline.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,hashlib,json,threading,time,datetime
import requests,numpy as np,pandas as pd,xarray as xr
from numcodecs import get_codec

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data'/'cmems'
BASE='https://s3.waw3-1.cloudferro.com/oceanbench-bucket/public/glorys_1degree_1993_2019'
SOURCE='https://github.com/mercator-ocean/oceanbench/blob/main/oceanbench/datasets/reference.py'
SITES={'south_china_sea':(115.5,20.5),'kuroshio':(130.5,25.5),'north_atlantic':(-44.5,40.5)}
VARS={'eastward_sea_water_velocity':'uo','northward_sea_water_velocity':'vo'}
local=threading.local()

def fetch(url,path):
    if path.exists():return path.read_bytes()
    if not hasattr(local,'session'):local.session=requests.Session()
    for attempt in range(4):
        try:
            response=local.session.get(url,timeout=(15,90));response.raise_for_status()
            path.parent.mkdir(parents=True,exist_ok=True)
            temporary=path.with_suffix(path.suffix+'.partial');temporary.write_bytes(response.content)
            temporary.replace(path);return response.content
        except requests.RequestException:
            if attempt==3:raise
            time.sleep(2**attempt)

def decode(raw,spec):
    assert spec['zarr_format']==2 and spec['filters'] is None
    assert spec['compressor']['id']=='blosc' and spec['order']=='C'
    return np.frombuffer(get_codec(spec['compressor']).decode(raw),dtype=np.dtype(spec['dtype']))

def download(cache):
    OUT.mkdir(parents=True,exist_ok=True);metadata=OUT/'public_archive_metadata';metadata.mkdir(exist_ok=True)
    all_requests=[];site_arrays={name:[] for name in SITES};site_details=None
    dates=pd.date_range('2018-01-01','2019-12-01',freq='MS')
    for month in dates:
        key=month.strftime('%Y%m');base=f'{BASE}/{key}.zarr'
        raw=fetch(base+'/.zmetadata',metadata/(key+'_zmetadata.json'));m=json.loads(raw)['metadata']
        assert m['.zattrs']['source']=='MERCATOR GLORYS12V1'
        coordinates={}
        for axis in ['time','depth','latitude','longitude']:
            url=base+'/'+axis+'/0';raw=fetch(url,metadata/(key+'_'+axis+'.bin'))
            spec=m[axis+'/.zarray'];coordinates[axis]=decode(raw,spec)[:spec['shape'][0]]
            all_requests.append({'url':url,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)})
        lat=coordinates['latitude'];lon=coordinates['longitude'];depth=coordinates['depth']
        below=np.flatnonzero(depth<=150)[-1];above=np.flatnonzero(depth>=150)[0]
        assert below!=above
        time_attrs=m['time/.zattrs'];units=time_attrs['units'];unit,origin=units.split(' since ')
        times=pd.Timestamp(origin)+pd.to_timedelta(coordinates['time'],unit={'hours':'h','days':'D','seconds':'s'}[unit])
        assert len(times)==month.days_in_month and times[0].date()==month.date()
        cells={name:(int(np.argmin(abs(lat-y))),int(np.argmin(abs(lon-x)))) for name,(x,y) in SITES.items()}
        current_details={name:{'requested_longitude':xy[0],'requested_latitude':xy[1],
            'actual_longitude':float(lon[cells[name][1]]),'actual_latitude':float(lat[cells[name][0]]),
            'native_depths_m':[float(depth[below]),float(depth[above])]} for name,xy in SITES.items()}
        if site_details is not None:assert current_details==site_details
        site_details=current_details;values={v:np.empty((len(times),2,len(SITES))) for v in VARS}
        tasks=[(v,j,d,k) for v in VARS for j,d in enumerate([below,above]) for k in range(len(times))]
        def chunk(task):
            variable,j,d,k=task;url=f'{base}/{variable}/{k}.{d}.0.0'
            spec=m[variable+'/.zarray'];assert spec['chunks']==[1,1,170,360]
            raw=fetch(url,cache/key/variable/f'{k}.{d}.bin')
            a=decode(raw,spec).reshape(170,360)
            picks=[float(a[cells[name]]) for name in SITES]
            return task,picks,{'url':url,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}
        with ThreadPoolExecutor(max_workers=8) as pool:
            for (variable,j,d,k),picks,entry in pool.map(chunk,tasks):
                values[variable][k,j]=picks;all_requests.append(entry)
        for i,name in enumerate(SITES):
            ds=xr.Dataset({short:(('time','depth'),values[v][:,:,i],{'units':m[v+'/.zattrs']['units'],
                'standard_name':m[v+'/.zattrs'].get('standard_name',v)}) for v,short in VARS.items()},
                coords={'time':times,'depth':depth[[below,above]]},attrs=m['.zattrs'])
            site_arrays[name].append(ds)
        print('Downloaded public CMEMS archive',key,'days',len(times),flush=True)
        (OUT/'public_download_requests.json').write_text(json.dumps(all_requests,indent=2),encoding='utf-8')
    frames=[];details=[]
    for name,datasets in site_arrays.items():
        ds=xr.concat(datasets,dim='time');info=site_details[name];path=OUT/(name+'.nc')
        ds.attrs.update({'public_archive_source':BASE,'longitude':info['actual_longitude'],'latitude':info['actual_latitude'],
                        'archive_resolution_degrees':1.0,'extraction_note':'Two original depth layers retained; no time interpolation'})
        ds.to_netcdf(path,engine='h5netcdf',encoding={v:{'zlib':True,'complevel':4} for v in ['uo','vo']})
        sub=ds.interp(depth=150.);frame=pd.DataFrame({'site':name,'time':sub.time.values,
            'u_m_s':sub.uo.values,'v_m_s':sub.vo.values})
        frame['valid']=np.isfinite(frame.u_m_s)&np.isfinite(frame.v_m_s)
        frame['split']=np.where(pd.to_datetime(frame.time).dt.year==2018,'train','test')
        frame['longitude']=info['actual_longitude'];frame['latitude']=info['actual_latitude'];frame['depth_m']=150.
        assert len(frame)==730 and frame.valid.all() and frame.time.is_unique
        assert np.all(np.diff(frame.time.to_numpy())==np.timedelta64(1,'D'))
        frames.append(frame);details.append(dict(site=name,**info,analysis_depth_m=150.,
            interpolation='provider 1-degree derivative, exact horizontal cell selection; linear depth interpolation',
            valid_train=365,valid_test=365,sha256_netcdf=hashlib.sha256(path.read_bytes()).hexdigest()))
    combined=pd.concat(frames,ignore_index=True);csv=OUT/'currents_150m.csv'
    combined.to_csv(csv,index=False,float_format='%.12g')
    manifest={'provider':'Mercator Ocean International / E.U. Copernicus Marine Service via public OceanBench EDITO archive',
        'product_id':'GLOBAL_MULTIYEAR_PHY_001_030','historical_product_alias':'GLOBAL_REANALYSIS_PHY_001_030',
        'product_doi':'10.48670/moi-00021','dataset_id':'OceanBench GLORYS12V1 1-degree historical derivative',
        'dataset_version':'Upstream version not specified in public archive; exact request hashes are frozen',
        'data_type':'GLORYS12V1 ocean reanalysis derivative; not direct in-situ measurements',
        'native_product_resolution_degrees':1/12,'downloaded_archive_resolution_degrees':1.,
        'variables':{'uo':'eastward sea water velocity [m s-1]','vo':'northward sea water velocity [m s-1]'},
        'download_date':datetime.date.today().isoformat(),'train_period':'2018-01-01/2018-12-31','test_period':'2019-01-01/2019-12-31',
        'sites':details,'sha256_csv':hashlib.sha256(csv.read_bytes()).hexdigest(),
        'source_url':'https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description',
        'public_access_documentation':SOURCE,'public_archive_url':BASE,
        'original_manuscript_link':'The original names CMEMS only; it does not identify a specific product, date, station or raw file.',
        'acknowledgement':'This study has been conducted using E.U. Copernicus Marine Service Information: https://doi.org/10.48670/moi-00021. Public historical subsets are distributed by Mercator Ocean International through OceanBench.',
        'limits':'Daily 1-degree derivative (~100 km). Represents large-scale environmental forcing; cannot resolve 100 m target turbulence or calibrate sensor residual positioning errors.'}
    (OUT/'provenance.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(combined.groupby(['site','split']).valid.agg(['size','sum']).to_string(),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,default=ROOT/'.download_cache')
    args=p.parse_args();download(args.cache)
