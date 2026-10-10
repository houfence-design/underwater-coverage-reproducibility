"""Download the explicitly identified CMEMS product mentioned generically in the original.

Run `copernicusmarine login` locally first. Never put credentials in this file.
The frozen NetCDF subsets and CSV extraction will be distributed for offline
reproduction under the original Copernicus Marine data terms.
"""
from pathlib import Path
import argparse,hashlib,json,datetime
import copernicusmarine as cm
import xarray as xr
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
# Optional native-resolution alternative. It must not overwrite the frozen
# public-archive records that underlie the manuscript's 2018/2019 replay.
OUT=ROOT/'data'/'cmems_native_optional'
DATASET='cmems_mod_glo_phy_my_0.083deg_P1D-m'
VERSION='202311'
SITES={'south_china_sea':(115.,20.),'kuroshio':(130.,25.),'north_atlantic':(-45.,40.)}

def requests():
    return {name:dict(dataset_id=DATASET,dataset_version=VERSION,
        variables=['uo','vo'],minimum_longitude=lon,maximum_longitude=lon,
        minimum_latitude=lat,maximum_latitude=lat,
        minimum_depth=100.,maximum_depth=200.,
        start_datetime='2022-01-01T00:00:00',end_datetime='2023-12-31T23:59:59',
        coordinates_selection_method='nearest',service='arco-time-series',
        output_filename=name+'.nc',file_format='netcdf',
        netcdf_compression_level=1,disable_progress_bar=True)
        for name,(lon,lat) in SITES.items()}

def extract():
    frames=[];details=[]
    for name,(lon,lat) in SITES.items():
        path=OUT/(name+'.nc')
        ds=xr.open_dataset(path)
        assert ds.uo.attrs['units'] in ['m s-1','m/s','m s^-1']
        assert ds.vo.attrs['units'] in ['m s-1','m/s','m s^-1']
        assert ds.sizes['longitude']==ds.sizes['latitude']==1
        depths=ds.depth.values
        assert depths.min()<=150<=depths.max()
        # Coordinate-only selection and linear depth interpolation. No smoothing,
        # calibration, replacement or interpolation across missing dates.
        sub=ds[['uo','vo']].interp(depth=150.).squeeze(['longitude','latitude'])
        frame=pd.DataFrame({'site':name,'time':pd.to_datetime(sub.time.values),
            'u_m_s':sub.uo.values,'v_m_s':sub.vo.values})
        frame['valid']=np.isfinite(frame.u_m_s)&np.isfinite(frame.v_m_s)
        frame['split']=np.where(frame.time.dt.year==2022,'train','test')
        frame['longitude']=float(ds.longitude.values[0]);frame['latitude']=float(ds.latitude.values[0])
        frame['depth_m']=150.
        assert frame.time.is_unique
        assert len(frame)==730
        assert np.all(np.diff(frame.time.to_numpy())==np.timedelta64(1,'D'))
        frames.append(frame)
        details.append({'site':name,'requested_longitude':lon,'requested_latitude':lat,
            'actual_longitude':float(ds.longitude.values[0]),'actual_latitude':float(ds.latitude.values[0]),
            'native_depths_m':depths.tolist(),'analysis_depth_m':150.,
            'interpolation':'linear in depth; nearest horizontal grid cell',
            'valid_train':int(frame.loc[frame.split=='train','valid'].sum()),
            'valid_test':int(frame.loc[frame.split=='test','valid'].sum()),
            'sha256_netcdf':hashlib.sha256(path.read_bytes()).hexdigest(),
            'dataset_attributes':dict(ds.attrs)})
    all_data=pd.concat(frames,ignore_index=True)
    all_data.to_csv(OUT/'currents_150m.csv',index=False,float_format='%.12g')
    manifest={'provider':'E.U. Copernicus Marine Service / Mercator Ocean International',
        'product_id':'GLOBAL_MULTIYEAR_PHY_001_030','product_doi':'10.48670/moi-00021',
        'dataset_id':DATASET,'dataset_version':VERSION,'data_type':'GLORYS12V1 ocean reanalysis; not direct in-situ measurements',
        'variables':{'uo':'eastward sea water velocity [m s-1]','vo':'northward sea water velocity [m s-1]'},
        'download_date':'2026-10-09','train_period':'2022-01-01/2022-12-31','test_period':'2023-01-01/2023-12-31',
        'copernicusmarine_version':cm.__version__,'sites':details,
        'sha256_csv':hashlib.sha256((OUT/'currents_150m.csv').read_bytes()).hexdigest(),
        'source_url':'https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description',
        'acknowledgement':'This study has been conducted using E.U. Copernicus Marine Service Information: https://doi.org/10.48670/moi-00021',
        'limits':'Daily means on an approximately 8 km grid; do not resolve turbulent variability within the 100 m target or establish sensor localization errors.'}
    (OUT/'provenance.json').write_text(json.dumps(manifest,indent=2,default=str),encoding='utf-8')
    print(all_data.groupby(['site','split']).valid.agg(['size','sum']).to_string(),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--extract-only',action='store_true');parser.add_argument('--requests-only',action='store_true');args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    reqs=requests();(OUT/'download_requests.json').write_text(json.dumps(reqs,indent=2),encoding='utf-8')
    if args.requests_only:return
    if not args.extract_only:
        catalogue=cm.describe(dataset_id=DATASET,disable_progress_bar=True)
        (OUT/'catalogue.json').write_text(catalogue.model_dump_json(indent=2),encoding='utf-8')
        for name,request in reqs.items():
            if (OUT/(name+'.nc')).exists():continue
            response=cm.subset(**request,output_directory=OUT)
            (OUT/(name+'_download_response.json')).write_text(response.model_dump_json(indent=2),encoding='utf-8')
    extract()

if __name__=='__main__':main()
