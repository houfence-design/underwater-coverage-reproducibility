"""One entry point for offline verification, figures or complete reruns."""
from pathlib import Path
import argparse,subprocess,sys,zipfile
ROOT=Path(__file__).resolve().parents[1]

def unpack_if_needed():
    if (ROOT/'data'/'summary.json').exists():return
    archive=ROOT/'frozen_data.zip'
    if not archive.exists():raise FileNotFoundError('Download frozen_data.zip from the same repository release')
    with zipfile.ZipFile(archive) as z:
        for member in z.infolist():
            destination=(ROOT/member.filename).resolve()
            if not destination.is_relative_to(ROOT.resolve()):raise ValueError('Unsafe archive path')
            if member.is_dir():continue
            destination.parent.mkdir(parents=True,exist_ok=True)
            destination.write_bytes(z.read(member))

def command(name,*args):
    print('Executing',name,flush=True)
    subprocess.run([sys.executable,str(ROOT/'code'/name),*args],cwd=ROOT,check=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['verify','figures','full'],default='verify')
    p.add_argument('--redownload-public-currents',action='store_true')
    args=p.parse_args();unpack_if_needed()
    if args.redownload_public_currents:command('download_public_cmems.py')
    if args.mode=='full':
        command('run_experiments.py');command('run_ocean_experiments.py')
    command('test_formulas.py');command('verify_results.py');command('verify_ocean_results.py')
    if args.mode in ['figures','full']:
        command('make_figures.py');command('prepare_tables.py')
        command('make_ocean_figures.py');command('prepare_ocean_tables.py');command('build_formula_inventory.py')
    print('Requested reproduction mode completed',flush=True)
if __name__=='__main__':main()
