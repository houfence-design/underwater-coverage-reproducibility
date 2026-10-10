"""Compile the full multi-file journal project with an installed TeX Live."""
from pathlib import Path
import subprocess,shutil,json
R=Path(__file__).resolve().parents[1];P=R/'paper';logs=P/'build_logs';logs.mkdir(exist_ok=True)
bin=Path('D:/texlive/2026/bin/windows')
def run(program,args,tag):
    exe=shutil.which(program) or str(bin/(program+'.exe'))
    result=subprocess.run([exe,*args],cwd=P,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (logs/(tag+'.txt')).write_bytes(result.stdout)
    if result.returncode:raise RuntimeError(f'{program} failed: inspect {logs/(tag+".txt")}')
for name in ['main','supplementary','response_to_reviewers','formula_code_companion']:
    if not (P/(name+'.tex')).exists():continue
    run('pdflatex',['-interaction=nonstopmode','-halt-on-error',name+'.tex'],name+'_pass1')
    if name=='main':run('bibtex',[name],name+'_bibtex')
    for k in [2,3]:run('pdflatex',['-interaction=nonstopmode','-halt-on-error',name+'.tex'],name+f'_pass{k}')
    print('Compiled',name+'.pdf',flush=True)
