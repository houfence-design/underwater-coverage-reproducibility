"""Check actual bytes against the delivered source/archive manifest."""
from pathlib import Path
import json,hashlib,sys
ROOT=Path(__file__).resolve().parents[1]

def main():
    entries=json.loads((ROOT/'PUBLIC_FILE_MANIFEST.json').read_text(encoding='utf-8'))
    seen=set();problems=[]
    for entry in entries:
        name=entry['path'];path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT.resolve()) or name in seen:
            problems.append({'path':name,'problem':'unsafe or duplicate path'});continue
        seen.add(name)
        if not path.is_file():
            problems.append({'path':name,'problem':'missing file'});continue
        raw=path.read_bytes()
        if hashlib.sha256(raw).hexdigest()!=entry['sha256']:
            problems.append({'path':name,'problem':'SHA-256 mismatch'})
        if len(raw)!=entry['bytes']:
            problems.append({'path':name,'problem':'byte-length mismatch'})
    print(json.dumps({'entries':len(entries),'passed':len(entries)-len({x['path'] for x in problems}),'problems':problems},ensure_ascii=False,indent=2))
    return 1 if problems else 0

if __name__=='__main__':sys.exit(main())
