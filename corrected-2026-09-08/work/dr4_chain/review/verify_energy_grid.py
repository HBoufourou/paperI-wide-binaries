"""Check the scalar-grid values against every relevant saved action integral."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--model',type=Path,required=True)
    parser.add_argument('--energy-root',type=Path);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();folder=args.energy_root or args.model.parent/'energies'
    with np.load(args.model,allow_pickle=False) as z:
        x=z['log_r'];mu=z['mu_nodes'];q=z['q_nodes'];a=z['A_values'];meta=json.loads(str(z['metadata']))
    values=np.full(a.shape,np.nan);seen=np.zeros(a.shape[1:],bool);records=[];failures=[]
    for path in sorted(folder.glob('*.json')):
        row=json.loads(path.read_text());iq=int(np.argmin(abs(q-row['q'])));ir=int(np.argmin(abs(x-np.log(row['r']))));im=int(np.argmin(abs(mu-row['mu'])))
        if max(abs(q[iq]-row['q']),abs(x[ir]-np.log(row['r'])),abs(mu[im]-row['mu']))>1e-11:continue
        if seen[iq,ir,im]:failures.append('duplicate/'+path.name)
        seen[iq,ir,im]=True
        reduced=row['q']/(1+row['q'])**2
        if abs(row['reduced_mass']/reduced-1)>1e-14:failures.append('reduced_mass/'+path.name)
        if row['source_sha256']!=meta['source_sha256']:failures.append('source_hash/'+path.name)
        values[:,iq,ir,im]=-np.exp(x[ir])*np.array([row['Newton_U'],row['U']])/reduced
        records.append({'file':path.name,'sha256':sha(path)})
    err=float(np.max(abs(values-a))) if seen.all() else None
    out={'model_sha256':sha(args.model),'source_sha256':sha(__file__),
         'grid_points_checked':len(records),'expected_grid_points':int(seen.size),
         'missing_points':np.argwhere(~seen).tolist(),'failures':failures,
         'A_max_absolute_difference':err,'passed':bool(seen.all() and not failures and err<1e-12),
         'scope':'Reconstruction from saved integrals, not repeated numerical integration.',
         'records':records}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='records'}),flush=True)

if __name__=='__main__':main()
