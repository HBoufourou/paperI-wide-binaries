"""Single-process resumable energy grid; no heldout-dependent tuning."""
import os
for _name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[_name]='1'
from pathlib import Path
import argparse,json,hashlib,time
import numpy as np
from interaction_energy import interaction_energy
from interaction_model import write_model
HERE=Path(__file__).resolve().parent


def run(refined=False):
    qnodes=np.array([.1,.3,1.]);nr,nmu=(33,17) if refined else (17,9)
    x=np.linspace(np.log(.03),np.log(30.),nr);mu=np.linspace(-1,1,nmu)
    values=np.empty((2,3,nr,nmu));cache=HERE/'results'/'energies';cache.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();done=0
    source_sha=hashlib.sha256((HERE/'interaction_energy.py').read_bytes()).hexdigest()
    for iq,q in enumerate(qnodes):
        for ir,xx in enumerate(x):
            for im,mm in enumerate(mu):
                key=f'q{q:g}_x{xx:.12f}_mu{mm:.12f}'
                output=cache/(key+'.json')
                if output.exists():
                    result=json.loads(output.read_text())
                    expected=dict(q=float(q),r=float(np.exp(xx)),mu=float(mm),epsilon=.00125,
                                  level='base',outer_factor=64.,weak_threshold=.01)
                    if result.get('source_sha256')!=source_sha:
                        raise RuntimeError(f'Stale energy source: {output}')
                    for field,value in expected.items():
                        got=result.get(field)
                        compatible=(got==value) if isinstance(value,str) else (got is not None and abs(got-value)<=1e-12*max(1,abs(value)))
                        if not compatible:raise RuntimeError(f'Energy cache parameter mismatch {field}: {output}')
                else:
                    result=interaction_energy(float(q),float(np.exp(xx)),float(mm))
                    result['source_sha256']=source_sha
                    output.write_text(json.dumps(result,indent=2),encoding='utf-8')
                values[:,iq,ir,im]=[-np.exp(xx)*result['Newton_Psi'],-np.exp(xx)*result['Psi']]
                done+=1
                if done%27==0:print(json.dumps(dict(done=done,total=3*nr*nmu,seconds=time.perf_counter()-start)),flush=True)
    name='refined' if refined else 'initial'
    metadata=dict(stage=name,external_newtonian=1.,epsilon=.00125,source_sha256=hashlib.sha256((HERE/'interaction_energy.py').read_bytes()).hexdigest(),
       origin='Psi(infinity)=0',coefficient_basis='unscaled powers (lnr-lnri)^p (mu-muj)^q',
       physical_scope='two prescribed finite Plummer sources; discrete mass ratios; numeric validation is pointwise',
       units='G=Mtotal=a0=1',radius_min=.03,radius_max=30.)
    path=HERE/'results'/f'model_{name}.npz';write_model(path,x,mu,qnodes,values,metadata)
    print(json.dumps(dict(complete=True,path=str(path),seconds=time.perf_counter()-start)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--refined',action='store_true');args=p.parse_args();run(args.refined)
