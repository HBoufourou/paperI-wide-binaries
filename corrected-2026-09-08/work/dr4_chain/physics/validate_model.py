"""Frozen development/heldout force checks and energy convergence controls."""
import os
for _name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[_name]='1'
from pathlib import Path
import sys,json,argparse,hashlib,shutil
import numpy as np
from interaction_model import InteractionModel
from interaction_energy import interaction_energy
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'dr4_binary'/'physics'))
from binary_qumond import solve_case
HERE=Path(__file__).resolve().parent
DEVELOPMENT=[(.055,-.83),(.18,.28),(.55,-.37),(1.4,.68),(4.7,-.55),(17,.1)]
FINAL=[(.045,.61),(.12,-.22),(.41,.87),(.9,-.71),(2.8,.19),(11.5,-.46),(24,.79)]

def vec(r,mu):return r*np.array([np.sqrt(1-mu*mu),0,mu])
def relative(a,b):return float(np.linalg.norm(np.array(a)-b)/np.linalg.norm(b))

def direct(q,r,mu,role):
    path=HERE/'results'/'forces'/f'{role}_q{q:g}_r{r:g}_mu{mu:g}.json';path.parent.mkdir(parents=True,exist_ok=True)
    source_path=HERE.parents[1]/'dr4_binary'/'physics'/'binary_qumond.py'
    source_sha=hashlib.sha256(source_path.read_bytes()).hexdigest()
    if path.exists():
        out=json.loads(path.read_text())
        if out.get('producer_sha256')!=source_sha:raise RuntimeError(f'Stale force producer: {path}')
        case=out['case']
        if abs(case['q']-q)>1e-12 or abs(case['separation']-r)>1e-12 or abs(case['theta_degrees']-np.degrees(np.arccos(mu)))>1e-10:
            raise RuntimeError(f'Force cache parameters differ: {path}')
        return out
    out=solve_case(q,r,float(np.degrees(np.arccos(mu))),level='fine')
    out['producer_sha256']=source_sha
    path.write_text(json.dumps(out,indent=2),encoding='utf-8');return out

def check_forces(model,points,role):
    rows=[]
    for q in [.1,.3,1.]:
        for r,mu in points:
            data=direct(q,r,mu,role);m=np.array(data['masses']);f=np.array(data['newton_force'])
            truth={'qumond':np.array(data['relative_acceleration']),'newton':f[1]/m[1]-f[0]/m[0]}
            for gravity in ['newton','qumond']:
                pos=vec(r,mu);acc=model.acceleration(pos,q,gravity)
                row=dict(q=q,r=r,mu=mu,gravity=gravity,relative_force_error=relative(acc,truth[gravity]),
                         interpolated_acceleration=acc.tolist(),reference_acceleration=truth[gravity].tolist())
                if role=='final':
                    step=1e-5*r
                    fd=np.array([-(model.psi(pos+step*axis,q,gravity)-model.psi(pos-step*axis,q,gravity))/(2*step) for axis in np.eye(3)])
                    h=1e-4*r
                    jac=np.column_stack([(model.acceleration(pos+h*axis,q,gravity)-model.acceleration(pos-h*axis,q,gravity))/(2*h) for axis in np.eye(3)])
                    row['finite_difference_gradient_error']=relative(fd,acc)
                    row['jacobian_asymmetry']=float(np.linalg.norm(jac-jac.T)/np.linalg.norm(jac))
                rows.append(row)
            print(json.dumps(dict(role=role,q=q,r=r,mu=mu,error=rows[-1]['relative_force_error'])),flush=True)
    return rows

def energy_controls():
    rows=[]
    for q in [.3,1.]:
        for r in [.05,1.,20.]:
            saved=[]
            for name,kwargs in [('base',{}),('reference',{'level':'reference'}),('outer',{'outer_factor':128.})]:
                path=HERE/'results'/'energy_controls'/f'q{q:g}_r{r:g}_{name}.json';path.parent.mkdir(parents=True,exist_ok=True)
                if path.exists():
                    out=json.loads(path.read_text())
                    if out.get('producer_sha256')!=hashlib.sha256((HERE/'interaction_energy.py').read_bytes()).hexdigest():
                        raise RuntimeError(f'Stale energy control: {path}')
                else:
                    out=interaction_energy(q,r,-.6,**kwargs)
                    out['producer_sha256']=hashlib.sha256((HERE/'interaction_energy.py').read_bytes()).hexdigest()
                    path.write_text(json.dumps(out,indent=2))
                saved.append(out)
            a,b,c=saved
            rows.append(dict(q=q,r=r,mu=-.6,mesh_relative=abs(a['Psi']/b['Psi']-1),outer_relative=abs(a['Psi']/c['Psi']-1)))
    return rows

def run(stage,model_path):
    model=InteractionModel.load(model_path,require_validation=False)
    rows=check_forces(model,DEVELOPMENT if stage=='development' else FINAL,stage)
    maximum=max(row['relative_force_error'] for row in rows)
    out=dict(stage=stage,model_path=str(model_path),max_force_error=maximum,rows=rows)
    if stage=='development':out['refine_required']=maximum>.005
    else:
        controls=energy_controls();out['energy_controls']=controls
        old=[]
        for q in [.1,.3,1.]:
            p=HERE.parents[1]/'dr4_binary'/'physics'/'results'/'jobs'/f'reference_E1_q{q:g}_d1_t45_reference.json'
            d=json.loads(p.read_text());acc=model.acceleration(vec(1.,1/np.sqrt(2)),q)
            old.append(dict(q=q,relative_error=relative(acc,d['relative_acceleration'])))
        out['old_reference_comparisons']=old
        out['passed']=bool(all(row['relative_force_error']<.01 and row['finite_difference_gradient_error']<1e-5 and row['jacobian_asymmetry']<1e-4 for row in rows)
            and all(row['mesh_relative']<.001 and row['outer_relative']<.001 for row in controls))
        out['scope']='pointwise force/energy validation; no uniform physical accuracy certificate'
        target=HERE/'results'/'interaction_model.npz';shutil.copyfile(model_path,target)
        out['model_sha256']=hashlib.sha256(target.read_bytes()).hexdigest()
        out['evaluator_sha256']=hashlib.sha256((HERE/'interaction_model.py').read_bytes()).hexdigest()
        out['cmax']={g:[float(model.cmax(q,g)) for q in [.1,.3,1.]] for g in ['newton','qumond']}
        out['cmin']={g:[float(model.cmin(q,g)) for q in [.1,.3,1.]] for g in ['newton','qumond']}
        out['outer_amax']={g:[float(model.outer_amax(q,g)) for q in [.1,.3,1.]] for g in ['newton','qumond']}
        (HERE/'validation.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    name=f'{stage}_{Path(model_path).stem}.json';(HERE/'results'/name).write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k not in ['rows','energy_controls']},indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','final']);p.add_argument('model_path');a=p.parse_args();run(a.stage,Path(a.model_path))
