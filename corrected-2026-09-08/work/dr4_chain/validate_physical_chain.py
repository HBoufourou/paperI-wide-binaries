"""Execute coupled trajectory/observer validation after final potential validation."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from population import generate_library,QS
from orbits import propagate,units,tangent_basis,project
from joint_observer import observe_pairs,cadence,SCENARIOS

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'physics'))
from interaction_model import InteractionModel

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save_json(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')

def main():
    implementation_files=['OBSERVER_VALIDATION.json','ORBIT_IMPLEMENTATION_VALIDATION.json','POPULATION_IMPLEMENTATION_VALIDATION.json']
    for name in implementation_files:
        proof=json.loads((HERE/name).read_text(encoding='utf-8'))
        assert proof['passed'],('Implementation check failed',name)
        for source,digest in proof['source_sha256'].items():
            assert sha(HERE/source)==digest,('Stale implementation proof',name,source)
    model=InteractionModel.load()
    cmax=float(max(np.max(model.cmax(QS,g)) for g in ('newton','qumond')))
    out=HERE/'trajectory_validation';out.mkdir(exist_ok=True)
    records=[];files=[]
    def check(name,value,threshold):
        row=dict(name=name,value=float(value),threshold=threshold,passed=bool(np.isfinite(value) and value<threshold));records.append(row)
        if not row['passed']:raise AssertionError(row)
    for gi,gravity in enumerate(('newton','qumond')):
        potential=lambda x,q:model.potential(x,q,gravity)
        acceleration=lambda x,q:model.acceleration(x,q,gravity)
        p,meta=generate_library(512,268000101+gi*1000,potential,cmax,'comp40')
        path=out/(gravity+'_initial.npz');np.savez_compressed(path,**p);files.append(path)
        basis=tangent_basis(p['ra_deg'],p['dec_deg']);inner={k.removeprefix('inner_'):v for k,v in p.items() if k.startswith('inner_')}
        for scenario in SCENARIOS:
            t,angles,n3=cadence(scenario)
            x1,v1=propagate(p['position'],p['velocity'],p['q'],p['Mtrue_total'],t,acceleration,substeps=1)
            x2,v2=propagate(p['position'],p['velocity'],p['q'],p['Mtrue_total'],t,acceleration,substeps=2)
            energy=.5*np.sum(v1*v1,axis=2)+model.potential(x1,p['q'][:,None],gravity)
            check(gravity+'/'+scenario+'/energy',np.max(abs((energy-p['energy'][:,None])/p['energy'][:,None])),1e-7)
            check(gravity+'/'+scenario+'/Lz',np.max(abs(np.cross(x1,v1)[:,:,2]-p['lz'][:,None])),1e-10)
            check(gravity+'/'+scenario+'/position_refinement_AU',np.max(np.linalg.norm(x2-x1,axis=2)*p['rM_au'][:,None]),1e-5)
            xy=project(x1-p['position'][:,None,:],basis)[:,:,:2]*p['rM_au'][:,None,None]
            observed=observe_pairs(None,p['distance_pc'],p['pm_err_masyr'],inner,scenario,268000201+gi*1000,
                p['ra_deg'],p['dec_deg'],outer_relative_xy_au=xy,mass_fraction=p['mass_fraction'])
            for b in (34,66):
                residual=observed[f'pm{b}_masyr']-observed[f'outer{b}_masyr']-observed[f'inner_bias{b}_masyr']-observed[f'noise{b}_masyr']
                check(gravity+'/'+scenario+f'/joint_decomposition{b}',np.max(abs(residual)),1e-10)
            assert all(np.isfinite(v).all() for v in observed.values())
            cross=observed['cov_cross_relative_masyr2'];delta=observed['cov_delta_relative_masyr2']
            reconstructed=observed['cov34_masyr2']+observed['cov66_masyr2']-cross-cross.swapaxes(-1,-2)
            check(gravity+'/'+scenario+'/shared_covariance',np.max(np.max(abs(delta-reconstructed),axis=(1,2))/np.max(abs(delta),axis=(1,2))),1e-9)
            path=out/(gravity+'_'+scenario+'.npz');np.savez_compressed(path,positions=x1,velocities=v1,positions_refined=x2,**observed);files.append(path)
        candidates=np.flatnonzero((p['barrier_rlo']>.2)&(np.linalg.norm(p['position'],axis=1)<4))
        if len(candidates)<24:raise RuntimeError('Insufficient declared long-trajectory validation states')
        ids=candidates[:24];initial=p['position'][ids];velocity=p['velocity'][ids];q=p['q'][ids]
        total=np.ones(len(ids));t=np.linspace(0,2,2001)*float(units(1.)['tM_year'])
        paths=[];velocities=[];errors=[]
        for substeps in (1,2):
            x,v=propagate(initial,velocity,q,total,t,acceleration,substeps=substeps)
            e=.5*np.sum(v*v,axis=2)+model.potential(x,q[:,None],gravity)
            errors.append(float(np.max(abs((e-p['energy'][ids,None])/p['energy'][ids,None]))))
            paths.append(x);velocities.append(v)
        check(gravity+'/long_energy',errors[1],5e-4)
        check(gravity+'/long_Lz',np.max(abs(np.cross(paths[1],velocities[1])[:,:,2]-p['lz'][ids,None])),1e-9)
        check(gravity+'/long_refinement',errors[1]/errors[0] if errors[0]>1e-9 else 0,.4)
        path=out/(gravity+'_long.npz');np.savez_compressed(path,ids=ids,times_dimensionless=np.linspace(0,2,2001),
            positions_coarse=paths[0],positions_fine=paths[1],velocities_coarse=velocities[0],velocities_fine=velocities[1]);files.append(path)
        print(json.dumps(dict(gravity=gravity,trajectory_controls_passed=True,common_cmax=cmax)),flush=True)
    result=dict(status='passed_for_declared_domain',passed=True,checks=records,
        scalar_model_sha256=sha(model.path),physics_validation_sha256=sha(HERE/'physics/validation.json'),
        common_cmax=cmax,implementation_proof_sha256={name:sha(HERE/name) for name in implementation_files},
        files=[dict(path=p.relative_to(HERE).as_posix(),sha256=sha(p)) for p in files],
        source_sha256={p.name:sha(p) for p in [HERE/'population.py',HERE/'orbits.py',HERE/'joint_observer.py',HERE/'physics/interaction_model.py',
            HERE/'POPULATION_PROTOCOL.md',HERE/'AMENDMENT_01_PHOTOMETRIC_ANCHOR.md',HERE/'TRAJECTORY_PROTOCOL.md',Path(__file__)]},
        orbit_parent='f(E) times invariant E,Lz restrictions in a validated scalar interpolant',
        selected_catalogue_stationary=False,astrophysical_population_validated=False,
        continuum_qumond_error_bound=False,real_gaia_data=False,sky_ready=False)
    save_json(HERE/'physical_validation.json',result)
    print(json.dumps(dict(passed=True,checks=len(records),production_conditionally_supported=True,sky_ready=False)),flush=True)

if __name__=='__main__':main()
