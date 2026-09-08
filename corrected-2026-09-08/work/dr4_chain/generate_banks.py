"""Generate fresh complete coupled libraries only after physical validation."""
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np
from population import generate_library,QS
from orbits import propagate,tangent_basis,project
from joint_observer import observe_pairs,cadence,SCENARIOS,kepler_eccentric_anomaly

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'physics'))
from interaction_model import InteractionModel

GRAVITIES=('Newton','QUMOND');ROLES=('fitting','calibration','test')
COMPONENTS=('pure','comp40','comp80','comp60','comp160')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def write(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')

def declaration():
    model=InteractionModel.load();validation=read(HERE/'physical_validation.json')
    assert validation['passed'] and validation['status']=='passed_for_declared_domain'
    assert validation['scalar_model_sha256']==sha(model.path)
    for name,digest in validation['source_sha256'].items():
        path=HERE/'physics'/name if name=='interaction_model.py' else HERE/name
        assert sha(path)==digest,('Changed validated source',name)
    records=[];streams=[]
    for gi,gravity in enumerate(GRAVITIES):
        for ri,role in enumerate(ROLES):
            for ci,component in enumerate(COMPONENTS):
                if ci>=3 and role!='test':continue
                ps=260000001+gi*1000000+ri*10000+ci*100
                os=263000001+gi*1000000+ri*10000+ci*100
                records.append(dict(gravity=gravity,role=role,component=component,
                    prior_au=None if component=='pure' else int(component[4:]),
                    population_seed=ps,observer_seed=os,initial_N=20000,max_N=40000,
                    folder='banks/'+gravity.lower()+'_'+role+'_'+component))
                streams.extend([ps,os])
    assert len(records)==22 and len(streams)==len(set(streams))==44 and max(streams)<268000000
    sourcefiles=['generate_banks.py','population.py','orbits.py','joint_observer.py','physics/interaction_model.py',
        'POPULATION_PROTOCOL.md','AMENDMENT_01_PHOTOMETRIC_ANCHOR.md','TRAJECTORY_PROTOCOL.md','OBSERVER_PROTOCOL.md']
    result=dict(records=records,streams=streams,common_cmax=validation['common_cmax'],
        physical_validation_sha256=sha(HERE/'physical_validation.json'),model_sha256=sha(model.path),
        source_sha256={name:sha(HERE/name) for name in sourcefiles},
        empirical_source_sha256=sha(HERE.parent/'reanalysis/data/observation_arrays.npz'),
        empirical_selection_source_sha256=sha(HERE.parent/'originality/dr4/run_comparison.py'),
        scenarios=list(SCENARIOS),scope='Conditional complete forward model, not actual Gaia population validation')
    return model,result

def later_resolution_fraction(p,inner,t):
    active=inner['active']
    if not active.any():return dict(unweighted_fraction=0.,weighted_fraction={k:0. for k in ('weight_beta1','weight_beta1p6','weight_beta1p3')})
    a=inner['a_au'];e=inner['e'];period=np.sqrt(a**3/inner['mass_host_true'])
    anomaly=kepler_eccentric_anomaly(inner['M0'][:,None]+2*np.pi*t[None,:]/period[:,None],e[:,None])
    plane=np.stack((a[:,None]*(np.cos(anomaly)-e[:,None]),a[:,None]*np.sqrt(1-e[:,None]**2)*np.sin(anomaly)),axis=2)
    projected=np.einsum('nij,ntj->nti',inner['rotation'][:,:2,:2],plane)
    resolved=np.any(np.linalg.norm(projected,axis=2)>=p['distance_pc'][:,None],axis=1)
    return dict(unweighted_fraction=float(np.mean(resolved[active])),
        weighted_fraction={k:float(np.sum(p[k]*resolved)/np.sum(p[k])) for k in ('weight_beta1','weight_beta1p6','weight_beta1p3')})

def run(freeze_only=False):
    model,spec=declaration();path=HERE/'GENERATION_PROTOCOL.json'
    if path.exists():assert read(path)==spec
    else:write(path,spec)
    if freeze_only:print('GENERATION FROZEN',sha(path),flush=True);return
    manifest_path=HERE/'bank_manifest.json'
    manifest=read(manifest_path) if manifest_path.exists() else dict(schema_version=1,purpose='physics_production',
        production_authorized=False,physical_validation=dict(status='passed_for_declared_domain',path='physical_validation.json',sha256=spec['physical_validation_sha256']),
        physical_scope=spec['scope'],generation_protocol_sha256=sha(path),scenarios=list(SCENARIOS),records=[])
    assert manifest['generation_protocol_sha256']==sha(path)
    for item in spec['records']:
        key=(item['gravity'],item['role'],item['component'])
        previous=next((r for r in manifest['records'] if (r['gravity'],r['role'],r['component'])==key),None)
        if previous:
            assert sha(HERE/previous['population_path'])==previous['sha256']['population']
            for s in SCENARIOS:assert sha(HERE/previous['observations'][s])==previous['sha256']['observations'][s]
            continue
        folder=HERE/item['folder'];folder.mkdir(parents=True,exist_ok=True)
        if list(folder.glob('*.npz')):raise RuntimeError('Incomplete existing bank; preserve and inspect it before any fresh rerun: '+str(folder))
        start=time.perf_counter();gravity=item['gravity'].lower()
        potential=lambda x,q:model.potential(x,q,gravity)
        acceleration=lambda x,q:model.acceleration(x,q,gravity)
        p,metadata=generate_library(item['initial_N'],item['population_seed'],potential,spec['common_cmax'],item['component'])
        initial_support=metadata['weights']
        if min(w['ess'] for w in initial_support.values())<5000:
            p,metadata=generate_library(item['max_N'],item['population_seed'],potential,spec['common_cmax'],item['component'])
        if min(w['ess'] for w in metadata['weights'].values())<5000:raise RuntimeError('Insufficient ESS after declared maximum library size')
        metadata['initial_support']=initial_support;metadata['source_specification']=item
        population_path=folder/'population.npz';np.savez_compressed(population_path,**p)
        inner={k.removeprefix('inner_'):v for k,v in p.items() if k.startswith('inner_')}
        basis=tangent_basis(p['ra_deg'],p['dec_deg']);observations={};hashes={};trajectory_controls={}
        print('POPULATION',*key,len(p['q']),flush=True)
        for scenario in SCENARIOS:
            t,angle,n3=cadence(scenario)
            x,v=propagate(p['position'],p['velocity'],p['q'],p['Mtrue_total'],t,acceleration)
            energy=.5*np.sum(v*v,axis=2)+model.potential(x,p['q'][:,None],gravity)
            relative_energy=float(np.max(abs((energy-p['energy'][:,None])/p['energy'][:,None])))
            lz_error=float(np.max(abs(np.cross(x,v)[:,:,2]-p['lz'][:,None])))
            assert relative_energy<1e-7 and lz_error<1e-10,(key,scenario,relative_energy,lz_error)
            xy=project(x-p['position'][:,None,:],basis)[:,:,:2]*p['rM_au'][:,None,None]
            obs=observe_pairs(None,p['distance_pc'],p['pm_err_masyr'],inner,scenario,item['observer_seed'],p['ra_deg'],p['dec_deg'],
                outer_relative_xy_au=xy,mass_fraction=p['mass_fraction'])
            assert all(np.isfinite(z).all() for z in obs.values())
            for baseline in (34,66):
                assert np.allclose(obs[f'pm{baseline}_masyr'],obs[f'outer{baseline}_masyr']+obs[f'inner_bias{baseline}_masyr']+obs[f'noise{baseline}_masyr'],rtol=1e-10,atol=1e-10)
            op=folder/(scenario+'.npz');np.savez_compressed(op,**obs)
            observations[scenario]=op.relative_to(HERE).as_posix();hashes[scenario]=sha(op)
            trajectory_controls[scenario]=dict(max_relative_energy_change=relative_energy,max_absolute_Lz_change=lz_error,
                later_resolution=later_resolution_fraction(p,inner,t),all_states_inside_potential_domain=True,
                no_later_resolution_selection_applied=True)
            print('OBSERVED',*key,scenario,round(time.perf_counter()-start,1),flush=True)
        metadata.update(trajectory_controls=trajectory_controls,elapsed_seconds=time.perf_counter()-start,
            inner_rotation_coordinates='local sky tangent basis',outer_coordinates='synthetic ICRS; external +z',
            selected_catalogue_stationary=False,parent_invariant_restricted=True)
        record=dict(gravity=item['gravity'],role=item['role'],component=item['component'],prior_au=item['prior_au'],
            population_path=population_path.relative_to(HERE).as_posix(),observations=observations,
            sha256=dict(population=sha(population_path),observations=hashes),metadata=metadata)
        write(folder/'manifest.json',record);manifest['records'].append(record);write(manifest_path,manifest)
    assert len(manifest['records'])==22
    manifest['production_authorized']=True;write(manifest_path,manifest)
    print('BANKS COMPLETE',sha(manifest_path),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--freeze-only',action='store_true');run(parser.parse_args().freeze_only)
