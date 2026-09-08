"""New physical realizations using byte-checked, unchanged upstream functions."""
from pathlib import Path
import argparse,hashlib,json,os,sys,time
import numpy as np

HERE=Path(__file__).resolve().parent
UPSTREAM=HERE.parent/'dr4_chain'
sys.path.insert(0,str(UPSTREAM));sys.path.insert(0,str(UPSTREAM/'physics'))
from population import generate_library
from orbits import propagate,tangent_basis,project
from joint_observer import observe_pairs,cadence,SCENARIOS
from interaction_model import InteractionModel
from generate_banks import later_resolution_fraction


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')


def relative(path):return Path(os.path.relpath(path,HERE)).as_posix()


def freeze():
    physical=read(UPSTREAM/'physical_validation.json');old=read(UPSTREAM/'bank_manifest.json')
    assert physical['passed'] and old['production_authorized'] and len(old['records'])==22
    model=InteractionModel.load()
    assert sha(model.path)==physical['scalar_model_sha256']
    for name,digest in physical['source_sha256'].items():
        path=UPSTREAM/'physics'/name if name=='interaction_model.py' else UPSTREAM/name
        assert sha(path)==digest,('Changed validated source',name)
    oldspec=read(UPSTREAM/'GENERATION_PROTOCOL.json')
    for name,digest in oldspec['source_sha256'].items():assert sha(UPSTREAM/name)==digest
    cases=[]
    components=('pure','comp40','comp80','comp60','comp160')
    for ri in range(2):
        for gi,gravity in enumerate(('Newton','QUMOND')):
            for rolei,role in enumerate(('calibration','test')):
                for ci,component in enumerate(components):
                    if role=='calibration' and ci>2:continue
                    offset=1000000*ri+100000*gi+10000*rolei+100*ci
                    cases.append(dict(gravity=gravity,role=role,realization=ri+1,component=component,
                        prior_au=None if component=='pure' else int(component[4:]),population_seed=300000001+offset,
                        observer_seed=303000001+offset,N_initial=20000,N_max=40000,
                        folder=f'banks/r{ri+1}_{gravity.lower()}_{role}_{component}'))
    seeds=[c[k] for c in cases for k in ('population_seed','observer_seed')]
    assert len(cases)==32 and len(seeds)==len(set(seeds))==64 and max(seeds)<305000000
    assert not(set(seeds)&set(oldspec['streams']))
    sources=['population.py','orbits.py','joint_observer.py','generate_banks.py','physics/interaction_model.py',
        'POPULATION_PROTOCOL.md','AMENDMENT_01_PHOTOMETRIC_ANCHOR.md','TRAJECTORY_PROTOCOL.md','OBSERVER_PROTOCOL.md']
    spec=dict(protocol_sha256=sha(HERE/'PROTOCOL.md'),producer_sha256=sha(Path(__file__)),
        upstream_source_sha256={relative(UPSTREAM/s):sha(UPSTREAM/s) for s in sources},
        model_sha256=sha(model.path),physical_validation_sha256=sha(UPSTREAM/'physical_validation.json'),
        old_manifest_sha256=sha(UPSTREAM/'bank_manifest.json'),empirical_source_sha256=oldspec['empirical_source_sha256'],
        common_cmax=physical['common_cmax'],scenarios=list(SCENARIOS),cases=cases,seeds=seeds,
        scope='Two fresh calibration/test realizations conditional on unchanged old fitting libraries and physical model.',
        complete_astrophysical_population=False,sky_ready=False)
    path=HERE/'GENERATION_PROTOCOL.json'
    if path.exists():assert read(path)==spec, 'Frozen generation inputs changed'
    else:write(path,spec)
    return model,spec,old


def historical_records(old):
    records=[]
    for rec in old['records']:
        if rec['role'] not in ('fitting','calibration'):continue
        row=dict(rec);row['realization']=0;row['origin']='historical_fixed_input'
        row['population_path']=relative(UPSTREAM/rec['population_path'])
        row['observations']={s:relative(UPSTREAM/p) for s,p in rec['observations'].items()}
        assert sha(HERE/row['population_path'])==row['sha256']['population']
        for s,p in row['observations'].items():assert sha(HERE/p)==row['sha256']['observations'][s]
        records.append(row)
    assert len(records)==12
    return records


def run(freeze_only=False):
    model,spec,old=freeze()
    if freeze_only:print('FRESH GENERATION FROZEN',sha(HERE/'GENERATION_PROTOCOL.json'),flush=True);return
    path=HERE/'bank_manifest.json'
    if path.exists():manifest=read(path)
    else:
        manifest=dict(schema_version=1,purpose='final_methodological_strengthening',production_authorized=False,
            generation_protocol_sha256=sha(HERE/'GENERATION_PROTOCOL.json'),scenarios=list(SCENARIOS),
            physical_validation=dict(path=relative(UPSTREAM/'physical_validation.json'),sha256=spec['physical_validation_sha256'],status='passed_for_declared_domain'),
            physical_scope=spec['scope'],records=historical_records(old))
        write(path,manifest)
    assert manifest['generation_protocol_sha256']==sha(HERE/'GENERATION_PROTOCOL.json')
    for case in spec['cases']:
        key=tuple(case[k] for k in ('realization','gravity','role','component'))
        previous=next((r for r in manifest['records'] if tuple(r[k] for k in ('realization','gravity','role','component'))==key),None)
        if previous:
            assert sha(HERE/previous['population_path'])==previous['sha256']['population']
            for s,p in previous['observations'].items():assert sha(HERE/p)==previous['sha256']['observations'][s]
            continue
        folder=HERE/case['folder'];folder.mkdir(parents=True,exist_ok=True)
        if list(folder.glob('*.npz')):raise RuntimeError('Incomplete fresh library; preserve and inspect: '+str(folder))
        start=time.perf_counter();gravity=case['gravity'].lower()
        potential=lambda x,q:model.potential(x,q,gravity)
        acceleration=lambda x,q:model.acceleration(x,q,gravity)
        p,metadata=generate_library(case['N_initial'],case['population_seed'],potential,spec['common_cmax'],case['component'])
        initial_support=metadata['weights']
        if min(w['ess'] for w in initial_support.values())<5000:
            p,metadata=generate_library(case['N_max'],case['population_seed'],potential,spec['common_cmax'],case['component'])
        assert min(w['ess'] for w in metadata['weights'].values())>=5000
        metadata['initial_support']=initial_support;metadata['source_specification']=case
        metadata['additional_beta_support']={}
        for beta in (1.15,1.45):
            w=p['weight_beta1']*(-p['energy'])**(beta-1);w=w/w.sum()
            metadata['additional_beta_support'][str(beta)]=dict(ess=float(1/np.sum(w*w)),max_weight=float(w.max()))
        pp=folder/'population.npz';np.savez_compressed(pp,**p)
        inner={k.removeprefix('inner_'):v for k,v in p.items() if k.startswith('inner_')}
        basis=tangent_basis(p['ra_deg'],p['dec_deg']);observations={};hashes={};trajectory={}
        print('FRESH POPULATION',*key,len(p['q']),flush=True)
        for scenario in SCENARIOS:
            t,angle,n3=cadence(scenario)
            x,v=propagate(p['position'],p['velocity'],p['q'],p['Mtrue_total'],t,acceleration)
            energy=.5*np.sum(v*v,axis=2)+model.potential(x,p['q'][:,None],gravity)
            emax=float(np.max(abs((energy-p['energy'][:,None])/p['energy'][:,None])))
            lmax=float(np.max(abs(np.cross(x,v)[:,:,2]-p['lz'][:,None])))
            assert emax<1e-7 and lmax<1e-10,(key,scenario,emax,lmax)
            xy=project(x-p['position'][:,None,:],basis)[:,:,:2]*p['rM_au'][:,None,None]
            obs=observe_pairs(None,p['distance_pc'],p['pm_err_masyr'],inner,scenario,case['observer_seed'],p['ra_deg'],p['dec_deg'],
                outer_relative_xy_au=xy,mass_fraction=p['mass_fraction'])
            assert all(np.isfinite(z).all() for z in obs.values())
            for baseline in (34,66):assert np.allclose(obs[f'pm{baseline}_masyr'],obs[f'outer{baseline}_masyr']+obs[f'inner_bias{baseline}_masyr']+obs[f'noise{baseline}_masyr'],rtol=1e-10,atol=1e-10)
            op=folder/(scenario+'.npz');np.savez_compressed(op,**obs)
            observations[scenario]=relative(op);hashes[scenario]=sha(op)
            trajectory[scenario]=dict(max_relative_energy_change=emax,max_absolute_Lz_change=lmax,
                later_resolution=later_resolution_fraction(p,inner,t),all_states_inside_potential_domain=True,
                no_later_resolution_selection_applied=True)
            print('FRESH OBSERVED',*key,scenario,round(time.perf_counter()-start,1),flush=True)
        metadata.update(trajectory_controls=trajectory,elapsed_seconds=time.perf_counter()-start,
            selected_catalogue_stationary=False,parent_invariant_restricted=True,
            inner_rotation_coordinates='local sky tangent basis',outer_coordinates='synthetic ICRS; external +z')
        record=dict(gravity=case['gravity'],role=case['role'],realization=case['realization'],component=case['component'],
            prior_au=case['prior_au'],population_path=relative(pp),observations=observations,origin='fresh_physical_generation',
            sha256=dict(population=sha(pp),observations=hashes),metadata=metadata)
        write(folder/'manifest.json',record);manifest['records'].append(record);write(path,manifest)
    assert len(manifest['records'])==44
    manifest['production_authorized']=True;write(path,manifest)
    print('FRESH BANKS COMPLETE',sha(path),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--freeze-only',action='store_true');run(parser.parse_args().freeze_only)
