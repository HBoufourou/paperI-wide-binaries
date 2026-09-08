"""Independent final-library verification; no production imports.

Reuses the unchanged independent upstream polynomial/population checker.
Verifies all accepted initial rows and measured-file algebra; does not replay
rejected proposals or regenerate all trajectory/noise/GLS realizations.
"""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[_k]='1'
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
UP=ROOT.parent/'dr4_chain'
OLD_REVIEW=UP/'review'
sys.path.insert(0,str(OLD_REVIEW))
from verify_populations import Checks,check_population,load
from verify_interpolant import IndependentModel
GRAVITIES=('Newton','QUMOND')
SCENARIOS=('uniform_white','clustered_correlated','sparse_floor')
COMPONENTS=('pure','comp40','comp80','comp60','comp160')
WEIGHTS=(('weight_beta1',1.),('weight_beta1p3',1.3),('weight_beta1p6',1.6))


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def write(path,data):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')


def expected_cases():
    rows=[]
    for realization in (1,2):
        for gravity in GRAVITIES:
            for role in ('calibration','test'):
                for component in COMPONENTS[:3] if role=='calibration' else COMPONENTS:
                    offset=(realization-1)*1000000+GRAVITIES.index(gravity)*100000+(role=='test')*10000+COMPONENTS.index(component)*100
                    rows.append(dict(realization=realization,gravity=gravity,role=role,component=component,
                        prior_au=None if component=='pure' else int(component[4:]),population_seed=300000001+offset,
                        observer_seed=303000001+offset,N_initial=20000,N_max=40000,
                        folder=f'banks/r{realization}_{gravity.lower()}_{role}_{component}'))
    return rows


def key(item):return tuple(item[k] for k in ('realization','gravity','role','component'))


def verify_schedule(checks):
    spec=read(ROOT/'GENERATION_PROTOCOL.json');expected=expected_cases()
    checks.yes('schedule/exact32_case_set',{key(x):x for x in spec['cases']}=={key(x):x for x in expected} and len(spec['cases'])==32)
    seeds=[r[name] for r in expected for name in ('population_seed','observer_seed')]
    checks.yes('schedule/64_unique_seeds',len(seeds)==len(set(seeds))==64 and set(seeds)==set(spec['seeds']) and len(spec['seeds'])==64)
    checks.yes('schedule/reserved_physical_namespace',min(seeds)>=300000000 and max(seeds)<305000000)
    oldspec=read(UP/'GENERATION_PROTOCOL.json')
    checks.yes('schedule/disjoint_upstream_streams',not(set(seeds)&set(oldspec['streams'])))
    checks.yes('schedule/disjoint_reserved_old_and_new_stats_ranges',all(s>=300000000 and s<305000000 for s in seeds))
    checks.yes('schedule/scenarios',tuple(spec['scenarios'])==SCENARIOS)
    checks.yes('schedule/protocol_hash',sha(ROOT/'PROTOCOL.md')==spec['protocol_sha256'])
    checks.yes('schedule/producer_hash',sha(ROOT/'generate_fresh_banks.py')==spec['producer_sha256'])
    for name,digest in spec['upstream_source_sha256'].items():checks.yes('schedule/unchanged_upstream/'+name,sha(ROOT/name)==digest)
    checks.yes('schedule/model_hash',sha(UP/'physics/results/interaction_model.npz')==spec['model_sha256'])
    checks.yes('schedule/physical_validation_hash',sha(UP/'physical_validation.json')==spec['physical_validation_sha256'])
    checks.yes('schedule/old_manifest_hash',sha(UP/'bank_manifest.json')==spec['old_manifest_sha256'])
    empirical=ROOT.parent/'reanalysis/data/observation_arrays.npz'
    checks.yes('schedule/empirical_hash',sha(empirical)==spec['empirical_source_sha256'])
    physical=read(UP/'physical_validation.json')
    checks.yes('schedule/physical_scope',physical['passed'] and physical['common_cmax']==spec['common_cmax'] and spec['complete_astrophysical_population'] is False and spec['sky_ready'] is False)
    return spec


def check_observations(path,digest,p,scenario,label,checks):
    checks.yes(label+'/hash',sha(path)==digest)
    o=load(path);n=len(p['q'])
    checks.yes(label+'/all_arrays_finite',all(np.isfinite(v).all() for v in o.values()))
    required_shapes={}
    for b in (34,66):
        for name in ('pm','outer','inner_bias','noise','bias_noise'):required_shapes[f'{name}{b}_masyr']=(n,2)
        required_shapes[f'cov{b}_masyr2']=(n,2,2)
        required_shapes[f'pm_source{b}_masyr']=(n,2,2)
        for name in ('residual_score','acceleration_score'):required_shapes[f'{name}_source{b}']=(n,2)
    required_shapes.update(cov_cross_relative_masyr2=(n,2,2),cov_delta_relative_masyr2=(n,2,2),cov_delta_source_masyr2=(n,2,2,2))
    checks.yes(label+'/measured_shapes',all(k in o and o[k].shape==v for k,v in required_shapes.items()))
    for b in (34,66):
        checks.near(label+f'/PM_decomposition{b}',o[f'pm{b}_masyr'],o[f'outer{b}_masyr']+o[f'inner_bias{b}_masyr']+o[f'noise{b}_masyr'],atol=1e-10,rtol=1e-10)
        checks.near(label+f'/source_difference{b}',o[f'pm{b}_masyr'],o[f'pm_source{b}_masyr'][:,1]-o[f'pm_source{b}_masyr'][:,0],atol=1e-10,rtol=1e-10)
        checks.near(label+f'/bias_decomposition{b}',o[f'bias_noise{b}_masyr'],o[f'inner_bias{b}_masyr']+o[f'noise{b}_masyr'],atol=1e-10,rtol=1e-10)
        for name in ('residual_score_source','acceleration_score_source'):
            checks.yes(label+f'/{name}{b}_nonnegative',np.all(o[f'{name}{b}']>=0))
    for name in ('cov34_masyr2','cov66_masyr2','cov_delta_source_masyr2','cov_delta_relative_masyr2'):
        c=o[name]
        checks.near(label+'/'+name+'_symmetric',c,c.swapaxes(-1,-2),atol=1e-15,rtol=1e-10)
        checks.yes(label+'/'+name+'_positive',np.all(np.linalg.eigvalsh(c)>0))
    expected=o['cov34_masyr2']+o['cov66_masyr2']-o['cov_cross_relative_masyr2']-o['cov_cross_relative_masyr2'].swapaxes(-1,-2)
    dc=o['cov_delta_relative_masyr2'];relative=float(np.max(np.max(abs(dc-expected),axis=(1,2))/np.max(abs(dc),axis=(1,2))))
    checks.yes(label+'/shared_covariance_identity',relative<1e-9,dict(max_relative_error=relative))
    n3,nt=(24,44) if scenario=='sparse_floor' else (34,66)
    t=o['times_year']
    checks.yes(label+'/declared_epoch_count',int(o['n_epoch34'])==n3 and int(o['n_epoch66'])==nt and t.shape==(nt,) and o['scan_angle_rad'].shape==(nt,))
    checks.yes(label+'/nested_time_support',np.all(np.diff(t)>0) and t[0]>=0 and t[-1]<=66/12 and t[n3-1]<34/12 and t[n3]>=34/12)
    return dict(scenario=scenario,rows=n,shared_covariance_relative_error=relative,sha256=digest)


def check_metadata(item,p,summary,spec,checks):
    label='/'.join(map(str,key(item)));meta=item['metadata'];n=len(p['q']);source=meta['source_specification']
    if item['realization']:
        expected=next(x for x in expected_cases() if key(x)==key(item))
        checks.yes(label+'/exact_source_specification',source==expected)
        checks.yes(label+'/origin',item['origin']=='fresh_physical_generation')
        checks.yes(label+'/population_seed',meta['seed']==expected['population_seed'])
    else:
        gi=GRAVITIES.index(item['gravity']);ri=('fitting','calibration','test').index(item['role']);ci=COMPONENTS.index(item['component'])
        checks.yes(label+'/historical_seed',meta['seed']==source['population_seed']==260000001+gi*1000000+ri*10000+ci*100 and source['observer_seed']==263000001+gi*1000000+ri*10000+ci*100)
    checks.yes(label+'/component_activity',np.all(p['inner_active']==(item['component']!='pure')))
    checks.yes(label+'/prior_label',item['prior_au']==(None if item['component']=='pure' else int(item['component'][4:])))
    for name,beta in WEIGHTS:
        expected=summary['weights'][name];stored=meta['weights'][name]
        checks.near(label+'/ESS/'+name,np.array(expected['ess']),np.array(stored['ess']))
        checks.near(label+'/maximum_weight/'+name,np.array(expected['maximum_normalized_weight']),np.array(stored['maximum_normalized_weight']))
        for q in (.1,.3,1.):checks.near(label+'/q_fraction/'+name+'/'+str(q),np.array(expected['q_fractions'][str(q)]),np.array(stored['q_fractions'][str(q)]))
        first=p[name][:20000];ess=first.sum()**2/np.sum(first*first)
        checks.near(label+'/initial_ESS/'+name,np.array(ess),np.array(meta['initial_support'][name]['ess']))
    need=min(v['ess'] for v in meta['initial_support'].values())<5000
    checks.yes(label+'/declared_adaptive_size',n==(40000 if need else 20000) and meta['requested_rows']==n)
    checks.yes(label+'/final_ESS',min(v['ess'] for v in summary['weights'].values())>=5000)
    extra={}
    for beta in (1.15,1.45):
        # Independent reference uses the reconstructed potential/energy already
        # validated in check_population; explicit no predictor consumes these.
        w=p['weight_beta1']*(-p['energy'])**(beta-1);norm=w/w.sum()
        val=dict(ess=float(1/np.sum(norm*norm)),max_weight=float(norm.max()))
        checks.yes(label+'/additional_beta_positive/'+str(beta),np.isfinite(w).all() and np.all(w>0))
        if item['realization']:
            stored=meta['additional_beta_support'][str(beta)]
            for k in val:checks.near(label+'/additional_beta_'+k+'/'+str(beta),np.array(val[k]),np.array(stored[k]))
        extra[str(beta)]=val
    account=meta['proposal_account'];accepted=sum(x['accepted'] for x in account)
    checks.yes(label+'/proposal_account',accepted==meta['retained_after_final_batch'] and accepted-n==meta['final_batch_discarded'] and all(x['accepted']+x['resolved_rejected']+x['projected_final_rejected']==x['invariant_and_projected_accepted']<=x['proposals'] for x in account))
    for scenario in SCENARIOS:
        trace=meta['trajectory_controls'][scenario]
        checks.yes(label+'/reported_trajectory_threshold/'+scenario,trace['max_relative_energy_change']<1e-7 and trace['max_absolute_Lz_change']<1e-10 and trace['all_states_inside_potential_domain'] is True and trace['no_later_resolution_selection_applied'] is True)
    checks.yes(label+'/scope_metadata',meta['selected_catalogue_stationary'] is False and meta['parent_invariant_restricted'] is True and meta['inner_rotation_coordinates']=='local sky tangent basis' and meta['outer_coordinates']=='synthetic ICRS; external +z')
    return extra


def verify_banks(spec,checks):
    manifest=read(ROOT/'bank_manifest.json');records=manifest['records'];old=read(UP/'bank_manifest.json')
    oldrows={(x['gravity'],x['role'],x['component']):x for x in old['records'] if x['role'] in ('fitting','calibration')}
    expected={key(c) for c in expected_cases()}|{(0,*k) for k in oldrows}
    keys=[key(r) for r in records]
    checks.yes('banks/exact44_keys',set(keys)==expected and len(keys)==len(set(keys))==44)
    checks.yes('banks/complete_authorization',manifest['production_authorized'] is True)
    checks.yes('banks/generation_hash',manifest['generation_protocol_sha256']==sha(ROOT/'GENERATION_PROTOCOL.json'))
    checks.yes('banks/scenarios',tuple(manifest['scenarios'])==SCENARIOS)
    evidence=manifest['physical_validation'];checks.yes('banks/physical_provenance',sha(ROOT/evidence['path'])==evidence['sha256']==spec['physical_validation_sha256'] and evidence['status']=='passed_for_declared_domain')
    model=IndependentModel(UP/'physics/results/interaction_model.npz');model.empirical=load(ROOT.parent/'reanalysis/data/observation_arrays.npz')
    result=[]
    for item in records:
        label='/'.join(map(str,key(item)));pp=ROOT/item['population_path']
        checks.yes(label+'/population_hash',sha(pp)==item['sha256']['population'])
        if not item['realization']:
            olditem=oldrows[(item['gravity'],item['role'],item['component'])]
            checks.yes(label+'/historical_identity',item['origin']=='historical_fixed_input' and item['metadata']==olditem['metadata'] and item['sha256']==olditem['sha256'] and pp.resolve()==(UP/olditem['population_path']).resolve() and all((ROOT/item['observations'][s]).resolve()==(UP/olditem['observations'][s]).resolve() for s in SCENARIOS))
        p=load(pp);summary=check_population(p,model,item['gravity'].lower(),label,checks)
        summary['additional_beta_support']=check_metadata(item,p,summary,spec,checks)
        summary['observations']=[check_observations(ROOT/item['observations'][s],item['sha256']['observations'][s],p,s,label+'/'+s,checks) for s in SCENARIOS]
        summary.update(realization=item['realization'],gravity=item['gravity'],role=item['role'],component=item['component'])
        result.append(summary)
        print(json.dumps(dict(status='checking_complete_bank',bank=label,rows=len(p['q']),completed=len(result),expected=44)),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--schedule-only',action='store_true');parser.add_argument('--wait',action='store_true');parser.add_argument('--max-wait-seconds',type=float,default=7200);args=parser.parse_args()
    start=time.time();checks=Checks();spec=verify_schedule(checks)
    proof=dict(source_sha256=sha(__file__),independent_population_source_sha256=sha(OLD_REVIEW/'verify_populations.py'),independent_model_source_sha256=sha(OLD_REVIEW/'verify_interpolant.py'),generation_protocol_sha256=sha(ROOT/'GENERATION_PROTOCOL.json'))
    schedule=dict(**proof,stage='generation_schedule_only',complete_schedule=True,passed=all(r['passed'] for r in checks.rows),checks=checks.rows,production_complete=False)
    write(HERE/'generation_schedule_review.json',schedule)
    if args.schedule_only:
        print(json.dumps({k:v for k,v in schedule.items() if k!='checks'}),flush=True);return
    if not schedule['passed']:raise RuntimeError('Generation schedule/provenance review failed')
    last=None
    while True:
        try:m=read(ROOT/'bank_manifest.json');count=len(m['records']);ready=count==44 and m.get('production_authorized') is True
        except (FileNotFoundError,json.JSONDecodeError):count=0;ready=False
        if ready:break
        status=dict(**proof,stage='all44_banks',complete=False,status='awaiting_complete_authorized_manifest',records_available=count,expected_records=44)
        write(HERE/'physical_banks_verification.json',status)
        if count!=last:print(json.dumps(status),flush=True);last=count
        if not args.wait or time.time()-start>args.max_wait_seconds:return
        time.sleep(15)
    rows=verify_banks(spec,checks);fresh=[r for r in rows if r['realization']]
    out=dict(**proof,stage='all44_banks',complete=True,passed=all(r['passed'] for r in checks.rows),checks_count=len(checks.rows),checks=checks.rows,populations=rows,
        manifest_sha256=sha(ROOT/'bank_manifest.json'),fresh_banks=len(fresh),historical_banks=len(rows)-len(fresh),
        fresh_rows=sum(r['rows'] for r in fresh),total_accepted_rows=sum(r['rows'] for r in rows),fresh_observation_files=3*len(fresh),all_observation_files=3*len(rows),
        scope='All accepted initial phase-space rows, mass/flux/selection/weights, frozen provenance and measured-file algebra. No production imports; rejected proposals and all epoch trajectory/noise/GLS generation not independently replayed.',seconds=time.time()-start)
    write(HERE/'physical_banks_verification.json',out)
    print(json.dumps({k:v for k,v in out.items() if k not in ('checks','populations')}),flush=True)
    if not out['passed']:raise RuntimeError('Final physical-bank review contains failures')


if __name__=='__main__':main()
