"""Prospective five-variant simulation scoring; old sources remain unchanged."""
from pathlib import Path
import argparse,csv,hashlib,json,time
import numpy as np
import features as ftr

HERE=Path(__file__).resolve().parent
MODES=ftr.MODES;SCENARIOS=('uniform_white','clustered_correlated','sparse_floor')
GRAVITIES=('Newton','QUMOND');BETAS=(1.,1.3,1.6);EXTRA_BETAS=(1.15,1.45)
PRIORS=(40,80,60,160);FGRID=(0.,.075,.15,.225,.3)
MIXTURES=tuple((f,e) for f in FGRID for e in ((0.,) if f==0 else (0.,.5,1.)))
THETA=tuple((b,f,e) for b in BETAS for f,e in MIXTURES)
N=1000;REPEATS=1000;ALPHA=.01;RULES=('baseline','envelope')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
def load_npz(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}
def array_hash(a):
    a=np.asarray(a);a=np.ascontiguousarray(a.astype(a.dtype.newbyteorder('<'),copy=False))
    h=hashlib.sha256(json.dumps(dict(dtype=a.dtype.str,shape=a.shape),sort_keys=True).encode());h.update(a.tobytes());return h.hexdigest()
def weight_key(beta):return 'weight_beta'+str(beta).replace('.0','').replace('.','p')


def calibration_cases():
    for realization in (0,1,2):
        for hi,g in enumerate(GRAVITIES):
            for ti,(beta,f,eta) in enumerate(THETA):
                bi=BETAS.index(beta);mi=MIXTURES.index((f,eta))
                yield dict(gravity=g,hypothesis_index=hi,role='calibration',realization=realization,beta=beta,f=f,eta=eta,
                    prior=None,theta_index=ti,seed=320000001+1000000*realization+100000*hi+1000*bi+10*mi)


def test_cases():
    for realization in (1,2):
        for hi,g in enumerate(GRAVITIES):
            local=0
            for beta in BETAS+EXTRA_BETAS:
                combinations=[(0.,None)]+([(f,p) for f in (.075,.225,.3) for p in PRIORS] if beta in BETAS else [(.3,60),(.3,160)])
                for f,prior in combinations:
                    yield dict(gravity=g,hypothesis_index=hi,role='test',realization=realization,beta=beta,f=f,eta=None,prior=prior,
                        beta_family='fitted_grid' if beta in BETAS else 'outside_beta_grid',
                        companion_family='pure' if f==0 else ('fitted_prior' if prior in (40,80) else 'outside_inner_prior'),
                        local_case_index=local,seed=325000001+1000000*(realization-1)+100000*hi+1000*local)
                    local+=1
            assert local==45


def source_hashes():
    return {name:sha(HERE/name) for name in ('engine.py','features.py')}


def declare():
    path=HERE/'STATISTICAL_DECLARATION.json'
    cases=list(calibration_cases())+list(test_cases());seeds=[c['seed'] for c in cases]
    physical=[base+1000000*r+100000*h+10000*role+100*c for base in (300000001,303000001)
        for r in (0,1) for h in (0,1) for role in (0,1) for c in range(3 if role==0 else 5)]
    assert len(cases)==len(set(seeds))==414 and len(physical)==len(set(physical))==64 and not set(seeds)&set(physical)
    value=dict(protocol_sha256=sha(HERE.parent/'PROTOCOL.md'),source_sha256=source_hashes(),N=N,repeats=REPEATS,
        modes=list(MODES),rules=list(RULES),theta=[list(x) for x in THETA],calibration_cases=list(calibration_cases()),test_cases=list(test_cases()),
        physical_streams=physical,fixture_input_seed=310000010,fixture_unit_seeds=[310000001,310000002,310000003,310000004],
        fixture_case_seed_rule='310100001+1000*global_case_index',
        calibrated_views=702000,test_views=540000,prospective=True,old_results_are_development=True,
        count_hash_definition='SHA256(sorted-key JSON dtype/shape followed by contiguous little-endian array bytes)',
        envelope='max_library(min(1,2*min_statistic(max_theta(rank_p))))',
        conditional='Full-mixture surface minus its marginal-z surface BEFORE nuisance profiling; no n_z-fixed calibration.')
    if path.exists():
        if read(path)!=value:raise ValueError('Refuse changed statistical source/declaration')
    else:write(path,value)
    return value


def check_declaration():
    value=read(HERE/'STATISTICAL_DECLARATION.json')
    if value['source_sha256']!=source_hashes() or value['protocol_sha256']!=sha(HERE.parent/'PROTOCOL.md'):raise ValueError('Frozen statistical source/protocol changed')
    return value


def required_records():
    return {(g,role,r,c) for g in GRAVITIES for role,rs,cs in (
        ('fitting',(0,),('pure','comp40','comp80')),('calibration',(0,1,2),('pure','comp40','comp80')),
        ('test',(1,2),('pure','comp40','comp80','comp60','comp160')))
        for r in rs for c in cs}


def load_raw_banks(manifest_path,fixture=False):
    path=Path(manifest_path).resolve();manifest=read(path)
    if manifest.get('purpose')=='software_fixture':
        if not fixture:raise ValueError('Explicit fixture mode required')
    elif fixture or manifest.get('production_authorized') is not True:raise ValueError('Complete authorized physical manifest required')
    if tuple(manifest['scenarios'])!=SCENARIOS:raise ValueError('Scenario mismatch')
    banks={}
    for rec in manifest['records']:
        key=(rec['gravity'],rec['role'],int(rec['realization']),rec['component'])
        if key in banks:raise ValueError('Duplicate bank key')
        pp=path.parent/rec['population_path']
        if sha(pp)!=rec['sha256']['population']:raise ValueError('Population hash mismatch')
        pop=load_npz(pp);codes={k:[] for k in ('full66','full34')}
        for scenario in SCENARIOS:
            op=path.parent/rec['observations'][scenario]
            if sha(op)!=rec['sha256']['observations'][scenario]:raise ValueError('Observation hash mismatch')
            obs=load_npz(op)
            for mode in codes:codes[mode].append(ftr.measured_codes(pop,obs,mode))
        n=len(codes['full66'][0]);weights={}
        energy=np.asarray(pop['energy'],float)
        if energy.shape!=(n,) or not np.isfinite(energy).all() or np.any(energy>=0):raise ValueError('Invalid population energy for weights')
        base=np.asarray(pop['weight_beta1'],float)
        for beta in BETAS+EXTRA_BETAS:
            w=base*(-energy)**(beta-1)
            if beta in BETAS:
                stored=np.asarray(pop[weight_key(beta)],float)
                if not np.allclose(w,stored,rtol=1e-12,atol=1e-14):raise ValueError('Stored beta-weight identity mismatch')
                w=stored
            if w.shape!=(n,) or not np.isfinite(w).all() or np.any(w<0) or w.sum()<=0:raise ValueError('Invalid importance weight')
            weights[beta]=w
        banks[key]=dict(codes={k:np.asarray(v,dtype=np.uint16) for k,v in codes.items()},weights=weights,N=n,metadata=rec)
    if set(banks)!=required_records():raise ValueError('Exactly44 declared banks required')
    return banks,manifest


def component_probability(codes,weights,bins):
    w=np.asarray(weights,float);w=w/w.sum();ess=1/np.sum(w*w)
    p=np.bincount(codes,weights=w,minlength=bins);p/=p.sum()
    audit=dict(N=len(w),effective_rows=float(ess),max_weight=float(w.max()),nonempty_bins=int((p>0).sum()),smoothing_mass=float(1/(ess+1)))
    return (ess*p+1/bins)/(ess+1),p>0,audit


def mixture_parts(f,eta=None,prior=None):
    if f==0:return [('pure',1.)]
    if prior is not None:return [('pure',1-f),(f'comp{prior}',f)]
    return [(c,w) for c,w in [('pure',1-f),('comp40',f*(1-eta)),('comp80',f*eta)] if w>0]


def make_templates(banks):
    out={};audits=[]
    for mode in ('full66','full34'):
        comp={};raw=np.zeros((3,2048),bool)
        for hi,g in enumerate(GRAVITIES):
            for beta in BETAS:
                for cname in ('pure','comp40','comp80'):
                    bank=banks[g,'fitting',0,cname]
                    for si,s in enumerate(SCENARIOS):
                        p,occupied,a=component_probability(bank['codes'][mode][si],bank['weights'][beta],2048)
                        comp[si,hi,beta,cname]=p;raw[si]|=occupied
                        audits.append(dict(mode=mode,gravity=g,beta=beta,component=cname,scenario=s,**a))
        grid=np.empty((3,2,len(THETA),2048))
        for si in range(3):
            for hi in range(2):
                for ti,(beta,f,eta) in enumerate(THETA):
                    grid[si,hi,ti]=sum(w*comp[si,hi,beta,c] for c,w in mixture_parts(f,eta))
        out[mode]=grid;out['raw_'+mode]=raw
    for mode in ('no_diagnostics66','covariates_only66'):
        out[mode]=ftr.project_histogram(out['full66'],mode)
        out['raw_'+mode]=ftr.project_histogram(out['raw_full66'].astype(np.uint16),mode)>0
    out['raw_conditional66']=out['raw_full66']
    for mode in ('full66','full34','no_diagnostics66','covariates_only66'):
        assert np.all(out[mode]>0) and np.allclose(out[mode].sum(axis=-1),1)
    return out,audits


def surface(counts,probabilities):
    counts=np.asarray(counts)
    if counts.ndim!=3 or counts.shape[0]!=3 or np.any(counts<0) or not np.all(counts.sum(axis=-1)==N):raise ValueError('Invalid catalogue counts')
    values=[]
    for si in range(3):
        c=counts[si].astype(float);lr=np.zeros_like(c);np.log(c/N,out=lr,where=c>0)
        constant=2*np.sum(c*lr,axis=1)
        values.append((-2*c@np.log(probabilities[si].reshape(-1,c.shape[-1])).T+constant[:,None]).reshape(len(c),2,len(THETA)))
    return np.maximum(np.asarray(values),0.)


def profile(counts,template):
    full=surface(counts['full66'],template['full66'])
    cov=surface(counts['covariates_only66'],template['covariates_only66'])
    conditional=full-cov
    if np.min(conditional)<-1e-8:raise ValueError('Materially negative conditional deviance')
    values=np.stack((full,surface(counts['no_diagnostics66'],template['no_diagnostics66']),
        surface(counts['full34'],template['full34']),cov,np.maximum(conditional,0.)))
    D=values.min(axis=-1)
    return dict(surface=values,D=D,R=D-D[...,::-1],best=np.argmin(values,axis=-1).astype(np.uint8),
        near_best_count=np.sum(values<=D[...,None]+2,axis=-1).astype(np.uint8))


def sample(banks,case,repeats):
    parts=mixture_parts(case['f'],case['eta'],case['prior']);g=case['gravity'];role=case['role'];r=case['realization'];beta=case['beta']
    chosen=[banks[g,role,r,c] for c,_ in parts];sizes=[b['N'] for b in chosen];offset=np.r_[0,np.cumsum(sizes)]
    weights=np.concatenate([mass*(bank['weights'][beta]/bank['weights'][beta].sum()) for (_,mass),bank in zip(parts,chosen)])
    weights/=weights.sum();index=np.random.default_rng(case['seed']).choice(len(weights),size=(repeats,N),p=weights).astype(np.uint32)
    counts={}
    for mode in ('full66','full34'):
        codes=np.concatenate([b['codes'][mode] for b in chosen],axis=1)
        c=np.zeros((3,repeats,2048),np.uint16)
        for si in range(3):
            for j in range(repeats):c[si,j]=np.bincount(codes[si,index[j]],minlength=2048)
        counts[mode]=c
    for mode in ('no_diagnostics66','covariates_only66'):counts[mode]=ftr.project_histogram(counts['full66'],mode).astype(np.uint16)
    lineage=dict(seed=case['seed'],components=[c for c,_ in parts],component_mass=[w for _,w in parts],component_offsets=offset.tolist(),
        beta=beta,gravity=g,role=role,realization=r,N=N,repeats=repeats,all_variants_and_scenarios_paired=True)
    return counts,index,lineage


def decisions(joint):
    rejected=joint<=ALPHA;d=np.full(joint.shape[:-1],-1,np.int8)
    d[~rejected[...,0]&rejected[...,1]]=0;d[rejected[...,0]&~rejected[...,1]]=1
    reason=np.zeros(d.shape,np.int8);reason[rejected.all(axis=-1)]=1;reason[d>=0]=2
    return dict(p_joint=joint,rejected=rejected,decision=d,reason=reason)


def calibrated(stats,reference):
    # references [mode,calibration_library,scenario,hypothesis,theta,repeat,statistic]
    nr=reference.shape[-2];repeat=stats['D'].shape[2]
    marginal=np.zeros((5,3,3,repeat,2,2))
    for mi in range(5):
        for bi in range(3):
            for si in range(3):
                for hi in range(2):
                    for ai,name in enumerate(('D','R')):
                        v=stats[name][mi,si,:,hi];p=np.zeros(repeat)
                        for ti in range(len(THETA)):
                            null=reference[mi,bi,si,hi,ti,:,ai]
                            p=np.maximum(p,(1+nr-np.searchsorted(null,v,side='left'))/(nr+1))
                        marginal[mi,bi,si,:,hi,ai]=p
    library=np.minimum(1.,2*marginal.min(axis=-1))
    joint=np.stack((library[:,0],library.max(axis=1)))
    return dict(p_marginal=marginal,p_library=library,**decisions(joint))


def predict_catalogue(population,observations,scenario,template,reference,mode='full66',rule='envelope',declared_scope_confirmed=False):
    """Measured-only conditional API; unsupported data never force a verdict.

    Caller must supply the matching validated templates/references. This scope
    assertion does not validate a physical population or real Gaia catalogue.
    Only the selected mode's measured-field whitelist is requested.
    """
    try:
        if not declared_scope_confirmed:raise ValueError('Declared physical/population/observational scope must be confirmed')
        if scenario not in SCENARIOS or mode not in MODES or rule not in RULES:raise ValueError('Unknown scenario, mode or calibration rule')
        codes=ftr.measured_codes(population,observations,mode)
        if len(codes)!=N:raise ValueError('Only catalogues of1000 parents are calibrated')
        si=SCENARIOS.index(scenario);mi=MODES.index(mode);ri=RULES.index(rule)
        bins=ftr.BINS[mi];one=np.bincount(codes,minlength=bins);counts=np.tile(one,(3,1,1))
        if mode=='conditional66':
            values=surface(counts,template['full66'])-surface(ftr.project_histogram(counts,'covariates_only66'),template['covariates_only66'])
            if values.min()<-1e-8:raise ValueError('Invalid conditional surface')
            values=np.maximum(values,0.)
        else:values=surface(counts,template[mode])
        # Other axes are computation-only placeholders; only the requested
        # mode/scenario is returned or interpreted.
        d=np.zeros((5,3,1,2));r=np.zeros_like(d);d[mi]=values.min(axis=-1);r[mi]=d[mi]-d[mi,...,::-1]
        result=calibrated(dict(D=d,R=r),reference);decision=int(result['decision'][ri,mi,si,0]);reason=int(result['reason'][ri,mi,si,0])
        return dict(status='conditional',decision=GRAVITIES[decision] if decision>=0 else 'indeterminate',
            reason=('both_retained','both_rejected','one_retained')[reason],scenario=scenario,mode=mode,rule=rule,N=N,
            p_values=result['p_joint'][ri,mi,si,0].tolist(),best_nuisance_indices=values[si,0].argmin(axis=-1).tolist(),
            raw_fitting_empty_fraction=float(one[~template['raw_'+mode][si]].sum()/N),sky_ready=False,
            scope='Finite-library conditional model; population validity, exact n_z conditioning and real-Gaia calibration are not established.')
    except (KeyError,ValueError,TypeError,np.linalg.LinAlgError) as error:
        return dict(status='unsupported',decision='indeterminate',reason=str(error),sky_ready=False)


def freeze(manifest_path,out,fixture=False,fixture_repeats=None):
    declaration=check_declaration();out=Path(out)
    if (out/'EXECUTION_PROTOCOL.json').exists():raise ValueError('Refuse overwrite of frozen execution')
    banks,manifest=load_raw_banks(manifest_path,fixture)
    repeats=REPEATS if fixture_repeats is None else int(fixture_repeats)
    if repeats<2 or (repeats!=REPEATS and not fixture):raise ValueError('Repeat override is software-only')
    cases_cal=[dict(x) for x in declaration['calibration_cases']];cases_test=[dict(x) for x in declaration['test_cases']]
    if fixture:
        for i,case in enumerate(cases_cal+cases_test):case['seed']=310100001+1000*i
    out.mkdir(parents=True,exist_ok=True);folder=out/'codes'
    if folder.exists():raise ValueError('Existing incomplete code cache requires inspection')
    folder.mkdir();records=[];support=[]
    for i,(key,bank) in enumerate(sorted(banks.items())):
        p=folder/f'bank{i:02d}.npz';np.savez_compressed(p,**bank['codes'],**{weight_key(b):w for b,w in bank['weights'].items()})
        records.append(dict(key=list(key),path=p.relative_to(out).as_posix(),sha256=sha(p),N=bank['N'],metadata=bank['metadata']))
        for beta,w in bank['weights'].items():
            w=w/w.sum();support.append(dict(key=list(key),beta=beta,N=len(w),effective_rows=float(1/np.sum(w*w)),max_weight=float(w.max())))
    template,audits=make_templates(banks);np.savez_compressed(out/'templates.npz',**template)
    execution=dict(statistical_declaration_sha256=sha(HERE/'STATISTICAL_DECLARATION.json'),source_sha256=source_hashes(),protocol_sha256=sha(HERE.parent/'PROTOCOL.md'),
        bank_manifest_sha256=sha(manifest_path),bank_manifest_path=Path(manifest_path).resolve().as_posix(),purpose=manifest.get('purpose','physics_production'),fixture=fixture,
        repeats=repeats,N=N,calibration_cases=cases_cal,test_cases=cases_test,codes=records,library_support=support,theta=[list(t) for t in THETA])
    write(out/'EXECUTION_PROTOCOL.json',execution)
    write(out/'template_manifest.json',dict(template_sha256=sha(out/'templates.npz'),execution_protocol_sha256=sha(out/'EXECUTION_PROTOCOL.json'),component_support=audits))
    return execution


def check_execution(manifest_path,out):
    check_declaration();out=Path(out);e=read(out/'EXECUTION_PROTOCOL.json')
    if e['statistical_declaration_sha256']!=sha(HERE/'STATISTICAL_DECLARATION.json') or e['source_sha256']!=source_hashes() or e['bank_manifest_sha256']!=sha(manifest_path):raise ValueError('Changed frozen execution source or manifest')
    tm=read(out/'template_manifest.json')
    if tm['execution_protocol_sha256']!=sha(out/'EXECUTION_PROTOCOL.json') or tm['template_sha256']!=sha(out/'templates.npz'):raise ValueError('Template provenance mismatch')
    banks={}
    for row in e['codes']:
        p=out/row['path']
        if sha(p)!=row['sha256']:raise ValueError('Cached category/weight hash mismatch')
        a=load_npz(p);banks[tuple(row['key'])]=dict(codes={m:a[m] for m in ('full66','full34')},weights={b:a[weight_key(b)] for b in BETAS+EXTRA_BETAS},N=row['N'])
    return e,banks,load_npz(out/'templates.npz')


def calibrate(manifest_path,out):
    out=Path(out);e,banks,template=check_execution(manifest_path,out);folder=out/'calibration'
    if folder.exists():raise ValueError('Refuse calibration overwrite')
    folder.mkdir();refs=np.empty((5,3,3,2,len(THETA),e['repeats'],2));records=[]
    for ci,case in enumerate(e['calibration_cases']):
        counts,index,lineage=sample(banks,case,e['repeats']);stats=profile(counts,template)
        bi,hi,ti=case['realization'],case['hypothesis_index'],case['theta_index']
        refs[:,bi,:,hi,ti,:,0]=stats['D'][...,hi];refs[:,bi,:,hi,ti,:,1]=stats['R'][...,hi]
        p=folder/f'case{ci:03d}.npz';np.savez_compressed(p,parent_index=index,**stats)
        records.append(dict(case,path=p.name,sha256=sha(p),lineage=lineage,count_hashes={k:array_hash(v) for k,v in counts.items()}))
        print('CALIBRATED',ci,bi,hi,ti,flush=True)
    np.savez_compressed(folder/'reference.npz',sorted=np.sort(refs,axis=5))
    write(folder/'manifest.json',dict(records=records,reference_sha256=sha(folder/'reference.npz'),execution_protocol_sha256=sha(out/'EXECUTION_PROTOCOL.json'),
        catalogue_views=3*len(records)*e['repeats'],score_views=5*3*len(records)*e['repeats'],purpose=e['purpose']))


def wilson(k,n):
    z=1.959963984540054;p=k/n;den=1+z*z/n;mid=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return float(mid-half),float(mid+half)


def metrics(case,stats,result,counts,template):
    rows=[];hi=case['hypothesis_index']
    for ri,rule in enumerate(RULES):
        for mi,mode in enumerate(MODES):
            for si,s in enumerate(SCENARIOS):
                d=result['decision'][ri,mi,si];n=len(d);c=counts['full66' if mode=='conditional66' else mode]
                row=dict(case,rule=rule,mode=mode,scenario=s,N_catalogues=n,N_parent=N,
                    correct=int(np.sum(d==hi)),wrong=int(np.sum((d>=0)&(d!=hi))),indeterminate=int(np.sum(d<0)),
                    both_retained=int(np.sum(result['reason'][ri,mi,si]==0)),both_rejected=int(np.sum(result['reason'][ri,mi,si]==1)),
                    true_retained=int(np.sum(~result['rejected'][ri,mi,si,:,hi])),
                    raw_fitting_empty_fraction_mean=float(c[si,:,~template['raw_'+mode][si]].sum()/n/N),
                    near_best_true_mean=float(stats['near_best_count'][mi,si,:,hi].mean()),near_best_other_mean=float(stats['near_best_count'][mi,si,:,1-hi].mean()))
                for label in ('correct','wrong','indeterminate','true_retained'):
                    row[label+'_rate']=row[label]/n;row[label+'_wilson_lo'],row[label+'_wilson_hi']=wilson(row[label],n)
                row['classification_count_pass']=row['correct_rate']>=.95 and row['wrong_rate']<=.01
                row['retention_count_pass']=row['true_retained_rate']>=.99
                row['overall_count_pass']=row['classification_count_pass'] and row['retention_count_pass']
                row['classification_confidence_pass']=row['correct_wilson_lo']>=.95 and row['wrong_wilson_hi']<=.01
                row['retention_confidence_pass']=row['true_retained_wilson_lo']>=.99
                row['overall_confidence_pass']=row['classification_confidence_pass'] and row['retention_confidence_pass']
                rows.append(row)
    return rows


def paired(case,result):
    rows=[];hi=case['hypothesis_index']
    comparisons=[('information',rule,'full66',rule,m) for rule in RULES for m in MODES[1:]]
    comparisons += [('calibration','envelope',m,'baseline',m) for m in MODES]
    for kind,rulea,modea,ruleb,modeb in comparisons:
        for si,s in enumerate(SCENARIOS):
            a=result['decision'][RULES.index(rulea),MODES.index(modea),si];b=result['decision'][RULES.index(ruleb),MODES.index(modeb),si]
            row=dict(case,comparison_type=kind,first_rule=rulea,first_mode=modea,second_rule=ruleb,second_mode=modeb,scenario=s,
                N_catalogues=len(a),same_decision=int(np.sum(a==b)),difference_direction='first_minus_second')
            for label,aa,bb in [('correct',a==hi,b==hi),('wrong',(a>=0)&(a!=hi),(b>=0)&(b!=hi)),('indeterminate',a<0,b<0)]:
                row[label+'_first_only']=int(np.sum(aa&~bb));row[label+'_second_only']=int(np.sum(bb&~aa));row[label+'_difference']=float(aa.mean()-bb.mean())
            rows.append(row)
    return rows


def test(manifest_path,out):
    out=Path(out);e,banks,template=check_execution(manifest_path,out);cm=read(out/'calibration/manifest.json')
    if cm['execution_protocol_sha256']!=sha(out/'EXECUTION_PROTOCOL.json') or cm['reference_sha256']!=sha(out/'calibration/reference.npz'):raise ValueError('Calibration provenance mismatch')
    reference=load_npz(out/'calibration/reference.npz')['sorted'];folder=out/'test'
    if folder.exists():raise ValueError('Refuse test overwrite')
    folder.mkdir();records=[];all_metrics=[];all_pairs=[]
    for ci,case in enumerate(e['test_cases']):
        counts,index,lineage=sample(banks,case,e['repeats']);stats=profile(counts,template);result=calibrated(stats,reference)
        p=folder/f'case{ci:03d}.npz';np.savez_compressed(p,parent_index=index,**stats,**result)
        records.append(dict(case,path=p.name,sha256=sha(p),lineage=lineage,count_hashes={k:array_hash(v) for k,v in counts.items()}))
        all_metrics.extend(metrics(case,stats,result,counts,template));all_pairs.extend(paired(case,result))
        print('TESTED',ci,case['realization'],case['gravity'],case['beta'],case['prior'],case['f'],flush=True)
    for name,rows in [('metrics.csv',all_metrics),('paired.csv',all_pairs)]:
        with (out/name).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write(folder/'manifest.json',dict(records=records,execution_protocol_sha256=sha(out/'EXECUTION_PROTOCOL.json'),calibration_manifest_sha256=sha(out/'calibration/manifest.json'),
        catalogue_views=3*len(records)*e['repeats'],score_views=5*3*len(records)*e['repeats'],decision_views=10*3*len(records)*e['repeats'],
        metric_rows=len(all_metrics),paired_rows=len(all_pairs),metrics_sha256=sha(out/'metrics.csv'),paired_sha256=sha(out/'paired.csv'),purpose=e['purpose'],sky_ready=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=('declare','freeze','calibrate','test'))
    parser.add_argument('--manifest',type=Path);parser.add_argument('--out',type=Path,default=HERE/'results')
    parser.add_argument('--software-fixture',action='store_true');parser.add_argument('--fixture-repeats',type=int)
    a=parser.parse_args();start=time.perf_counter()
    if a.stage=='declare':declare()
    else:
        if a.manifest is None:parser.error('--manifest required for execution')
        if a.stage=='freeze':freeze(a.manifest,a.out,a.software_fixture,a.fixture_repeats)
        elif a.stage=='calibrate':calibrate(a.manifest,a.out)
        else:test(a.manifest,a.out)
    print('COMPLETE',a.stage,'seconds',time.perf_counter()-start,flush=True)
