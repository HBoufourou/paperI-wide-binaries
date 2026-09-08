"""Frozen-bin conditional simulation decision; contains no orbit generator."""
from pathlib import Path
import argparse,hashlib,json,math,time
import numpy as np

HERE=Path(__file__).resolve().parent
SCENARIOS=('uniform_white','clustered_correlated','sparse_floor')
GRAVITIES=('Newton','QUMOND')
BETAS=(1.,1.6)
PRIORS=(40,80,60,160)
FTEST=(0.,.075,.225,.3)
FGRID=(0.,.075,.15,.225,.3)
ETAGRID=(0.,.5,1.)
MIXTURES=tuple((f,eta) for f in FGRID for eta in ((0.,) if f==0 else ETAGRID))
THETA=tuple((beta,f,eta) for beta in BETAS for f,eta in MIXTURES)
N=1000
REPEATS=1000
NBIN=2048
ALPHA=.01
WEIGHT_KEYS={1.:'weight_beta1',1.6:'weight_beta1p6',1.3:'weight_beta1p3'}
GM_SUN=1.3271244e20;AU=149597870700.;A0=1.2e-10
K=4.740470463533349
VC=np.sqrt(GM_SUN/AU)/1000
SCALAR=('distance_pc','Mphot1','Mphot2')
SOURCE_SCORES=('residual_score_source34','residual_score_source66','acceleration_score_source34','acceleration_score_source66')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
def load_npz(path):
    with np.load(path,allow_pickle=False) as f:return {k:f[k] for k in f.files}

def checked_covariance(cov,shape):
    cov=np.asarray(cov,float)
    if cov.shape!=shape or not np.isfinite(cov).all():raise ValueError('Invalid covariance shape or missing value')
    if not np.allclose(cov,cov.swapaxes(-1,-2),rtol=1e-10,atol=1e-15):raise ValueError('Asymmetric covariance')
    if np.any(np.linalg.eigvalsh(cov)<=0):raise ValueError('Nonpositive covariance')
    return cov

def qform(delta,cov):return np.einsum('...i,...i->...',delta,np.linalg.solve(cov,delta[...,None])[...,0])

def measured_features(population,observations):
    """Use only measured fields; ignore all additional/latent dictionary entries."""
    p={k:np.asarray(population[k],float) for k in SCALAR}
    if p['distance_pc'].ndim!=1:raise ValueError('Distance must be one dimensional')
    n=len(p['distance_pc'])
    if n==0 or any(v.shape!=(n,) or np.any(v<=0) or not np.isfinite(v).all() for v in p.values()):raise ValueError('Invalid measured covariates')
    r=np.asarray(population['rproj_kau'],float)
    if r.shape!=(n,2) or not np.isfinite(r).all():raise ValueError('rproj_kau must have shape (N,2)')
    s=np.linalg.norm(r,axis=1)
    if np.any(s<=0):raise ValueError('Nonpositive projected separation')
    mu={}
    for k,shape in [('pm34_masyr',(n,2)),('pm66_masyr',(n,2)),('pm_source34_masyr',(n,2,2)),('pm_source66_masyr',(n,2,2))]:
        mu[k]=np.asarray(observations[k],float)
        if mu[k].shape!=shape or not np.isfinite(mu[k]).all():raise ValueError('Invalid proper-motion measurement '+k)
    scores=[]
    for k in SOURCE_SCORES:
        a=np.asarray(observations[k],float)
        if a.shape!=(n,2) or np.any(a<0) or not np.isfinite(a).all():raise ValueError('Invalid source score '+k)
        scores.append(a)
    c66=checked_covariance(observations['cov66_masyr2'],(n,2,2))
    cds=checked_covariance(observations['cov_delta_source_masyr2'],(n,2,2,2))
    cdr=checked_covariance(observations['cov_delta_relative_masyr2'],(n,2,2))
    sd=qform(mu['pm_source66_masyr']-mu['pm_source34_masyr'],cds)
    rd=qform(mu['pm66_masyr']-mu['pm34_masyr'],cdr)
    diagnostic=np.max(np.column_stack(scores+[sd,rd[:,None]]),axis=1)
    m=p['Mphot1']+p['Mphot2'];rm_kau=np.sqrt(GM_SUN*m/A0)/(1000*AU)
    vnewton=VC*np.sqrt(m/(s*1000));factor=K*p['distance_pc']/1000/vnewton
    u=r/s[:,None];velocity=mu['pm66_masyr']*factor[:,None]
    parallel=np.einsum('ni,ni->n',velocity,u)
    perpendicular=velocity[:,1]*u[:,0]-velocity[:,0]*u[:,1]
    precision=np.sqrt(np.trace(c66,axis1=1,axis2=2)/2)*factor
    result=np.column_stack((s/rm_kau,parallel,perpendicular,diagnostic,precision))
    if not np.isfinite(result).all():raise ValueError('Nonfinite derived observable')
    return result

def encode_features(x):
    x=np.asarray(x,float)
    if x.ndim!=2 or x.shape[1]!=5 or not np.isfinite(x).all() or np.any(x[:,[0,3,4]]<0):raise ValueError('Invalid feature matrix')
    sep=np.searchsorted([.3,1.,3.],x[:,0],side='right')
    edges=np.array([-np.sqrt(2),-1.,-.5,0.,.5,1.,np.sqrt(2)])
    par=np.searchsorted(edges,x[:,1],side='right');perp=np.searchsorted(edges,x[:,2],side='right')
    diag=np.searchsorted([3.,10.,100.],x[:,3],side='right')
    noise=(x[:,4]>=.05).astype(int)
    return (((sep*8+par)*8+perp)*4+diag)*2+noise

def probabilities(codes,weights,smooth=False):
    weights=np.asarray(weights,float)
    if weights.shape!=(len(codes),) or np.any(weights<0) or not np.isfinite(weights).all() or weights.sum()<=0:raise ValueError('Invalid library weights')
    w=weights/weights.sum();ess=float(1/np.sum(w*w));p=np.bincount(codes,weights=w,minlength=NBIN)
    if len(p)!=NBIN:raise ValueError('Invalid category code')
    p/=p.sum()
    audit=dict(N=len(w),positive_weight_rows=int((w>0).sum()),effective_rows=ess,max_weight=float(w.max()),
        nonempty_bins=int((p>0).sum()),total_bins=NBIN,smoothing_mass=1/(ess+1) if smooth else 0.)
    if smooth:p=(ess*p+1/NBIN)/(ess+1)
    return p,audit

def mixture_weights(f,eta=0.,prior=None):
    if prior is not None:return [('pure',1-f),(f'comp{prior}',f)]
    return [('pure',1-f),('comp40',f*(1-eta)),('comp80',f*eta)]

def required_records():
    return {(g,r,c) for g in GRAVITIES for r in ('fitting','calibration','test') for c in
        (('pure','comp40','comp80','comp60','comp160') if r=='test' else ('pure','comp40','comp80'))}

def load_banks(manifest_path,allow_software_fixture=False):
    path=Path(manifest_path).resolve();manifest=read_json(path)
    software=manifest.get('purpose')=='software_fixture'
    if software and not allow_software_fixture:raise ValueError('Software-only bank requires explicit fixture mode')
    if not software and not (manifest.get('production_authorized') is True and manifest.get('physical_validation',{}).get('status')=='passed_for_declared_domain'):
        raise ValueError('Physical production has not been authorized for a validated declared domain')
    if not software:
        evidence=manifest['physical_validation'];vp=Path(evidence['path']);vp=vp if vp.is_absolute() else path.parent/vp
        if sha(vp)!=evidence['sha256'] or not manifest.get('physical_scope'):raise ValueError('Missing or mismatched physical-validation provenance')
    if tuple(manifest['scenarios'])!=SCENARIOS:raise ValueError('Scenario contract mismatch')
    banks={};audits=[]
    for rec in manifest['records']:
        key=(rec['gravity'],rec['role'],rec['component'])
        if key in banks:raise ValueError('Duplicate bank '+str(key))
        pp=Path(rec['population_path']);pp=pp if pp.is_absolute() else path.parent/pp
        hashes=rec['sha256']
        if sha(pp)!=hashes['population']:raise ValueError('Population hash mismatch')
        pop=load_npz(pp);allcodes=[]
        for s in SCENARIOS:
            op=Path(rec['observations'][s]);op=op if op.is_absolute() else path.parent/op
            if sha(op)!=hashes['observations'][s]:raise ValueError('Observation hash mismatch')
            observed=load_npz(op);features=measured_features(pop,observed);allcodes.append(encode_features(features))
        weights={beta:np.asarray(pop[wkey],float) for beta,wkey in WEIGHT_KEYS.items() if wkey in pop}
        for beta in BETAS:
            if beta not in weights:raise ValueError('Missing endpoint beta weights')
        if rec['role']=='test' and 1.3 not in weights:raise ValueError('Missing declared outside-grid beta weights')
        banks[key]=dict(codes=np.asarray(allcodes,dtype=np.uint16),weights=weights,N=len(allcodes[0]),metadata=rec)
        for beta,w in weights.items():
            for si,s in enumerate(SCENARIOS):
                _,audit=probabilities(allcodes[si],w)
                audits.append(dict(gravity=key[0],role=key[1],component=key[2],beta=beta,scenario=s,**audit))
    if set(banks)!=required_records():raise ValueError('Need exactly22 declared gravity/role/component banks')
    return banks,manifest,audits

def templates(banks):
    component={};audit=[];raw_union=np.zeros((3,NBIN),bool)
    for hi,g in enumerate(GRAVITIES):
        for beta in BETAS:
            for component_name in ('pure','comp40','comp80'):
                bank=banks[(g,'fitting',component_name)]
                for si,s in enumerate(SCENARIOS):
                    p,a=probabilities(bank['codes'][si],bank['weights'][beta],True)
                    component[(si,hi,beta,component_name)]=p
                    raw,_=probabilities(bank['codes'][si],bank['weights'][beta],False);raw_union[si]|=raw>0
                    audit.append(dict(gravity=g,beta=beta,component=component_name,scenario=s,**a))
    grid=np.empty((3,2,len(THETA),NBIN))
    for si in range(3):
        for hi in range(2):
            for ti,(beta,f,eta) in enumerate(THETA):
                grid[si,hi,ti]=sum(w*component[(si,hi,beta,c)] for c,w in mixture_weights(f,eta))
    assert np.all(grid>0) and np.allclose(grid.sum(axis=-1),1)
    return dict(probability=grid,raw_union=raw_union),audit

def profile(counts,template):
    counts=np.asarray(counts)
    if counts.ndim!=3 or counts.shape[0]!=3 or counts.shape[2]!=NBIN or np.any(counts<0) or not np.all(counts.sum(axis=2)==N):raise ValueError('Counts must be (3,R,2048), each total1000')
    profiles=[]
    for si in range(3):
        c=counts[si].astype(float);logratio=np.zeros_like(c);np.log(c/N,out=logratio,where=c>0)
        constant=2*np.sum(c*logratio,axis=1)
        cross=-2*c@np.log(template['probability'][si].reshape(-1,NBIN)).T
        profiles.append((cross+constant[:,None]).reshape(len(c),2,len(THETA)))
    surface=np.maximum(np.asarray(profiles),0.)
    best=np.argmin(surface,axis=-1);D=surface.min(axis=-1);R=D-D[:,:,::-1]
    return dict(D=D,R=R,best=best,near_best_count=np.sum(surface<=D[:,:,:,None]+2,axis=-1),surface=surface)

def sample_counts(banks,gravity,role,beta,f,eta,seed,prior=None,repeats=REPEATS):
    parts=[(c,w) for c,w in mixture_weights(f,eta,prior) if w>0]
    sizes=[banks[(gravity,role,c)]['N'] for c,_ in parts];offset=np.r_[0,np.cumsum(sizes)]
    flatweights=np.concatenate([mix*(banks[(gravity,role,c)]['weights'][beta]/banks[(gravity,role,c)]['weights'][beta].sum()) for c,mix in parts])
    rng=np.random.default_rng(seed);index=rng.choice(len(flatweights),size=(repeats,N),p=flatweights/flatweights.sum()).astype(np.uint32)
    counts=np.zeros((3,repeats,NBIN),np.uint16)
    for si in range(3):
        codes=np.concatenate([banks[(gravity,role,c)]['codes'][si] for c,_ in parts])
        for i in range(repeats):counts[si,i]=np.bincount(codes[index[i]],minlength=NBIN)
    lineage=dict(seed=seed,components=[c for c,_ in parts],component_mass=[w for _,w in parts],component_offsets=offset.tolist(),
        beta=beta,gravity=gravity,role=role,N=N,repeats=repeats,scenarios_paired=True)
    return counts,index,lineage

def composite_pvalues(stats,calibration):
    # calibration sorted arrays: [scenario,true_hypothesis,nuisance,repeat,statistic].
    nrepeat=calibration.shape[3];out=np.zeros((3,stats['D'].shape[1],2,2))
    for si in range(3):
        for hi in range(2):
            for ai,name in enumerate(('D','R')):
                values=stats[name][si,:,hi]
                p=np.zeros(len(values))
                for ti in range(len(THETA)):
                    null=calibration[si,hi,ti,:,ai]
                    tails=nrepeat-np.searchsorted(null,values,side='left')
                    p=np.maximum(p,(1+tails)/(nrepeat+1))
                out[si,:,hi,ai]=p
    joint=np.minimum(1.,2*np.min(out,axis=-1));rejected=joint<=ALPHA
    decision=np.full(joint.shape[:2],-1,np.int8)
    decision[~rejected[:,:,0]&rejected[:,:,1]]=0
    decision[rejected[:,:,0]&~rejected[:,:,1]]=1
    reason=np.full(joint.shape[:2],0,np.int8) #0 both retained;1both rejected;2classified
    reason[rejected.all(axis=2)]=1;reason[decision>=0]=2
    return dict(p_marginal=out,p_joint=joint,rejected=rejected,decision=decision,reason=reason)

def infer_observed(population,observations,scenario,template,calibration,declared_scope_confirmed=False):
    """Measured-only single-scenario catalogue interface; never force a verdict."""
    try:
        if not declared_scope_confirmed:raise ValueError('Physical/population/observation scope must be explicitly confirmed')
        if scenario not in SCENARIOS:raise ValueError('Unknown observation scenario')
        codes=encode_features(measured_features(population,observations))
        if len(codes)!=N:raise ValueError('Only N=1000 is calibrated')
        si=SCENARIOS.index(scenario)
        # Other repeated rows are computational placeholders, never interpreted.
        count=np.bincount(codes,minlength=NBIN);counts=np.tile(count,(3,1,1))
        stats=profile(counts,template);p=composite_pvalues(stats,calibration);decision=int(p['decision'][si,0]);reason=int(p['reason'][si,0])
        return dict(status='conditional',decision=GRAVITIES[decision] if decision>=0 else 'indeterminate',
            reason=('both_retained','both_rejected','one_retained')[reason],p_values=p['p_joint'][si,0].tolist(),
            best_nuisance_indices=stats['best'][si,0].tolist(),near_best_count=stats['near_best_count'][si,0].tolist(),
            raw_fitting_empty_fraction=float(count[~template['raw_union'][si]].sum()/N),N=N,
            scope='Rank calibration for the calibrating empirical laws; independent-library transfer and physical/real-Gaia validity are not guaranteed.')
    except (KeyError,ValueError,np.linalg.LinAlgError) as error:return dict(status='unsupported',decision='indeterminate',reason=str(error))

def calibration_cases():
    for hi,g in enumerate(GRAVITIES):
        for ti,(beta,f,eta) in enumerate(THETA):
            bi=BETAS.index(beta);mi=MIXTURES.index((f,eta));seed=270000001+1000000*hi+10000*bi+100*mi
            yield dict(gravity=g,hypothesis_index=hi,beta=beta,f=f,eta=eta,theta_index=ti,seed=seed,prior=None,role='calibration')

def test_cases():
    for hi,g in enumerate(GRAVITIES):
        for bi,beta in enumerate((1.,1.6,1.3)):
            for pi,prior in enumerate(PRIORS):
                for fi,f in enumerate(FTEST):
                    if beta==1.3 and (prior not in (60,160) or f not in (0.,.3)):continue
                    seed=274000001+3000000*hi+10000*bi+1000*pi+100*fi
                    yield dict(gravity=g,hypothesis_index=hi,beta=beta,f=f,eta=0.,seed=seed,prior=prior,role='test',
                        family='endpoint' if prior in (40,80) else 'outside_inner_prior',beta_family='endpoint' if beta in BETAS else 'outside_beta_grid')

def freeze_input(manifest_path,out,fixture=False,fixture_repeats=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);banks,manifest,audits=load_banks(manifest_path,fixture)
    repeats=REPEATS if fixture_repeats is None else int(fixture_repeats)
    if repeats!=REPEATS and not (fixture and manifest.get('purpose')=='software_fixture'):raise ValueError('Only software fixtures can change repeats')
    if repeats<2:raise ValueError('At least two fixture repeats required')
    declaration=dict(protocol_sha256=sha(HERE/'PROTOCOL.md'),engine_sha256=sha(__file__),bank_manifest_sha256=sha(manifest_path),
        bank_manifest_path=Path(manifest_path).resolve().as_posix(),purpose=manifest.get('purpose','physics_production'),
        N=N,repeats=repeats,theta=[list(x) for x in THETA],calibration_cases=list(calibration_cases()),test_cases=list(test_cases()),
        library_support=audits,physical_scope=manifest.get('physical_scope'),fixture=fixture)
    if fixture:
        for i,case in enumerate(declaration['calibration_cases']):case['seed']=279100001+1000*i
        for i,case in enumerate(declaration['test_cases']):case['seed']=279300001+1000*i
    seeds=[r['seed'] for r in declaration['calibration_cases']+declaration['test_cases']]
    assert len(seeds)==len(set(seeds))==124
    path=out/'EXECUTION_PROTOCOL.json'
    if path.exists():raise ValueError('Refuse to replace frozen execution protocol')
    write_json(path,declaration);template,audit=templates(banks)
    np.savez_compressed(out/'templates.npz',**template)
    write_json(out/'template_manifest.json',dict(template_sha256=sha(out/'templates.npz'),execution_protocol_sha256=sha(path),component_support=audit))
    return declaration

def check_execution(manifest_path,out,fixture=False):
    out=Path(out);spec=read_json(out/'EXECUTION_PROTOCOL.json')
    if spec['engine_sha256']!=sha(__file__) or spec['protocol_sha256']!=sha(HERE/'PROTOCOL.md') or spec['bank_manifest_sha256']!=sha(manifest_path):raise ValueError('Frozen source/protocol/input mismatch')
    tm=read_json(out/'template_manifest.json')
    if sha(out/'templates.npz')!=tm['template_sha256']:raise ValueError('Template hash mismatch')
    banks,_,_=load_banks(manifest_path,fixture)
    return spec,banks,load_npz(out/'templates.npz')

def run_calibration(manifest_path,out,fixture=False):
    out=Path(out);spec,banks,template=check_execution(manifest_path,out,fixture);folder=out/'calibration'
    if folder.exists():raise ValueError('Refuse to overwrite calibration')
    folder.mkdir();cal=np.empty((3,2,len(THETA),spec['repeats'],2));records=[]
    for case in spec['calibration_cases']:
        counts,index,lineage=sample_counts(banks,case['gravity'],'calibration',case['beta'],case['f'],case['eta'],case['seed'],repeats=spec['repeats'])
        stats=profile(counts,template);hi=case['hypothesis_index'];ti=case['theta_index']
        cal[:,hi,ti,:,0]=stats['D'][:,:,hi];cal[:,hi,ti,:,1]=stats['R'][:,:,hi]
        path=folder/f'h{hi}_theta{ti}.npz';np.savez_compressed(path,counts=counts,parent_index=index,D=stats['D'],R=stats['R'],best=stats['best'])
        records.append(dict(case,path=path.name,sha256=sha(path),lineage=lineage));print('CALIBRATED',hi,ti,flush=True)
    sortedcal=np.sort(cal,axis=3);np.savez_compressed(folder/'reference.npz',sorted=sortedcal)
    write_json(folder/'manifest.json',dict(records=records,reference_sha256=sha(folder/'reference.npz'),execution_protocol_sha256=sha(out/'EXECUTION_PROTOCOL.json'),views=3*len(records)*spec['repeats'],purpose=spec['purpose']))

def wilson(k,n):
    z=1.959963984540054;p=k/n;den=1+z*z/n;m=(p+z*z/(2*n))/den;h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return float(m-h),float(m+h)

def metric_row(case,si,result,stats,counts,template):
    hi=case['hypothesis_index'];d=result['decision'][si];n=len(d);correct=int((d==hi).sum());wrong=int(((d>=0)&(d!=hi)).sum())
    row=dict(case,scenario=SCENARIOS[si],N_catalogues=n,N_parent=N,correct=correct,wrong=wrong,indeterminate=int((d<0).sum()),
        both_retained=int((result['reason'][si]==0).sum()),both_rejected=int((result['reason'][si]==1).sum()),
        true_retained=int((~result['rejected'][si,:,hi]).sum()),raw_fitting_empty_fraction_mean=float(counts[si,:,~template['raw_union'][si]].sum()/n/N),
        near_best_true_mean=float(stats['near_best_count'][si,:,hi].mean()),near_best_other_mean=float(stats['near_best_count'][si,:,1-hi].mean()))
    for key in ('correct','wrong','indeterminate','true_retained'):
        row[key+'_rate']=row[key]/n;row[key+'_wilson_lo'],row[key+'_wilson_hi']=wilson(row[key],n)
    row['count_target_passed']=row['correct_rate']>=.95 and row['wrong_rate']<=.01
    row['confidence_target_passed']=row['correct_wilson_lo']>=.95 and row['wrong_wilson_hi']<=.01
    return row

def run_test(manifest_path,out,fixture=False):
    out=Path(out);spec,banks,template=check_execution(manifest_path,out,fixture);cm=read_json(out/'calibration/manifest.json')
    if sha(out/'calibration/reference.npz')!=cm['reference_sha256']:raise ValueError('Calibration reference hash mismatch')
    cal=load_npz(out/'calibration/reference.npz')['sorted'];folder=out/'test'
    if folder.exists():raise ValueError('Refuse to overwrite test')
    folder.mkdir();records=[];metrics=[];paired=[]
    for ci,case in enumerate(spec['test_cases']):
        counts,index,lineage=sample_counts(banks,case['gravity'],'test',case['beta'],case['f'],0.,case['seed'],prior=case['prior'],repeats=spec['repeats'])
        stats=profile(counts,template);result=composite_pvalues(stats,cal)
        path=folder/f'case{ci:03d}.npz';np.savez_compressed(path,counts=counts,parent_index=index,**stats,**result)
        records.append(dict(case,path=path.name,sha256=sha(path),lineage=lineage))
        for si in range(3):metrics.append(metric_row(case,si,result,stats,counts,template))
        for si in (1,2):
            a=result['decision'][si];b=result['decision'][0];hi=case['hypothesis_index'];ac=a==hi;bc=b==hi
            paired.append(dict(case,scenario=SCENARIOS[si],reference_scenario=SCENARIOS[0],correct_rate_difference=float(ac.mean()-bc.mean()),
                correct_only_new=int((ac&~bc).sum()),correct_only_reference=int((bc&~ac).sum()),same_decision=int((a==b).sum())))
        print('TESTED',ci,case['gravity'],case['beta'],case['prior'],case['f'],flush=True)
    import csv
    for name,rows in [('metrics.csv',metrics),('paired_cadences.csv',paired)]:
        with (out/name).open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write_json(folder/'manifest.json',dict(records=records,execution_protocol_sha256=sha(out/'EXECUTION_PROTOCOL.json'),calibration_manifest_sha256=sha(out/'calibration/manifest.json'),
        views=3*len(records)*spec['repeats'],metrics_sha256=sha(out/'metrics.csv'),paired_sha256=sha(out/'paired_cadences.csv'),
        scientific_status='software_only' if spec['fixture'] else 'conditional_library_experiment',purpose=spec['purpose'],sky_ready=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=('freeze','calibrate','test'));parser.add_argument('--manifest',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path);parser.add_argument('--software-fixture',action='store_true');parser.add_argument('--fixture-repeats',type=int)
    args=parser.parse_args();start=time.perf_counter()
    if args.stage=='freeze':freeze_input(args.manifest,args.out,args.software_fixture,args.fixture_repeats)
    else:
        if args.fixture_repeats is not None:raise ValueError('Fixture repeats are set only when freezing')
        {'calibrate':run_calibration,'test':run_test}[args.stage](args.manifest,args.out,args.software_fixture)
    print('COMPLETE',args.stage,'seconds',time.perf_counter()-start,flush=True)
