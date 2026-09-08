"""Independent numerical reconstruction of the prospective five-mode experiment.

No import of production statistics, features, observation or physics modules.
Reads raw saved observables; reproduces the declared RNG only for catalogue
indices. Conditional deviances use the within-stratum multinomial formula,
not subtraction of separately profiled or unprofiled production surfaces.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_key]='1'
import argparse,csv,hashlib,json,time,traceback
from pathlib import Path
import numpy as np

G=('Newton','QUMOND');S=('uniform_white','clustered_correlated','sparse_floor')
MODES=('full66','no_diagnostics66','full34','covariates_only66','conditional66')
RULES=('baseline','envelope');BETA=(1.,1.3,1.6);EXTRA=(1.15,1.45)
MIX=tuple((f,e) for f in (0.,.075,.15,.225,.3) for e in ((0.,) if f==0 else (0.,.5,1.)))
THETA=tuple((b,f,e) for b in BETA for f,e in MIX)
N=1000;NB=2048
Z=np.array([2*(k//512)+k%2 for k in range(NB)])
NO_DIAG=np.array([2*(k//8)+k%2 for k in range(NB)])

def js(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def npz(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k] for k in z.files}
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for a in iter(lambda:f.read(4*1024*1024),b''):h.update(a)
    return h.hexdigest()
def wk(b):return {1.:'weight_beta1',1.3:'weight_beta1p3',1.6:'weight_beta1p6',1.15:'weight_beta1p15',1.45:'weight_beta1p45'}[b]
def ah(a):
    a=np.ascontiguousarray(a.astype(a.dtype.newbyteorder('<'),copy=False))
    h=hashlib.sha256(json.dumps({'dtype':a.dtype.str,'shape':a.shape},sort_keys=True).encode())
    h.update(a.tobytes());return h.hexdigest()

class Checks:
    def __init__(self):self.rows=[];self.maxima={};self.elements=0
    def yes(self,name,ok,detail=None):
        self.rows.append(dict(name=name,passed=bool(ok),detail=detail))
        if not ok:print('FAILED',name,detail,flush=True)
    def equal(self,name,a,b):
        a=np.asarray(a);b=np.asarray(b);self.elements+=a.size
        good=a.shape==b.shape and np.array_equal(a,b)
        self.yes(name,good,dict(elements=int(a.size),mismatches=int(np.count_nonzero(a!=b)) if a.shape==b.shape else None))
    def close(self,name,a,b,group='numeric',atol=2e-7,rtol=2e-11):
        a=np.asarray(a);b=np.asarray(b);self.elements+=a.size
        good=a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all()
        err=float(abs(a-b).max()) if good and a.size else 0.
        self.maxima[group]=max(err,self.maxima.get(group,0.))
        self.yes(name,good and np.allclose(a,b,atol=atol,rtol=rtol),dict(elements=int(a.size),max_absolute_error=err))

def cases():
    cal=[];test=[]
    for r in range(3):
        for h,g in enumerate(G):
            for t,(b,f,e) in enumerate(THETA):
                cal.append(dict(gravity=g,hypothesis_index=h,role='calibration',realization=r,beta=b,f=f,eta=e,prior=None,theta_index=t,
                    seed=320000001+1000000*r+100000*h+1000*BETA.index(b)+10*MIX.index((f,e))))
    for r in (1,2):
        for h,g in enumerate(G):
            j=0
            for b in BETA+EXTRA:
                selections=[(0.,None)]
                selections += [(f,p) for f in (.075,.225,.3) for p in (40,80,60,160)] if b in BETA else [(.3,60),(.3,160)]
                for f,p in selections:
                    test.append(dict(gravity=g,hypothesis_index=h,role='test',realization=r,beta=b,f=f,eta=None,prior=p,
                        beta_family='fitted_grid' if b in BETA else 'outside_beta_grid',companion_family='pure' if f==0 else ('fitted_prior' if p in (40,80) else 'outside_inner_prior'),
                        local_case_index=j,seed=325000001+1000000*(r-1)+100000*h+1000*j));j+=1
    return cal,test

def quadratic(d,c):
    det=c[...,0,0]*c[...,1,1]-c[...,0,1]*c[...,1,0]
    if np.any(det<=0):raise ValueError('Nonpositive covariance determinant')
    return (c[...,1,1]*d[...,0]**2+c[...,0,0]*d[...,1]**2-(c[...,0,1]+c[...,1,0])*d[...,0]*d[...,1])/det

def observed_codes(p,o,baseline):
    # Full34 is deliberately called with a strict four-key observations dict.
    # No energy, true mass, gravity, companion label or future observation read.
    sep=np.linalg.norm(p['rproj_kau'],axis=1);mass=p['Mphot1']+p['Mphot2']
    unit=p['rproj_kau']/sep[:,None]
    rm=np.sqrt(1.3271244e20*mass/1.2e-10)/(149597870700.*1000)
    vk=np.sqrt(1.3271244e20/149597870700.)/1000*np.sqrt(mass/(1000*sep))
    fac=4.740470463533349*p['distance_pc']/(1000*vk)
    pm=o[f'pm{baseline}_masyr'];cov=o[f'cov{baseline}_masyr2']
    par=fac*(pm[:,0]*unit[:,0]+pm[:,1]*unit[:,1])
    perp=fac*(pm[:,1]*unit[:,0]-pm[:,0]*unit[:,1])
    precision=fac*np.sqrt((cov[:,0,0]+cov[:,1,1])/2)
    names=[f'{a}_source{b}' for a in ('residual_score','acceleration_score') for b in ((34,) if baseline==34 else (34,66))]
    diag=[o[k] for k in names]
    if baseline==66:
        diag+=[quadratic(o['pm_source66_masyr']-o['pm_source34_masyr'],o['cov_delta_source_masyr2']),
            quadratic(pm-o['pm34_masyr'],o['cov_delta_relative_masyr2'])[:,None]]
    diagnostic=np.max(np.column_stack(diag),axis=1)
    v_edges=(-np.sqrt(2),-1.,-.5,0.,.5,1.,np.sqrt(2))
    cat=(np.digitize(sep/rm,(.3,1.,3.)),np.digitize(par,v_edges),np.digitize(perp,v_edges),
        np.digitize(diagnostic,(3.,10.,100.)),(precision>=.05).astype(int))
    return np.ravel_multi_index(cat,(4,8,8,4,2)).astype(np.uint16)

def project(a,mode):
    # Explicit grouping of bin IDs, independent of production reshape/sum axes.
    mapping=Z if mode=='covariates_only66' else NO_DIAG
    return np.stack([a[...,mapping==i].sum(axis=-1) for i in range(int(mapping.max())+1)],axis=-1)

def parts(c):
    f=c['f'];p=c['prior'];e=c['eta']
    if f==0:return [('pure',1.)]
    if p is not None:return [('pure',1-f),(f'comp{p}',f)]
    return [(k,w) for k,w in [('pure',1-f),('comp40',f*(1-e)),('comp80',f*e)] if w>0]

def reconstruct_banks(root,result,manifest,execution,ch,hashes):
    raw=js(manifest);ch.yes('authorized bank manifest',raw.get('production_authorized') is True or execution['fixture'])
    ch.equal('scenario order',raw['scenarios'],S)
    records={tuple([x['gravity'],x['role'],int(x['realization']),x['component']]):x for x in raw['records']}
    required={(g,role,r,k) for g in G for role,rs,ks in [('fitting',(0,),('pure','comp40','comp80')),('calibration',(0,1,2),('pure','comp40','comp80')),('test',(1,2),('pure','comp40','comp80','comp60','comp160'))] for r in rs for k in ks}
    ch.yes('44 unique required physical records',len(raw['records'])==44 and set(records)==required)
    ch.yes('44 unique cached records',len(execution['codes'])==44 and {tuple(r['key']) for r in execution['codes']}==required)
    banks={};support=[];raw_count=0
    for j,row in enumerate(execution['codes']):
        key=tuple(row['key']);r=records[key];prefix=f'bank{j:02d}'
        ch.yes(prefix+' metadata identity',row['metadata']==r)
        pp=manifest.parent/r['population_path'];ph=sha(pp);hashes[pp.resolve().as_posix()]=ph
        ch.equal(prefix+' population hash',ph,r['sha256']['population']);p=npz(pp)
        cache_path=result/row['path'];hashes[cache_path.resolve().as_posix()]=sha(cache_path)
        ch.equal(prefix+' cache hash',hashes[cache_path.resolve().as_posix()],row['sha256']);cache=npz(cache_path)
        n=len(p['energy']);ch.yes(prefix+' N',n==row['N']);raw_count+=n
        ch.yes(prefix+' finite negative energy',np.isfinite(p['energy']).all() and (p['energy']<0).all())
        weights={};codes={'full66':[],'full34':[]}
        for b in BETA+EXTRA:
            expected=p['weight_beta1']*np.power(-p['energy'],b-1)
            if b in BETA:
                ch.close(prefix+f' beta{b} identity',expected,p[wk(b)],'weights',1e-14,1e-12);expected=p[wk(b)]
            ch.equal(prefix+f' beta{b} cached weights',expected,cache[wk(b)])
            ch.yes(prefix+f' beta{b} nonnegative weights',np.isfinite(expected).all() and np.all(expected>=0) and expected.sum()>0)
            weights[b]=expected;total=expected.sum();wn=expected/total
            support.append(dict(key=list(key),beta=b,N=n,effective_rows=float(1/np.dot(wn,wn)),max_weight=float(wn.max())))
        measured={k:p[k] for k in ('distance_pc','Mphot1','Mphot2','rproj_kau')}
        for si,sc in enumerate(S):
            op=manifest.parent/r['observations'][sc];oh=sha(op);hashes[op.resolve().as_posix()]=oh
            ch.equal(prefix+' '+sc+' observation hash',oh,r['sha256']['observations'][sc]);o=npz(op)
            only34={k:o[k] for k in ('pm34_masyr','cov34_masyr2','residual_score_source34','acceleration_score_source34')}
            for baseline,oo in ((66,o),(34,only34)):
                mode='full'+str(baseline);code=observed_codes(measured,oo,baseline)
                ch.equal(prefix+' '+sc+' '+mode+' raw measured codes',code,cache[mode][si]);codes[mode].append(code)
        banks[key]=dict(codes={k:np.array(v) for k,v in codes.items()},weights=weights,N=n)
        print('REBUILT BANK',j+1,44,flush=True)
    compare_rows(ch,'library support',execution['library_support'],support)
    return banks,raw_count

def templates(banks):
    out={};audits=[]
    for mode in ('full66','full34'):
        comp={};raw=np.zeros((3,NB),bool)
        for h,g in enumerate(G):
            for b in BETA:
                for k in ('pure','comp40','comp80'):
                    bank=banks[g,'fitting',0,k];w=bank['weights'][b];total=w.sum();ess=total**2/np.dot(w,w)
                    for s,sc in enumerate(S):
                        p=np.bincount(bank['codes'][mode][s],weights=w,minlength=NB)/total;p/=p.sum()
                        raw[s]|=p>0;comp[s,h,b,k]=(ess*p+1/NB)/(ess+1)
                        audits.append(dict(mode=mode,gravity=g,beta=b,component=k,scenario=sc,N=len(w),effective_rows=float(ess),max_weight=float(w.max()/total),nonempty_bins=int((p>0).sum()),smoothing_mass=float(1/(ess+1))))
        out[mode]=np.stack([np.stack([np.stack([sum(m*comp[s,h,b,k] for k,m in parts(dict(f=f,eta=e,prior=None))) for b,f,e in THETA]) for h in range(2)]) for s in range(3)])
        out['raw_'+mode]=raw
    for mode in ('no_diagnostics66','covariates_only66'):
        out[mode]=project(out['full66'],mode);out['raw_'+mode]=project(out['raw_full66'].astype(int),mode)>0
    out['raw_conditional66']=out['raw_full66']
    return out,audits

def sample(banks,c,repeats):
    split=parts(c);selected=[banks[c['gravity'],c['role'],c['realization'],name] for name,_ in split]
    w=np.concatenate([m*(b['weights'][c['beta']]/b['weights'][c['beta']].sum()) for (_,m),b in zip(split,selected)])
    w/=w.sum();index=np.random.default_rng(c['seed']).choice(len(w),size=(repeats,N),p=w).astype(np.uint32)
    offset=(np.arange(repeats)*NB)[:,None];counts={}
    for mode in ('full66','full34'):
        joined=np.concatenate([b['codes'][mode] for b in selected],axis=1)
        counts[mode]=np.stack([np.bincount((joined[s,index]+offset).ravel(),minlength=repeats*NB).reshape(repeats,NB).astype(np.uint16) for s in range(3)])
    for mode in ('no_diagnostics66','covariates_only66'):counts[mode]=project(counts['full66'],mode).astype(np.uint16)
    lineage=dict(seed=c['seed'],components=[k for k,_ in split],component_mass=[v for _,v in split],component_offsets=np.r_[0,np.cumsum([b['N'] for b in selected])].tolist(),
        beta=c['beta'],gravity=c['gravity'],role=c['role'],realization=c['realization'],N=N,repeats=repeats,all_variants_and_scenarios_paired=True)
    return counts,index,lineage

def profile(counts,t):
    values=[]
    for mode in MODES:
        if mode=='conditional66':
            c=counts['full66'].astype(float);nz=counts['covariates_only66'].astype(float)
            logs=np.log(t['full66']/t['covariates_only66'][...,Z])
            entropy=np.sum(c*np.log(np.maximum(c,1)),axis=-1)-np.sum(nz*np.log(np.maximum(nz,1)),axis=-1)
        else:
            c=counts[mode].astype(float);logs=np.log(t[mode])
            entropy=np.sum(c*np.log(np.maximum(c,1)),axis=-1)-N*np.log(N)
        # Three ordinary matrix products avoid a slow shared-axis einsum.
        # The independently derived conditional log density stays unchanged.
        cross=np.stack([(c[s]@logs[s].reshape(78,-1).T).reshape(len(c[s]),2,39) for s in range(3)])
        values.append(np.maximum(2*(entropy[...,None,None]-cross),0))
    surface=np.stack(values);d=surface.min(axis=-1)
    return dict(surface=surface,D=d,R=d-d[...,::-1],best=surface.argmin(axis=-1).astype(np.uint8),near_best_count=np.sum(surface<=d[...,None]+2,axis=-1).astype(np.uint8))

def calibrated(stat,ref):
    nr=ref.shape[-2];repeat=stat['D'].shape[2];marginal=np.empty((5,3,3,repeat,2,2))
    for m in range(5):
        for b in range(3):
            for s in range(3):
                for h in range(2):
                    for a,name in enumerate(('D','R')):
                        # Integer tail counts are maximised before converting to p.
                        tails=np.array([nr-np.searchsorted(ref[m,b,s,h,t,:,a],stat[name][m,s,:,h],side='left') for t in range(39)])
                        marginal[m,b,s,:,h,a]=(tails.max(axis=0)+1)/(nr+1)
    library=np.minimum(2*np.minimum(marginal[...,0],marginal[...,1]),1)
    joint=np.stack([library[:,0],np.maximum.reduce(library,axis=1)])
    accepted=joint>.01;decision=np.full(joint.shape[:-1],-1,np.int8)
    decision[accepted[...,0]&~accepted[...,1]]=0;decision[accepted[...,1]&~accepted[...,0]]=1
    reason=np.where(accepted.sum(axis=-1)==0,1,np.where(accepted.sum(axis=-1)==2,0,2)).astype(np.int8)
    return dict(p_marginal=marginal,p_library=library,p_joint=joint,rejected=~accepted,decision=decision,reason=reason)

def wilson(k,n):
    z=1.959963984540054;den=n+z*z;center=(k+z*z/2)/den;half=z/den*np.sqrt(k*(1-k/n)+z*z/4)
    return center-half,center+half

def metrics(c,stats,p,counts,t):
    rows=[];h=c['hypothesis_index']
    for ri,rule in enumerate(RULES):
        for mi,mode in enumerate(MODES):
            for si,sc in enumerate(S):
                d=p['decision'][ri,mi,si];n=len(d);hist=counts['full66' if mode=='conditional66' else mode][si]
                row=dict(c,rule=rule,mode=mode,scenario=sc,N_catalogues=n,N_parent=N,correct=int((d==h).sum()),wrong=int(((d>=0)&(d!=h)).sum()),indeterminate=int((d==-1).sum()),
                    both_retained=int((p['reason'][ri,mi,si]==0).sum()),both_rejected=int((p['reason'][ri,mi,si]==1).sum()),true_retained=int((p['p_joint'][ri,mi,si,:,h]>.01).sum()),
                    raw_fitting_empty_fraction_mean=float(hist[:,~t['raw_'+mode][si]].sum()/n/N),near_best_true_mean=float(stats['near_best_count'][mi,si,:,h].mean()),near_best_other_mean=float(stats['near_best_count'][mi,si,:,1-h].mean()))
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

def paired(c,p):
    comparisons=[('information',r,'full66',r,m) for r in RULES for m in MODES[1:]]+[('calibration','envelope',m,'baseline',m) for m in MODES]
    rows=[];h=c['hypothesis_index']
    for kind,ra,ma,rb,mb in comparisons:
        for si,sc in enumerate(S):
            a=p['decision'][RULES.index(ra),MODES.index(ma),si];b=p['decision'][RULES.index(rb),MODES.index(mb),si]
            row=dict(c,comparison_type=kind,first_rule=ra,first_mode=ma,second_rule=rb,second_mode=mb,scenario=sc,N_catalogues=len(a),same_decision=int((a==b).sum()),difference_direction='first_minus_second')
            for label,aa,bb in [('correct',a==h,b==h),('wrong',(a>=0)&(a!=h),(b>=0)&(b!=h)),('indeterminate',a<0,b<0)]:
                row[label+'_first_only']=int((aa&~bb).sum());row[label+'_second_only']=int((bb&~aa).sum());row[label+'_difference']=float((aa.astype(int)-bb.astype(int)).mean())
            rows.append(row)
    return rows

def envelope_effects(c,p):
    rows=[];h=c['hypothesis_index']
    for m,mode in enumerate(MODES):
        for s,sc in enumerate(S):
            a=p['decision'][1,m,s];b=p['decision'][0,m,s]
            row={k:c[k] for k in ('realization','gravity','beta','f','prior','seed','beta_family','companion_family')}
            row.update(mode=mode,scenario=sc,N_catalogues=len(a),difference_direction='envelope_minus_baseline')
            for label,aa,bb in [('correct',a==h,b==h),('wrong',(a>=0)&(a!=h),(b>=0)&(b!=h)),('indeterminate',a<0,b<0),('true_retained',~p['rejected'][1,m,s,:,h],~p['rejected'][0,m,s,:,h])]:
                row[label+'_gained']=int((aa&~bb).sum());row[label+'_lost']=int((bb&~aa).sum())
                row[label+'_difference']=float(np.mean(aa.astype(int)-bb.astype(int)))
            row['baseline_double_rejected_to_correct']=int(((p['reason'][0,m,s]==1)&(a==h)).sum())
            row['baseline_double_rejected_to_wrong']=int(((p['reason'][0,m,s]==1)&(a>=0)&(a!=h)).sum())
            row['baseline_correct_to_both_retained']=int(((b==h)&(p['reason'][1,m,s]==0)).sum())
            rows.append(row)
    return rows

def verify_descriptive(result,metrics_rows,pair_rows,effects,library_support,component_support,ch,hashes,fixture):
    """Regroup independently reconstructed rows, not production aggregates."""
    artifact_path=result/'report_artifacts.json'
    if not artifact_path.exists():raise FileNotFoundError('Descriptive report_artifacts.json not yet available; complete report verification requires this artifact.')
    artifact=js(artifact_path);hashes[artifact_path.resolve().as_posix()]=sha(artifact_path)
    # CSV writer semantics matter only for identifiers serialized inside JSON.
    data=[{k:('' if v is None else str(v)) for k,v in r.items()} for r in metrics_rows]
    for r in data:r['evaluation_family']='outside_beta_grid' if r['beta_family']=='outside_beta_grid' else ('outside_inner_prior' if r['companion_family']=='outside_inner_prior' else 'fitted_family')
    gates=('classification_count_pass','retention_count_pass','overall_count_pass','classification_confidence_pass','retention_confidence_pass','overall_confidence_pass')
    def summary(keys):
        grouped={}
        for r in data:grouped.setdefault(tuple(r[k] for k in keys),[]).append(r)
        rows=[];worst=[]
        for key,items in grouped.items():
            row=dict(zip(keys,key));row['cells']=len(items)
            for gate in gates:row[gate]=sum(x[gate]=='True' for x in items)
            for label in ('correct','wrong','true_retained'):
                rates=np.array([float(x[label+'_rate']) for x in items]);j=int(np.argmax(rates) if label=='wrong' else np.argmin(rates));r=items[j]
                row['worst_'+label]=float(r[label+'_rate'])
                row['worst_'+label+'_identity']=json.dumps({k:r[k] for k in ('mode','rule','realization','gravity','beta','f','prior','scenario','seed')},sort_keys=True)
                for tail in ('lo','hi'):row['worst_'+label+'_wilson_'+tail]=float(r[label+'_wilson_'+tail])
                worst.append(dict(extreme=label,**r))
            for reason in ('both_retained','both_rejected'):row['maximum_'+reason+'_rate']=max(int(r[reason])/int(r['N_catalogues']) for r in items)
            for value in ('raw_fitting_empty_fraction_mean','near_best_true_mean'):row['maximum_'+value]=max(float(r[value]) for r in items)
            rows.append(row)
        return rows,worst
    methods,worst=summary(('mode','rule'));families,_=summary(('mode','rule','evaluation_family'));realizations,_=summary(('mode','rule','realization'))
    failures=[r for r in data if r['overall_confidence_pass']=='False'];countfails=[r for r in data if r['overall_count_pass']=='False']
    effect_summary=[]
    for mode in MODES:
        selected=[r for r in effects if r['mode']==mode];row=dict(mode=mode,cells=len(selected))
        for label in ('correct','wrong','indeterminate','true_retained'):
            delta=np.array([r[label+'_difference'] for r in selected])
            row[label+'_increased_cells']=int((delta>1e-14).sum());row[label+'_unchanged_cells']=int((abs(delta)<=1e-14).sum());row[label+'_decreased_cells']=int((delta<-1e-14).sum())
            row[label+'_minimum_difference']=float(delta.min());row[label+'_maximum_difference']=float(delta.max())
            row[label+'_paired_gained']=sum(r[label+'_gained'] for r in selected);row[label+'_paired_lost']=sum(r[label+'_lost'] for r in selected)
        row['double_rejected_to_wrong_total']=sum(r['baseline_double_rejected_to_wrong'] for r in selected);effect_summary.append(row)
    information=[];grouped={};keys=('first_rule','first_mode','second_rule','second_mode')
    for r in pair_rows:
        if r['comparison_type']=='information':grouped.setdefault(tuple(r[k] for k in keys),[]).append(r)
    for key,selected in grouped.items():
        row=dict(zip(keys,key));row['cells']=len(selected)
        for label in ('correct','wrong','indeterminate'):
            delta=np.array([r[label+'_difference'] for r in selected])
            row[label+'_increased_cells']=int((delta>1e-14).sum());row[label+'_unchanged_cells']=int((abs(delta)<=1e-14).sum());row[label+'_decreased_cells']=int((delta<-1e-14).sum())
            row[label+'_minimum_difference']=float(delta.min());row[label+'_maximum_difference']=float(delta.max())
            for suffix in ('first_only','second_only'):row[label+'_'+suffix+'_total']=sum(r[label+'_'+suffix] for r in selected)
        information.append(row)
    library_rows=[]
    for r in library_support:
        g,role,realization,component=r['key'];library_rows.append(dict(gravity=g,role=role,realization=realization,component=component,**{k:v for k,v in r.items() if k!='key'}))
    support_summary=dict(minimum_effective_rows=min(r['effective_rows'] for r in library_rows),maximum_parent_weight=max(r['max_weight'] for r in library_rows),library_beta_rows=len(library_rows),fitted_component_rows=len(component_support))
    tables={'method_summary.csv':methods,'family_summary.csv':families,'realization_summary.csv':realizations,'worst_cells.csv':worst,'confidence_failures.csv':failures,'count_failures.csv':countfails,'envelope_effects.csv':effects,'envelope_summary.csv':effect_summary,
        'information_summary.csv':information,'library_support.csv':library_rows,'fitting_component_support.csv':component_support}
    for name,expected in tables.items():
        path=result/name;hashes[path.resolve().as_posix()]=sha(path)
        with path.open(encoding='utf-8',newline='') as f:actual=list(csv.DictReader(f))
        compare_rows(ch,'descriptive '+name,actual,expected)
    for name,expected in [('methods',methods),('by_family',families),('by_realization',realizations),('envelope',effect_summary),('information',information)]:compare_rows(ch,'artifact '+name,artifact[name],expected)
    compare_rows(ch,'artifact support',[artifact['support']],[support_summary])
    ch.yes('descriptive counts and status',artifact['complete'] is True and artifact['physical_cells']==540 and artifact['metric_rows']==5400 and artifact['paired_rows']==7020 and artifact['confidence_failure_rows']==len(failures) and artifact['count_failure_rows']==len(countfails) and artifact['scientific_status']==('software_only' if fixture else 'conditional_development_experiment'))
    ch.equal('descriptive execution provenance',artifact['execution_protocol_sha256'],sha(result/'EXECUTION_PROTOCOL.json'));ch.equal('descriptive test provenance',artifact['test_manifest_sha256'],sha(result/'test/manifest.json'))
    if not fixture:
        writer=result.parent/'report_results.py';hashes[writer.resolve().as_posix()]=sha(writer)
        ch.equal('descriptive writer identity',artifact['summary_source_sha256'],sha(writer))
    for name in ('metrics.csv','paired.csv'):ch.equal('descriptive input hash '+name,artifact['input_sha256'][name],sha(result/name))
    ch.yes('envelope never loses true retention, independently counted',all(r['true_retained_lost']==0 for r in effects))
    return {name:len(rows) for name,rows in tables.items()}

def compare_rows(ch,label,actual,expected):
    ch.yes(label+' row count',len(actual)==len(expected),dict(actual=len(actual),expected=len(expected)))
    for i,(a,b) in enumerate(zip(actual,expected)):
        bad=[];maximum=0.
        if set(a)!=set(b):bad.append('column names')
        for k,v in b.items():
            x=a.get(k)
            if isinstance(v,(float,np.floating)):
                try:err=abs(float(x)-v);maximum=max(maximum,err);ok=np.isfinite(float(x)) and np.isclose(float(x),v,atol=3e-12,rtol=3e-12)
                except (TypeError,ValueError):ok=False
            elif isinstance(v,(bool,np.bool_)):ok=x==bool(v) if isinstance(x,bool) else x==str(bool(v))
            elif v is None:ok=x is None or x==''
            elif isinstance(v,list):ok=x==v
            elif isinstance(v,str):
                try:
                    xf=float(x);vf=float(v);err=abs(xf-vf);maximum=max(maximum,err)
                    ok=np.isfinite(xf) and np.isfinite(vf) and np.isclose(xf,vf,atol=3e-12,rtol=3e-12)
                except (TypeError,ValueError):ok=x==v
            else:ok=str(x)==str(v)
            if not ok:bad.append(k)
        ch.maxima['descriptive_tables']=max(ch.maxima.get('descriptive_tables',0),maximum)
        ch.yes(f'{label} row{i:04d}',not bad,dict(mismatching_fields=bad,max_absolute_error=maximum))

def self_test(ch):
    # Deterministic counterexamples; no scientific random streams consumed.
    p0=np.array([[.81,.09],[.05,.05]]);p1=np.array([[.02,.08],[.18,.72]])
    mixed=.7*p0+.3*p1
    ch.yes('conditioning cannot precede component mixing',np.max(abs(mixed/mixed.sum(1,keepdims=True)-(.7*p0/p0.sum(1,keepdims=True)+.3*p1/p1.sum(1,keepdims=True))))>.1)
    stat={'D':np.full((5,3,1,2),5.),'R':np.full((5,3,1,2),5.)};ref=np.ones((5,3,3,2,39,4,2))
    ref[:,0,...,1]=10;ref[:,1,...,0]=10;ref[:,2,...,1]=10
    p=calibrated(stat,ref);ch.close('envelope order counterexample',p['p_joint'],.4*np.ones_like(p['p_joint']),atol=0,rtol=0)
    wrong=np.minimum(1,2*p['p_marginal'].max(axis=1).min(axis=-1));ch.yes('wrong envelope order distinguished',np.all(wrong==1))
    for label,k in [('correct',950),('retained',990),('wrong',10)]:
        lo,hi=wilson(k,1000);ch.yes('Wilson is not point rate '+label,lo<k/1000<hi)

def verify(root,result,manifest,ch,hashes,fixture):
    decl_path=root/'stats/STATISTICAL_DECLARATION.json';ex_path=result/'EXECUTION_PROTOCOL.json'
    e=js(ex_path);decl=js(decl_path);tm=js(result/'template_manifest.json');calm=js(result/'calibration/manifest.json');testm=js(result/'test/manifest.json')
    for path in (decl_path,ex_path,root/'PROTOCOL.md',manifest,result/'template_manifest.json',result/'templates.npz',result/'calibration/manifest.json',result/'calibration/reference.npz',result/'test/manifest.json',result/'metrics.csv',result/'paired.csv',root/'stats/engine.py',root/'stats/features.py'):
        hashes[path.resolve().as_posix()]=sha(path)
    if not fixture:hashes[(root/'GENERATION_PROTOCOL.json').resolve().as_posix()]=sha(root/'GENERATION_PROTOCOL.json')
    current={name:sha(root/'stats'/name) for name in ('engine.py','features.py')}
    ch.yes('frozen production source identities',current==decl['source_sha256']==e['source_sha256'])
    ch.equal('root protocol hash',sha(root/'PROTOCOL.md'),decl['protocol_sha256']);ch.equal('execution root protocol hash',sha(root/'PROTOCOL.md'),e['protocol_sha256'])
    ch.equal('declaration chain',sha(decl_path),e['statistical_declaration_sha256']);ch.equal('bank manifest chain',sha(manifest),e['bank_manifest_sha256'])
    ch.equal('template execution chain',sha(ex_path),tm['execution_protocol_sha256']);ch.equal('template file chain',sha(result/'templates.npz'),tm['template_sha256'])
    for stage,m in [('calibration',calm),('test',testm)]:ch.equal(stage+' execution chain',sha(ex_path),m['execution_protocol_sha256'])
    ch.equal('test calibration chain',sha(result/'calibration/manifest.json'),testm['calibration_manifest_sha256'])
    ch.equal('calibration reference hash',sha(result/'calibration/reference.npz'),calm['reference_sha256'])
    for name in ('metrics','paired'):ch.equal(name+' hash',sha(result/(name+'.csv')),testm[name+'_sha256'])
    expected_purpose='software_fixture' if fixture else 'final_methodological_strengthening'
    ch.yes('explicit fixture policy',e['fixture']==fixture and e['purpose']==expected_purpose and (fixture or e['repeats']==1000))
    ch.yes('stage purpose identities',calm['purpose']==testm['purpose']==e['purpose']==js(manifest)['purpose'])
    ch.equal('N parents',e['N'],N);ch.equal('mode order',decl['modes'],MODES);ch.equal('rule order',decl['rules'],RULES)
    ch.equal('39 nuisance grid',e['theta'],THETA);ch.equal('declared nuisance grid',decl['theta'],THETA)
    calcases,testcases=cases();ch.yes('234+180 independently derived declaration cases',decl['calibration_cases']==calcases and decl['test_cases']==testcases)
    physical=[base+1000000*r+100000*h+10000*role+100*k for base in (300000001,303000001) for r in (0,1) for h in (0,1) for role in (0,1) for k in range(3 if role==0 else 5)]
    science=[c['seed'] for c in calcases+testcases];fixtures=[310000010,310000001,310000002,310000003,310000004]+[310100001+1000*i for i in range(414)]
    ch.equal('physical 64 seed declaration',decl['physical_streams'],physical)
    if not fixture:
        generation=js(root/'GENERATION_PROTOCOL.json')
        ch.yes('physical64 streams match frozen generation schedule',len(generation['seeds'])==64 and set(generation['seeds'])==set(physical))
    ch.yes('all scientific/physical/fixture seeds disjoint',len(set(science))==414 and len(set(physical))==64 and len(set(fixtures))==419 and len(set(science+physical+fixtures))==897)
    if fixture:
        for i,c in enumerate(calcases+testcases):c['seed']=310100001+1000*i
    ch.yes('execution cases and seed schedule',e['calibration_cases']==calcases and e['test_cases']==testcases)
    ch.yes('complete stage records',len(calm['records'])==234 and len(testm['records'])==180)
    repeats=e['repeats'];ch.equal('calibration view counts',calm['catalogue_views'],702*repeats);ch.equal('test view counts',testm['catalogue_views'],540*repeats)
    ch.equal('calibration score views',calm['score_views'],3510*repeats);ch.equal('test score views',testm['score_views'],2700*repeats);ch.equal('decision views',testm['decision_views'],5400*repeats)
    ch.yes('scientific outputs do not claim sky readiness',testm['sky_ready'] is False)
    banks,nraw=reconstruct_banks(root,result,manifest,e,ch,hashes)
    t,audits=templates(banks);saved=npz(result/'templates.npz')
    ch.equal('template keys',sorted(t),sorted(saved))
    for name,value in t.items():
        if name.startswith('raw_'):ch.equal('template '+name,value,saved[name])
        else:
            ch.close('template '+name,value,saved[name],'templates',1e-14,2e-12)
            ch.close('template normalization '+name,value.sum(-1),np.ones(value.shape[:-1]),'template normalization',1e-14,1e-14)
            ch.yes('template positive '+name,(value>0).all())
    compare_rows(ch,'component support',tm['component_support'],audits)
    reference=np.empty((5,3,3,2,39,repeats,2));metrics_rows=[];pair_rows=[];effects=[]
    for stage,cases_,manifest_ in [('calibration',calcases,calm),('test',testcases,testm)]:
        if stage=='test':
            reference=np.sort(reference,axis=5);saved_ref=npz(result/'calibration/reference.npz')['sorted']
            ch.close('all sorted calibration references',reference,saved_ref,'calibration reference')
        for ci,c in enumerate(cases_):
            rec=manifest_['records'][ci];prefix=f'{stage}{ci:03d}'
            ch.yes(prefix+' case identity',all(rec.get(k)==v for k,v in c.items()))
            ch.equal(prefix+' expected case path',rec['path'],f'case{ci:03d}.npz')
            path=result/stage/rec['path'];hashed=sha(path);hashes[path.resolve().as_posix()]=hashed
            ch.equal(prefix+' file hash',hashed,rec['sha256']);a=npz(path)
            counts,index,lineage=sample(banks,c,repeats);ch.equal(prefix+' every parent index',index,a['parent_index'])
            ch.yes(prefix+' shared component lineage',rec['lineage']==lineage)
            for mode,count in counts.items():
                ch.equal(prefix+' '+mode+' canonical counts hash',ah(count),rec['count_hashes'][mode])
                ch.yes(prefix+' '+mode+' parent totals',np.all(count.sum(-1)==N))
            stats=profile(counts,t)
            for name in ('surface','D','R'):ch.close(prefix+' '+name,stats[name],a[name],name)
            for name in ('best','near_best_count'):ch.equal(prefix+' '+name,stats[name],a[name])
            if stage=='calibration':
                b,h,ti=c['realization'],c['hypothesis_index'],c['theta_index']
                reference[:,b,:,h,ti,:,0]=stats['D'][...,h];reference[:,b,:,h,ti,:,1]=stats['R'][...,h]
            else:
                p=calibrated(stats,reference)
                for name,value in p.items():ch.equal(prefix+' '+name,value,a[name])
                ch.yes(prefix+' envelope never decreases p',np.all(p['p_joint'][1]>=p['p_joint'][0]))
                ch.yes(prefix+' decision partitions',np.all((p['reason']==2)==(p['decision']>=0)))
                metrics_rows+=metrics(c,stats,p,counts,t);pair_rows+=paired(c,p);effects+=envelope_effects(c,p)
            del a,counts,index,stats
            if ci%10==0 or ci==len(cases_)-1:print('VERIFIED',stage,ci+1,len(cases_),flush=True)
    for name,rows in [('metrics',metrics_rows),('paired',pair_rows)]:
        with (result/(name+'.csv')).open(encoding='utf-8',newline='') as f:actual=list(csv.DictReader(f))
        compare_rows(ch,name,actual,rows);ch.equal(name+' manifest count',len(rows),testm['metric_rows' if name=='metrics' else 'paired_rows'])
    ch.yes('540 unique physical/scenario cells and5400 metric rows',len(metrics_rows)==5400 and len({(r['seed'],r['scenario']) for r in metrics_rows})==540)
    ch.yes('pure cases counted once per beta/gravity/realization',sum(c['f']==0 for c in testcases)==20 and all(c['prior'] is None for c in testcases if c['f']==0))
    descriptive=verify_descriptive(result,metrics_rows,pair_rows,effects,e['library_support'],audits,ch,hashes,fixture)
    return dict(raw_parent_rows=nraw,raw_bank_views=44*3,full_code_elements=nraw*3*2,calibration_cases=234,test_cases=180,catalogue_indices=414*repeats*N,
        calibration_catalogue_views=702*repeats,test_catalogue_views=540*repeats,decision_entries=5400*repeats,metric_rows=len(metrics_rows),paired_rows=len(pair_rows),
        descriptive_tables=descriptive,
        overall_count_pass=int(sum(r['overall_count_pass'] for r in metrics_rows)),overall_confidence_pass=int(sum(r['overall_confidence_pass'] for r in metrics_rows)))

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--results',type=Path)
    p.add_argument('--manifest',type=Path);p.add_argument('--out','--output',dest='out',type=Path);p.add_argument('--wait',action='store_true');p.add_argument('--software-fixture',action='store_true');p.add_argument('--self-test',action='store_true')
    a=p.parse_args();root=a.root.resolve();result=(a.results or root/'stats/results').resolve();manifest=(a.manifest or root/'bank_manifest.json').resolve()
    out=(a.out or root/'review/statistics_verification.json').resolve();ch=Checks();hashes={};start=time.monotonic();waited=0.;complete=False;failure=None;scope={}
    self_test(ch)
    if a.self_test:
        print(json.dumps(dict(passed=all(r['passed'] for r in ch.rows),checks_count=len(ch.rows),checks=ch.rows),indent=2));return
    while not (result/'test/manifest.json').exists():
        if not a.wait or time.monotonic()-start>14400:break
        if waited==0:print('WAITING explicitly for final test manifest; no scientific completion claimed.',flush=True)
        time.sleep(15);waited=time.monotonic()-start
    try:
        if not (result/'test/manifest.json').exists():raise FileNotFoundError('Final test manifest not available; full verification has not run.')
        scope=verify(root,result,manifest,ch,hashes,a.software_fixture);complete=True
    except Exception as e:
        failure=repr(e);traceback.print_exc();ch.yes('verification completed without exception',False,failure)
    changed=[]
    for path,h in hashes.items():
        if sha(path)!=h:changed.append(path)
    ch.yes('all verified input hashes unchanged at completion',not changed,changed)
    doc=dict(passed=complete and all(r['passed'] for r in ch.rows),complete=complete,software_fixture=a.software_fixture,checks_count=len(ch.rows),checks=ch.rows,
        numeric_elements_compared=ch.elements,max_absolute_errors=ch.maxima,scope=scope,input_sha256=hashes,verifier_sha256=sha(__file__),failure=failure,
        independence='No production imports. Reconstructed raw full34/full66 codes, independent2x2 inverse, grouping projections, weighted templates, all RNG indices/counts, direct within-stratum conditional deviance before profiling, integer rank tails, within-library combination then envelope, decisions and descriptive tables.',
        limits=['No physics or astrometric simulation is regenerated by this verifier; raw saved banks are inputs.','RNG reuses the documented NumPy generator/choice algorithm to verify exact stream identity.','Numerical verification is not scientific adequacy or a real-Gaia gravity verdict.','Finite fitting/calibration/test libraries remain conditional simulator laws; marginal Wilson intervals are not simultaneous guarantees.'],
        waited_seconds=waited,compute_seconds=time.monotonic()-start-waited)
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(doc,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:doc[k] for k in ('passed','complete','checks_count','max_absolute_errors','scope','failure','compute_seconds')},indent=2),flush=True)
    if not doc['passed']:raise SystemExit(1)

if __name__=='__main__':main()
