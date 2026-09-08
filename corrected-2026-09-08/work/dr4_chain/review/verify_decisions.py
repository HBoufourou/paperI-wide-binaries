"""Independent reconstruction of the complete coupled-library decision experiment.

No production imports. Rebuilds all observable category codes, weighted mixture
templates, catalogue indices/counts, profiled deviances, composite calibrated
probabilities, classifications and descriptive tables from saved inputs.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[_key]='1'
import argparse,csv,hashlib,json,math,time
from pathlib import Path
import numpy as np

GRAVITIES=('Newton','QUMOND');SCENARIOS=('uniform_white','clustered_correlated','sparse_floor')
BETAS=(1.,1.6);PRIORS=(40,80,60,160);FTEST=(0.,.075,.225,.3)
MIXES=tuple((f,e) for f in (0.,.075,.15,.225,.3) for e in ((0.,) if f==0 else (0.,.5,1.)))
THETA=tuple((b,f,e) for b in BETAS for f,e in MIXES)
WK={1.:'weight_beta1',1.6:'weight_beta1p6',1.3:'weight_beta1p3'}
N=1000;R=1000;B=2048


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for data in iter(lambda:f.read(4*1024*1024),b''):h.update(data)
    return h.hexdigest()


def load(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}


def js(path):return json.loads(Path(path).read_text(encoding='utf-8'))


class Checks:
    def __init__(self):self.rows=[];self.elements=0;self.maxima={}
    def yes(self,name,ok,detail=None):
        self.rows.append({'name':name,'passed':bool(ok),'detail':detail})
        if not ok:print(json.dumps({'FAILED':name,'detail':detail}),flush=True)
    def equal(self,name,a,b):
        a=np.asarray(a);b=np.asarray(b);self.elements+=a.size
        self.yes(name,a.shape==b.shape and np.array_equal(a,b),{'elements':int(a.size)})
    def close(self,name,a,b,group='numerical',atol=2e-7,rtol=2e-11):
        a=np.asarray(a);b=np.asarray(b);self.elements+=a.size
        good=a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all()
        err=float(np.max(abs(a-b))) if good and a.size else 0.
        self.maxima[group]=max(self.maxima.get(group,0.),err)
        self.yes(name,good and np.allclose(a,b,atol=atol,rtol=rtol),{'elements':int(a.size),'max_absolute_error':err})


def mahalanobis2(d,c):
    determinant=c[...,0,0]*c[...,1,1]-c[...,0,1]*c[...,1,0]
    if np.any(determinant<=0):raise ValueError('Nonpositive2x2 covariance')
    return (c[...,1,1]*d[...,0]**2-(c[...,0,1]+c[...,1,0])*d[...,0]*d[...,1]+c[...,0,0]*d[...,1]**2)/determinant


def categories(p,o):
    separation=np.linalg.norm(p['rproj_kau'],axis=1);mass=p['Mphot1']+p['Mphot2']
    rm=np.sqrt(1.3271244e20*mass/1.2e-10)/(149597870700.*1000)
    vk=np.sqrt(1.3271244e20/149597870700.)/1000*np.sqrt(mass/(1000*separation))
    factor=4.740470463533349*p['distance_pc']/(1000*vk)
    unit=p['rproj_kau']/separation[:,None];pm=o['pm66_masyr']
    parallel=factor*(unit[:,0]*pm[:,0]+unit[:,1]*pm[:,1])
    perpendicular=factor*(unit[:,0]*pm[:,1]-unit[:,1]*pm[:,0])
    precision=factor*np.sqrt((o['cov66_masyr2'][:,0,0]+o['cov66_masyr2'][:,1,1])/2)
    source=mahalanobis2(o['pm_source66_masyr']-o['pm_source34_masyr'],o['cov_delta_source_masyr2'])
    relative=mahalanobis2(o['pm66_masyr']-o['pm34_masyr'],o['cov_delta_relative_masyr2'])
    parts=[o[k] for k in ('residual_score_source34','residual_score_source66','acceleration_score_source34','acceleration_score_source66')]
    diagnostic=np.max(np.column_stack(parts+[source,relative[:,None]]),axis=1)
    edges=(-np.sqrt(2),-1.,-.5,0.,.5,1.,np.sqrt(2))
    inds=(np.digitize(separation/rm,(.3,1.,3.)),np.digitize(parallel,edges),np.digitize(perpendicular,edges),
          np.digitize(diagnostic,(3.,10.,100.)),(precision>=.05).astype(int))
    return np.ravel_multi_index(inds,(4,8,8,4,2)).astype(np.uint16)


def probability(code,w,smoothing=False):
    total=w.sum();raw=np.bincount(code,weights=w,minlength=B)/total
    raw/=raw.sum();ess=float(total*total/(w@w))
    return ((ess*raw+1/B)/(ess+1) if smoothing else raw),ess


def parts(f,eta,prior=None):
    raw=[('pure',1-f),(f'comp{prior}',f)] if prior is not None else [('pure',1-f),('comp40',f*(1-eta)),('comp80',f*eta)]
    return [(c,w) for c,w in raw if w>0]


def build_templates(banks):
    component={};union=np.zeros((3,B),bool)
    for hi,g in enumerate(GRAVITIES):
        for beta in BETAS:
            for name in ('pure','comp40','comp80'):
                bank=banks[(g,'fitting',name)]
                for s in range(3):
                    component[(s,hi,beta,name)],_=probability(bank['codes'][s],bank['weights'][beta],True)
                    raw,_=probability(bank['codes'][s],bank['weights'][beta]);union[s]|=raw>0
    template=np.empty((3,2,26,B))
    for si in range(3):
        for hi in range(2):
            for ti,(beta,f,eta) in enumerate(THETA):
                template[si,hi,ti]=sum(w*component[(si,hi,beta,c)] for c,w in parts(f,eta))
    return template,union


def scores(counts,template):
    c=counts.astype(float);entropy=np.sum(c*np.log(np.maximum(c,1)),axis=-1)-N*np.log(N)
    cross=np.einsum('srb,shtb->srht',c,np.log(template),optimize=True)
    surface=np.maximum(2*(entropy[:,:,None,None]-cross),0)
    minimum=surface.min(axis=3)
    return {'surface':surface,'D':minimum,'R':minimum-minimum[:,:,::-1],
            'best':surface.argmin(axis=3),'near_best_count':np.sum(surface<=minimum[:,:,:,None]+2,axis=3)}


def draw_counts(banks,case):
    pieces=parts(case['f'],case['eta'],case['prior']);g=case['gravity'];role=case['role'];beta=case['beta']
    weight=np.concatenate([mix*banks[(g,role,c)]['weights'][beta]/banks[(g,role,c)]['weights'][beta].sum() for c,mix in pieces])
    index=np.random.default_rng(case['seed']).choice(len(weight),size=(R,N),p=weight/weight.sum()).astype(np.uint32)
    counts=[]
    offset=np.arange(R,dtype=np.int64)[:,None]*B
    for si in range(3):
        codes=np.concatenate([banks[(g,role,c)]['codes'][si] for c,_ in pieces])
        counts.append(np.bincount((codes[index]+offset).ravel(),minlength=R*B).reshape(R,B))
    return np.asarray(counts),index


def decisions(stats,cal):
    marginal=np.empty((3,R,2,2))
    for si in range(3):
        for hi in range(2):
            for ai,name in enumerate(('D','R')):
                tails=np.array([(R+1-np.searchsorted(cal[si,hi,ti,:,ai],stats[name][si,:,hi],side='left'))/(R+1) for ti in range(26)])
                marginal[si,:,hi,ai]=tails.max(axis=0)
    joint=np.minimum(1,2*marginal.min(axis=3));reject=joint<=.01
    decision=np.where(reject[:,:,0]&~reject[:,:,1],1,np.where(~reject[:,:,0]&reject[:,:,1],0,-1)).astype(np.int8)
    reason=np.where(decision>=0,2,np.where(reject.all(axis=2),1,0)).astype(np.int8)
    return {'p_marginal':marginal,'p_joint':joint,'rejected':reject,'decision':decision,'reason':reason}


def wilson(k,n):
    z=1.959963984540054;rate=k/n;den=n+z*z
    center=(k+z*z/2)/den;half=z/den*np.sqrt(k*(1-rate)+z*z/4)
    return center-half,center+half


def metrics(case,counts,stat,result,union):
    rows=[];paired=[];hi=case['hypothesis_index']
    for si,s in enumerate(SCENARIOS):
        d=result['decision'][si]
        row=dict(case,scenario=s,N_catalogues=R,N_parent=N,
                 correct=int(np.sum(d==hi)),wrong=int(np.sum((d>=0)&(d!=hi))),indeterminate=int(np.sum(d<0)),
                 both_retained=int(np.sum(result['reason'][si]==0)),both_rejected=int(np.sum(result['reason'][si]==1)),
                 true_retained=int(np.sum(~result['rejected'][si,:,hi])),
                 raw_fitting_empty_fraction_mean=float(counts[si][:,~union[si]].sum()/(R*N)),
                 near_best_true_mean=float(stat['near_best_count'][si,:,hi].mean()),
                 near_best_other_mean=float(stat['near_best_count'][si,:,1-hi].mean()))
        for key in ('correct','wrong','indeterminate','true_retained'):
            row[key+'_rate']=row[key]/R;row[key+'_wilson_lo'],row[key+'_wilson_hi']=wilson(row[key],R)
        row['count_target_passed']=row['correct_rate']>=.95 and row['wrong_rate']<=.01
        row['confidence_target_passed']=row['correct_wilson_lo']>=.95 and row['wrong_wilson_hi']<=.01
        rows.append(row)
    for si in (1,2):
        new=result['decision'][si];ref=result['decision'][0];a=new==hi;b=ref==hi
        paired.append(dict(case,scenario=SCENARIOS[si],reference_scenario=SCENARIOS[0],correct_rate_difference=float(np.mean(a.astype(int)-b.astype(int))),
                           correct_only_new=int(np.sum(a&~b)),correct_only_reference=int(np.sum(b&~a)),same_decision=int(np.sum(new==ref))))
    return rows,paired


def compare_csv(path,expected,checks,label):
    observed=list(csv.DictReader(path.open(newline='',encoding='utf-8')))
    checks.yes(label+'/rows',len(observed)==len(expected))
    for i,(a,b) in enumerate(zip(observed,expected)):
        for k,v in b.items():
            name=f'{label}/{i}/{k}'
            if isinstance(v,(bool,np.bool_)):checks.yes(name,a[k]==str(bool(v)))
            elif isinstance(v,(int,float,np.number)):checks.close(name,float(a[k]),v,group='tables',atol=3e-12,rtol=3e-12)
            else:checks.yes(name,a[k]==('' if v is None else str(v)))


def declared_cases():
    cal=[];test=[]
    for hi,g in enumerate(GRAVITIES):
        for ti,(beta,f,eta) in enumerate(THETA):
            bi=BETAS.index(beta);mi=MIXES.index((f,eta))
            cal.append(dict(gravity=g,hypothesis_index=hi,beta=beta,f=f,eta=eta,theta_index=ti,
                            seed=270000001+1000000*hi+10000*bi+100*mi,prior=None,role='calibration'))
        for bi,beta in enumerate((1.,1.6,1.3)):
            for pi,prior in enumerate(PRIORS):
                for fi,f in enumerate(FTEST):
                    if beta==1.3 and (prior not in (60,160) or f not in (0.,.3)):continue
                    test.append(dict(gravity=g,hypothesis_index=hi,beta=beta,f=f,eta=0.,
                        seed=274000001+3000000*hi+10000*bi+1000*pi+100*fi,prior=prior,role='test',
                        family='endpoint' if prior in (40,80) else 'outside_inner_prior',beta_family='endpoint' if beta in BETAS else 'outside_beta_grid'))
    return cal,test


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path(__file__).parents[1])
    parser.add_argument('--out',type=Path,default=Path(__file__).with_name('decision_verification.json'));parser.add_argument('--wait',action='store_true')
    args=parser.parse_args();root=args.root.resolve();folder=root/'decision/results';start=time.time()
    while not (folder/'test/manifest.json').exists():
        if not args.wait:raise RuntimeError('Waiting for completed calibration and test manifests; no verification report written')
        if time.time()-start>7200:raise TimeoutError('Final decision artifacts absent after bounded wait')
        time.sleep(15)
    compute_start=time.time();c=Checks();inputs={}
    def record(path):inputs[str(Path(path).resolve().relative_to(root))]=sha(path)
    bank_manifest=root/'bank_manifest.json';spec=js(folder/'EXECUTION_PROTOCOL.json');record(bank_manifest);record(folder/'EXECUTION_PROTOCOL.json')
    c.yes('source_hash',sha(root/'decision/engine.py')==spec['engine_sha256'])
    c.yes('protocol_hash',sha(root/'decision/PROTOCOL.md')==spec['protocol_sha256'])
    c.yes('bank_manifest_hash',sha(bank_manifest)==spec['bank_manifest_sha256'])
    c.equal('nuisance_declaration',spec['theta'],THETA);c.yes('N_R',spec['N']==N and spec['repeats']==R)
    calcases,testcases=declared_cases();c.yes('calibration_case_declaration',spec['calibration_cases']==calcases)
    c.yes('test_case_declaration',spec['test_cases']==testcases)
    seeds=[r['seed'] for r in calcases+testcases];c.yes('124_decision_seeds',len(seeds)==len(set(seeds))==124)
    manifest=js(bank_manifest);banks={};codes_count=0
    for rec in manifest['records']:
        key=(rec['gravity'],rec['role'],rec['component']);pp=root/rec['population_path'];p=load(pp)
        c.yes('/'.join(key)+'/population_hash',sha(pp)==rec['sha256']['population']);record(pp)
        code=[]
        for si,s in enumerate(SCENARIOS):
            op=root/rec['observations'][s];c.yes('/'.join(key)+'/'+s+'/hash',sha(op)==rec['sha256']['observations'][s]);record(op)
            code.append(categories(p,load(op)));codes_count+=len(code[-1])
        banks[key]={'codes':np.asarray(code),'weights':{b:p[k] for b,k in WK.items()}}
    print(json.dumps({'reconstructed_observable_category_rows':codes_count,'views':len(banks)*3}),flush=True)
    template,union=build_templates(banks);saved_template=load(folder/'templates.npz');record(folder/'templates.npz')
    c.close('templates/probabilities',template,saved_template['probability'],group='templates',atol=2e-14,rtol=2e-12)
    c.equal('templates/raw_union',union,saved_template['raw_union'])
    independently_sorted=np.empty((3,2,26,R,2));calmanifest=js(folder/'calibration/manifest.json');record(folder/'calibration/manifest.json')
    c.yes('calibration_manifest_cases',len(calmanifest['records'])==52)
    for i,case in enumerate(calcases):
        entry=calmanifest['records'][i];c.yes(f'cal/{i}/case',all(entry[k]==v for k,v in case.items()))
        path=folder/'calibration'/entry['path'];c.yes(f'cal/{i}/hash',sha(path)==entry['sha256']);record(path);stored=load(path)
        counts,index=draw_counts(banks,case);c.equal(f'cal/{i}/indices',index,stored['parent_index']);c.equal(f'cal/{i}/counts',counts,stored['counts'])
        stat=scores(counts,template)
        for k in ('D','R'):c.close(f'cal/{i}/'+k,stat[k],stored[k],group='calibration_scores')
        c.equal(f'cal/{i}/best',stat['best'],stored['best'])
        hi=case['hypothesis_index'];ti=case['theta_index']
        for ai,k in enumerate(('D','R')):independently_sorted[:,hi,ti,:,ai]=np.sort(stat[k][:,:,hi],axis=1)
        if i%13==12:print(json.dumps({'calibration_cases_reconstructed':i+1}),flush=True)
    rp=folder/'calibration/reference.npz';record(rp);c.yes('calibration/reference_hash',sha(rp)==calmanifest['reference_sha256'])
    c.close('calibration/sorted_reference',independently_sorted,load(rp)['sorted'],group='calibration_scores')
    testmanifest=js(folder/'test/manifest.json');record(folder/'test/manifest.json');c.yes('test_manifest_cases',len(testmanifest['records'])==72)
    metric_rows=[];paired_rows=[];aggregate=[]
    for i,case in enumerate(testcases):
        entry=testmanifest['records'][i];c.yes(f'test/{i}/case',all(entry[k]==v for k,v in case.items()))
        path=folder/'test'/entry['path'];c.yes(f'test/{i}/hash',sha(path)==entry['sha256']);record(path);stored=load(path)
        counts,index=draw_counts(banks,case);c.equal(f'test/{i}/indices',index,stored['parent_index']);c.equal(f'test/{i}/counts',counts,stored['counts'])
        stat=scores(counts,template)
        for k in ('surface','D','R'):c.close(f'test/{i}/'+k,stat[k],stored[k],group='test_scores')
        for k in ('best','near_best_count'):c.equal(f'test/{i}/'+k,stat[k],stored[k])
        result=decisions(stat,independently_sorted)
        for k in ('p_marginal','p_joint'):c.close(f'test/{i}/'+k,result[k],stored[k],group='calibrated_probabilities',atol=0,rtol=0)
        for k in ('rejected','decision','reason'):c.equal(f'test/{i}/'+k,result[k],stored[k])
        rows,pairs=metrics(case,counts,stat,result,union);metric_rows.extend(rows);paired_rows.extend(pairs)
        aggregate.append({'case_index':i,'gravity':case['gravity'],'beta':case['beta'],'prior':case['prior'],'f':case['f'],
                          'correct':[r['correct'] for r in rows],'wrong':[r['wrong'] for r in rows],
                          'indeterminate':[r['indeterminate'] for r in rows]})
        if i%12==11:print(json.dumps({'test_cases_reconstructed':i+1}),flush=True)
    for filename,digestkey,rows,label in [('metrics.csv','metrics_sha256',metric_rows,'metrics'),('paired_cadences.csv','paired_sha256',paired_rows,'paired')]:
        path=folder/filename;record(path);c.yes(label+'/hash',sha(path)==testmanifest[digestkey]);compare_csv(path,rows,c,label)
    out={'complete':True,'passed':all(r['passed'] for r in c.rows),'checks_count':len(c.rows),'checks':c.rows,
         'numeric_elements_compared':c.elements,'max_absolute_errors':c.maxima,'input_sha256':inputs,
         'verifier_sha256':sha(__file__),'production_engine_sha256':spec['engine_sha256'],'production_protocol_sha256':spec['protocol_sha256'],
         'observable_views_reconstructed':66,'observable_rows_reconstructed':codes_count,
         'calibration_catalogue_views':156000,'test_catalogue_views':216000,'catalogue_indices_reconstructed':124000000,
         'metric_rows_reconstructed':len(metric_rows),'paired_rows_reconstructed':len(paired_rows),'case_counts':aggregate,
         'scope':'No production imports. Rebuilds all66 category views, weighted templates, all124M indices/counts, scores, p-values, decisions and tables. Does not regenerate raw physical/observer libraries or assert real-Gaia validity.',
         'wall_seconds_including_wait':time.time()-start,'compute_seconds':time.time()-compute_start}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('checks','input_sha256','case_counts')}),flush=True)

if __name__=='__main__':main()
