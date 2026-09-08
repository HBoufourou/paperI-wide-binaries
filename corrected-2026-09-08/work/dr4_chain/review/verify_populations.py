"""Reconstruct saved phase-space, flux, selection and trajectory checks.

No production module is imported. Uses only the independently evaluated stored
scalar polynomial. This verifies saved rows and traces, not every rejected draw
or the entire original random-number stream.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from verify_interpolant import IndependentModel,sha

GM=1.3271244e20;AU=149597870700.;A0=1.2e-10
WEIGHTS=(('weight_beta1',1.),('weight_beta1p6',1.6),('weight_beta1p3',1.3))


def load(path):
    with np.load(path,allow_pickle=False) as p:return {k:p[k] for k in p.files}


class Checks:
    def __init__(self):self.rows=[]
    def yes(self,name,value,detail=None):
        self.rows.append({'name':name,'passed':bool(value),'detail':detail})
    def near(self,name,a,b,atol=2e-10,rtol=2e-10):
        a=np.asarray(a);b=np.asarray(b)
        ok=a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all()
        error=float(np.max(abs(a-b))) if ok and a.size else 0.
        self.yes(name,ok and np.allclose(a,b,atol=atol,rtol=rtol),{'max_absolute_error':error,'elements':int(a.size)})


def kepler_bisection(mean,e):
    mean=np.mod(mean,2*np.pi);low=np.zeros_like(mean);high=np.full_like(mean,2*np.pi)
    for _ in range(64):
        mid=(low+high)/2;left=mid-e*np.sin(mid)<mean
        low=np.where(left,mid,low);high=np.where(left,high,mid)
    return (low+high)/2


def check_population(p,model,gravity,label,checks):
    n=len(p['q']);x=p['position'];v=p['velocity'];r=np.linalg.norm(x,axis=1)
    if hasattr(model,'empirical'):
        d=model.empirical;idx=p['empirical_row_index']
        mask=(d['d']<200)&(d['R_chance']<.01)&(d['ruwe1']<1.4)&(d['ruwe2']<1.4)
        mask&=(np.maximum(d['sigmax'],d['sigmay'])<.1)&d['has_rv']&np.isfinite(d['sigmaRV'])
        mask&=np.hypot(d['vx_flat'],d['vy_flat'])<=2.23/np.sqrt(d['Mtot'])
        mask&=((d['s_kau']>=.2)&(d['s_kau']<30))|((d['gN_a0']>=.03)&(d['gN_a0']<.3))
        checks.yes(label+'/empirical_row_selection',np.all(mask[idx]))
        checks.near(label+'/empirical_photometric_anchor',p['Mphot_total_anchor'],d['M1'][idx]+d['M2'][idx],atol=0,rtol=0)
        for target,source in [('distance_pc','d'),('ra_deg','ra1'),('dec_deg','dec1')]:
            checks.near(label+'/empirical_'+target,p[target],d[source][idx],atol=0,rtol=0)
        pm=np.stack([np.stack([d['mu1ra_err'][idx],d['mu1dec_err'][idx]],axis=1),
                     np.stack([d['mu2ra_err'][idx],d['mu2dec_err'][idx]],axis=1)],axis=1)
        checks.near(label+'/empirical_PM_errors',p['pm_err_masyr'],pm,atol=0,rtol=0)
    cmax=float(model.raw['cmax'].max())
    psi,_=model.evaluate(x,p['q'],gravity)
    energy=psi+.5*np.sum(v*v,axis=1);lz=x[:,0]*v[:,1]-x[:,1]*v[:,0]
    checks.near(label+'/energy',p['energy'],energy)
    checks.near(label+'/lz',p['lz'],lz)
    discriminant=cmax*cmax+2*energy*lz*lz
    checks.yes(label+'/admissible_energy',np.all((energy<=-.1+1e-12)&(energy>=-5-1e-12)&(energy<-cmax/30)&(discriminant>=0)))
    rlo=lz*lz/(cmax+np.sqrt(np.maximum(discriminant,0)))
    checks.near(label+'/barrier',p['barrier_rlo'],rlo)
    checks.yes(label+'/invariant_radius_domain',np.all((rlo>.035)&(r>=.03)&(r<=30)&(r>=rlo-1e-12)))
    # This barrier is defined for the encoded scalar model and a common Cmax.
    checks.yes(label+'/barrier_at_inner_boundary',np.all(energy<lz*lz/(2*.03**2)-cmax/.03))
    for key,beta in WEIGHTS:
        expected=r**3*(-psi)**3*(-energy)**(beta-1)
        checks.near(label+'/'+key,p[key],expected)
    total=p['Mtrue_total'];rm=np.sqrt(GM*total/A0)/AU
    checks.near(label+'/units_rM',p['rM_au'],rm)
    fractions=np.column_stack((1/(1+p['q']),p['q']/(1+p['q'])))
    checks.near(label+'/source_fractions',p['mass_fraction'],fractions)
    checks.near(label+'/source_masses',p['Mtrue_source'],total[:,None]*fractions)
    ra=np.deg2rad(p['ra_deg']);de=np.deg2rad(p['dec_deg'])
    east=-np.sin(ra)*x[:,0]+np.cos(ra)*x[:,1]
    north=-np.cos(ra)*np.sin(de)*x[:,0]-np.sin(ra)*np.sin(de)*x[:,1]+np.cos(de)*x[:,2]
    bary=np.column_stack((east,north))*rm[:,None]/1000
    checks.near(label+'/projected_barycentres',p['rproj_bary_kau'],bary)
    active=p['inner_active'].astype(bool);ids=np.flatnonzero(active)
    offset=np.zeros((n,2));photo=p['Mtrue_source'].copy()
    if len(ids):
        host=p['inner_host'][ids].astype(int);q=p['inner_q'][ids]
        inner_a=p['inner_a_au'][ids];e=p['inner_e'][ids];rot=p['inner_rotation'][ids]
        checks.yes(label+'/inner_parameter_domain',np.all((q>=.02)&(q<=1)&(e>=0)&(e<1)&(inner_a>=.01)&(inner_a<=300)&np.isin(host,[0,1])))
        checks.yes(label+'/inner_hierarchy',np.all(inner_a*(1+e)<rlo[ids]*rm[ids]/20))
        checks.near(label+'/rotation_orthogonality',np.einsum('nij,nkj->nik',rot,rot),np.tile(np.eye(3),(len(ids),1,1)))
        checks.near(label+'/rotation_determinants',np.linalg.det(rot),np.ones(len(ids)))
        anomaly=kepler_bisection(p['inner_M0'][ids],e)
        plane=np.column_stack((inner_a*(np.cos(anomaly)-e),inner_a*np.sqrt(1-e*e)*np.sin(anomaly)))
        xy=np.einsum('nij,nj->ni',rot[:,:2,:2],plane)
        checks.yes(label+'/unresolved_initial',np.all(np.linalg.norm(xy,axis=1)<p['distance_pc'][ids]))
        coefficient=q**4/(1+q**4)-q/(1+q)
        offset[ids]=(2*host-1)[:,None]*coefficient[:,None]*xy/1000
        source_mass=photo[ids,host]
        checks.near(label+'/inner_monopole_mass',p['inner_mass_host_true'][ids],source_mass)
        photo[ids,host]=source_mass*(1+q**4)**.25/(1+q)
    checks.near(label+'/photometric_masses',np.column_stack((p['Mphot1'],p['Mphot2'])),photo)
    if 'Mphot_total_anchor' in p:
        checks.near(label+'/photometric_anchor',p['Mphot1']+p['Mphot2'],p['Mphot_total_anchor'])
        factor=np.ones(n)
        if len(ids):
            flux=(1+q**4)**.25/(1+q)
            factor[ids]=1-fractions[ids,host]*(1-flux)
        checks.near(label+'/anchored_true_mass',p['Mtrue_total'],p['Mphot_total_anchor']/factor)
        if 'true_mass_flux_factor' in p:checks.near(label+'/true_mass_flux_factor',p['true_mass_flux_factor'],factor)
    checks.near(label+'/photocentre_offset',p['initial_photocentre_offset_kau'],offset)
    checks.near(label+'/observed_projected_separation',p['rproj_kau'],bary+offset)
    observed=np.linalg.norm(p['rproj_kau'],axis=1)
    checks.yes(label+'/projected_selection',np.all((observed>=2)&(observed<=30)))
    checks.yes(label+'/positive_empirical_covariates',np.all(p['distance_pc']>0)&np.all(p['pm_err_masyr']>0))
    return {'label':label,'rows':n,'companion_rows':len(ids),'minimum_rlo':float(rlo.min()),
            'weights':{k:{'ess':float(p[k].sum()**2/np.sum(p[k]**2)),
                          'maximum_normalized_weight':float(p[k].max()/p[k].sum()),
                          'q_fractions':{str(q):float(p[k][p['q']==q].sum()/p[k].sum()) for q in model.q}} for k,_ in WEIGHTS}}


def verify_trajectories(root,model,checks):
    manifest=json.loads((root/'physical_validation.json').read_text());rows=[]
    checks.yes('physical/model_hash',manifest['scalar_model_sha256']==sha(model.path))
    for item in manifest['files']:
        checks.yes('physical/file_hash/'+item['path'],sha(root/item['path'])==item['sha256'])
    recorded={r['name']:r for r in manifest['checks']}
    def compare(label,val):
        checks.near('trajectory/recount/'+label,np.array(float(recorded[label]['value'])),np.array(float(val)),atol=3e-12,rtol=1e-4)
        checks.yes('trajectory/criterion/'+label,val<recorded[label]['threshold'])
    folder=root/'trajectory_validation'
    for gravity in ('newton','qumond'):
        p=load(folder/(gravity+'_initial.npz'))
        rows.append(check_population(p,model,gravity,'trajectory_initial/'+gravity,checks))
        for scenario in ('uniform_white','clustered_correlated','sparse_floor'):
            t=load(folder/(gravity+'_'+scenario+'.npz'));xx=t['positions'];vv=t['velocities']
            energy=model.potential(xx,p['q'][:,None],gravity)+.5*np.sum(vv*vv,axis=2)
            compare(gravity+'/'+scenario+'/energy',np.max(abs((energy-p['energy'][:,None])/p['energy'][:,None])))
            compare(gravity+'/'+scenario+'/Lz',np.max(abs(np.cross(xx,vv)[:,:,2]-p['lz'][:,None])))
            compare(gravity+'/'+scenario+'/position_refinement_AU',np.max(np.linalg.norm(t['positions_refined']-xx,axis=2)*p['rM_au'][:,None]))
            for baseline in (34,66):
                residual=t[f'pm{baseline}_masyr']-t[f'outer{baseline}_masyr']-t[f'inner_bias{baseline}_masyr']-t[f'noise{baseline}_masyr']
                compare(gravity+'/'+scenario+f'/joint_decomposition{baseline}',np.max(abs(residual)))
            cov=t['cov34_masyr2']+t['cov66_masyr2']-t['cov_cross_relative_masyr2']-t['cov_cross_relative_masyr2'].swapaxes(-1,-2)
            delta=t['cov_delta_relative_masyr2']
            compare(gravity+'/'+scenario+'/shared_covariance',np.max(np.max(abs(delta-cov),axis=(1,2))/np.max(abs(delta),axis=(1,2))))
        long=load(folder/(gravity+'_long.npz'));ids=long['ids'];errors=[]
        for suffix in ('coarse','fine'):
            xx=long['positions_'+suffix];vv=long['velocities_'+suffix]
            energy=model.potential(xx,p['q'][ids,None],gravity)+.5*np.sum(vv*vv,axis=2)
            errors.append(float(np.max(abs((energy-p['energy'][ids,None])/p['energy'][ids,None]))))
            checks.yes('trajectory/long_radius/'+gravity+'/'+suffix,np.all((np.linalg.norm(xx,axis=2)>.03)&(np.linalg.norm(xx,axis=2)<30)))
            checks.near('trajectory/long_initial_position/'+gravity+'/'+suffix,xx[:,0],p['position'][ids])
            checks.near('trajectory/long_initial_velocity/'+gravity+'/'+suffix,vv[:,0],p['velocity'][ids])
        compare(gravity+'/long_energy',errors[1])
        compare(gravity+'/long_Lz',np.max(abs(np.cross(long['positions_fine'],long['velocities_fine'])[:,:,2]-p['lz'][ids,None])))
        compare(gravity+'/long_refinement',errors[1]/errors[0] if errors[0]>1e-9 else 0.)
    return rows


def verify_banks(root,model,checks):
    manifest_path=root/'bank_manifest.json';manifest=json.loads(manifest_path.read_text())
    rows=[];keys=[];seeds=[]
    expected={(g,r,c) for g in ('Newton','QUMOND') for r in ('fitting','calibration','test') for c in
              (('pure','comp40','comp80','comp60','comp160') if r=='test' else ('pure','comp40','comp80'))}
    for item in manifest['records']:
        key=(item['gravity'],item['role'],item['component']);keys.append(key);label='/'.join(key)
        path=Path(item['population_path']);path=path if path.is_absolute() else root/path
        checks.yes(label+'/sha256',sha(path)==item['sha256']['population'])
        p=load(path);rows.append(check_population(p,model,item['gravity'].lower(),label,checks))
        meta=item.get('metadata',{})
        rows[-1]['metadata']=meta
        gi=('Newton','QUMOND').index(key[0]);ri=('fitting','calibration','test').index(key[1]);ci=('pure','comp40','comp80','comp60','comp160').index(key[2])
        expected_seed=260000001+gi*1000000+ri*10000+ci*100
        observer_seed=263000001+gi*1000000+ri*10000+ci*100
        source=meta['source_specification']
        checks.yes(label+'/seeds',meta['seed']==expected_seed==source['population_seed'] and source['observer_seed']==observer_seed)
        seeds.extend([expected_seed,observer_seed])
        for k,_ in WEIGHTS:
            w=p[k];norm=w/w.sum();support=meta['weights'][k]
            checks.near(label+'/ESS/'+k,np.array(support['ess']),np.array(1/np.sum(norm*norm)))
            checks.near(label+'/max_weight/'+k,np.array(support['maximum_normalized_weight']),np.array(norm.max()))
            for qq in model.q:
                checks.near(label+'/weighted_q/'+k+'/'+str(qq),np.array(support['q_fractions'][str(qq)]),np.array(norm[p['q']==qq].sum()))
            initial=w[:20000];initial=initial/initial.sum()
            checks.near(label+'/initial_ESS/'+k,np.array(meta['initial_support'][k]['ess']),np.array(1/np.sum(initial*initial)))
        need_larger=min(v['ess'] for v in meta['initial_support'].values())<5000
        checks.yes(label+'/declared_library_size',len(p['q'])==(40000 if need_larger else 20000))
        checks.yes(label+'/minimum_final_ESS',all(v['ess']>=5000 for v in meta['weights'].values()))
        accepted=sum(r['accepted'] for r in meta['proposal_account'])
        checks.yes(label+'/proposal_account',accepted==meta['retained_after_final_batch'] and accepted-len(p['q'])==meta['final_batch_discarded'] and
                   all(r['accepted']+r['resolved_rejected']+r['projected_final_rejected']==r['invariant_and_projected_accepted']<=r['proposals'] for r in meta['proposal_account']))
        # All input hashes, including the three paired observer products, are
        # checked; the observer GLS reconstruction belongs to its own verifier.
        for scenario,op in item['observations'].items():
            op=Path(op);op=op if op.is_absolute() else root/op
            checks.yes(label+'/observation_hash/'+scenario,sha(op)==item['sha256']['observations'][scenario])
    checks.yes('bank_record_set',set(keys)==expected and len(keys)==len(expected))
    checks.yes('independent_role_seed_namespaces',len(seeds)==44 and len(set(seeds))==44 and max(seeds)<268000000)
    return rows


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,default=Path(__file__).parents[1])
    parser.add_argument('--out',type=Path,default=Path(__file__).with_name('population_verification.json'))
    parser.add_argument('--stage',choices=('trajectory','banks','all'),default='all')
    parser.add_argument('--wait',action='store_true')
    parser.add_argument('--max-wait-seconds',type=float,default=7200.)
    args=parser.parse_args();root=args.root.resolve();start=time.time()
    required=[root/'physics/results/interaction_model.npz',root.parent/'reanalysis/data/observation_arrays.npz']
    if args.stage in ('trajectory','all'):required.append(root/'physical_validation.json')
    if args.stage in ('banks','all'):required.append(root/'bank_manifest.json')
    last_missing=None
    while True:
        missing=[str(p) for p in required if not p.exists()]
        if args.stage in ('banks','all') and (root/'bank_manifest.json').exists():
            try:
                current=json.loads((root/'bank_manifest.json').read_text())
                if len(current.get('records',[]))!=22 or current.get('production_authorized') is not True:
                    missing.append('22 completed/authorized banks; currently '+str(len(current.get('records',[]))))
            except json.JSONDecodeError:missing.append('Bank manifest write is in progress')
        if not missing or not args.wait or time.time()-start>args.max_wait_seconds:break
        if missing!=last_missing:
            print(json.dumps({'status':'waiting_for_artifacts','complete':False,'missing':missing}),flush=True)
            last_missing=missing
        time.sleep(15)
    out={'source_sha256':sha(__file__),'independent_model_source_sha256':sha(Path(__file__).with_name('verify_interpolant.py')),
         'scope':'All saved accepted rows and saved trajectory invariant statistics; no production import; rejected proposals/RNG stream not replayed.',
         'stage':args.stage,'complete':False}
    if missing:out.update(status='waiting_for_artifacts',missing=missing)
    else:
        model=IndependentModel(required[0]);model.empirical=load(required[1]);checks=Checks();rows=[]
        if args.stage in ('trajectory','all'):rows.extend(verify_trajectories(root,model,checks))
        if args.stage in ('banks','all'):rows.extend(verify_banks(root,model,checks))
        out.update(complete=True,passed=all(r['passed'] for r in checks.rows),checks=checks.rows,
                   checks_count=len(checks.rows),populations=rows,total_accepted_rows=sum(r['rows'] for r in rows),
                   model_sha256=sha(model.path),empirical_source_sha256=sha(required[1]))
    out['seconds']=time.time()-start;args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('checks','populations')}),flush=True)

if __name__=='__main__':main()
