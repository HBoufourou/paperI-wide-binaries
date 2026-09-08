"""Run the fixed paired observation experiment described in PROTOCOL.md."""
from pathlib import Path
import argparse,hashlib,json,platform,time
import numpy as np
import pandas as pd
from epoch_operator import observe_pairs,SCENARIOS,random_rotations,kepler_eccentric_anomaly,K

HERE=Path(__file__).resolve().parent
DATA=HERE.parents[1]/'reanalysis'/'data'/'observation_arrays.npz'

def source_data():
    with np.load(DATA) as f:d={k:f[k] for k in f.files}
    m=(d['d']<200)&(d['R_chance']<.01)&(d['ruwe1']<1.4)&(d['ruwe2']<1.4)
    m&=(np.maximum(d['sigmax'],d['sigmay'])<.1)&d['has_rv']&np.isfinite(d['sigmaRV'])
    m&=np.hypot(d['vx_flat'],d['vy_flat'])<=2.23/np.sqrt(d['Mtot'])
    m&=((d['s_kau']>=.2)&(d['s_kau']<30))|((d['gN_a0']>=.03)&(d['gN_a0']<.3))
    return d,np.flatnonzero(m)

def sampled_covariates(n,seed):
    d,pool=source_data();rng=np.random.default_rng(seed);idx=rng.choice(pool,n,replace=True)
    errors=np.stack((np.stack((d['mu1ra_err'][idx],d['mu1dec_err'][idx]),axis=1),
                     np.stack((d['mu2ra_err'][idx],d['mu2dec_err'][idx]),axis=1)),axis=1)
    return {'row_index':idx,'distance_pc':d['d'][idx],'pm_err_masyr':errors,
        'ra_deg':d['ra1'][idx],'dec_deg':d['dec1'][idx],
        'Mphot1':d['M1'][idx],'Mphot2':d['M2'][idx],'s_kau':d['s_kau'][idx],
        'source_id1':d['source_id1'][idx],'source_id2':d['source_id2'][idx],
        'zone_V':(d['s_kau'][idx]>=.2)&(d['s_kau'][idx]<2),
        'zone_T':(d['s_kau'][idx]>=2)&(d['s_kau'][idx]<30),
        'zone_D':(d['gN_a0'][idx]>=.03)&(d['gN_a0'][idx]<.3)}

def unit_vectors(n,rng):
    x=rng.normal(size=(n,3));return x/np.linalg.norm(x,axis=1)[:,None]

def population(n=30000,seed=72401):
    p=sampled_covariates(n,seed);rng=np.random.default_rng(seed+1)
    r=2*rng.beta(2.5,1.5,n);rvec=unit_vectors(n,rng)*r[:,None]
    vv=unit_vectors(n,rng)*np.sqrt(2/r-1)[:,None]
    a=1000*p['s_kau']/np.linalg.norm(rvec[:,:2],axis=1)
    eccentricity=np.sqrt(np.maximum(0,1-np.sum(np.cross(rvec,vv)**2,axis=1)))
    active=rng.random(n)<.30;host=rng.integers(0,2,n)
    # A tiny subset of near-radial outer draws admits no inner orbit above .01AU.
    # Condition active systems on existence of the declared hierarchical support.
    bad=np.flatnonzero(active&(.342*a*(1-eccentricity)**2<=.011))
    while len(bad):
        rb=2*rng.beta(2.5,1.5,len(bad));rr=unit_vectors(len(bad),rng)*rb[:,None]
        vb=unit_vectors(len(bad),rng)*np.sqrt(2/rb-1)[:,None]
        vv[bad]=vb;a[bad]=1000*p['s_kau'][bad]/np.linalg.norm(rr[:,:2],axis=1)
        eccentricity[bad]=np.sqrt(np.maximum(0,1-np.sum(np.cross(rr,vb)**2,axis=1)))
        bad=bad[.342*a[bad]*(1-eccentricity[bad])**2<=.011]
    q=rng.uniform(.02,1,n);mapp=np.where(host==0,p['Mphot1'],p['Mphot2'])
    mh=mapp/(1+q**4)**.25;mc=q*mh
    inner={'active':active,'host':host,'q':q,'mass_host_true':mh+mc,
        'a_au':np.ones(n),'e':np.zeros(n),'M0':np.zeros(n),'rotation':np.tile(np.eye(3),(n,1,1))}
    todo=np.flatnonzero(active)
    for _ in range(10000):
        nn=len(todo)
        if not nn:break
        ai=10**rng.normal(np.log10(40),1.5,nn)
        ei=(-.4+np.sqrt(.16+2.4*rng.random(nn)))/1.2
        ma=rng.uniform(0,2*np.pi,nn);rot=random_rotations(nn,rng)
        E=kepler_eccentric_anomaly(ma,ei)
        xy=np.stack((ai*(np.cos(E)-ei),ai*np.sqrt(1-ei**2)*np.sin(E)),axis=1)
        xyz=np.einsum('nij,nj->ni',rot[:,:,:2],xy)
        good=(ai>=.01)&(ai<=300)&(ai<.342*a[todo]*(1-eccentricity[todo])**2)
        good&=np.linalg.norm(xyz[:,:2],axis=1)/p['distance_pc'][todo]<1
        ids=todo[good]
        for key,val in [('a_au',ai),('e',ei),('M0',ma),('rotation',rot)]:inner[key][ids]=val[good]
        todo=todo[~good]
    if len(todo):raise RuntimeError('inner eligibility rejection failed')
    mphoto=p['Mphot1']+p['Mphot2'];mtrue=mphoto+np.where(active,mh+mc-mapp,0.)
    wide=vv[:,:2]*29.784691831696804*np.sqrt(mtrue/a)[:,None]
    vc=29.784691831696804*np.sqrt(mphoto/(1000*p['s_kau']))
    p.update(Mphot_total=mphoto,Mtrue_total=mtrue,wide_v_kms=wide,
        true_wide_vx=wide[:,0]/vc,true_wide_vy=wide[:,1]/vc,
        vc_kms=vc,outer_a_au=a,outer_e=eccentricity,has_companion=active,
        period_year=np.where(active,np.sqrt(inner['a_au']**3/inner['mass_host_true']),np.nan))
    p.update({f'inner_{k}':v for k,v in inner.items()})
    return p,inner

def observation(p,inner,scenario,seed):
    return observe_pairs(p.get('wide_v_kms'),p['distance_pc'],p['pm_err_masyr'],inner,
        scenario,seed,p['ra_deg'],p['dec_deg'])

def wilson(k,n):
    if n==0:return [None,None]
    pp=k/n;z=1.959963984540054;den=1+z*z/n
    mid=(pp+z*z/(2*n))/den;half=z*np.sqrt(pp*(1-pp)/n+z*z/(4*n*n))/den
    return [float(mid-half),float(mid+half)]

METHODS={'34_residual':('residual_score34',34),'66_residual':('residual_score66',66),
         '66_acceleration':('acceleration_score66',66),'66_epoch_combined':('combined_score66',66)}

def calibrate(scenario,ntrain=20000,ntest=20000):
    a=sampled_covariates(ntrain,51001);b=sampled_covariates(ntest,51002)
    train=observation(a,None,scenario,61001);test=observation(b,None,scenario,61002)
    result={};rows=[]
    for method,(score,baseline) in METHODS.items():
        threshold=float(np.quantile(train[score],.99,method='higher'))
        k=int(np.sum(test[score]>threshold));ci=wilson(k,ntest)
        result[method]=threshold
        rows.append(dict(scenario=scenario,method=method,threshold=threshold,n_train=ntrain,n_test=ntest,
            false_flags=k,false_alarm=k/ntest,false_alarm_wilson_lo=ci[0],false_alarm_wilson_hi=ci[1]))
    # Null normalized residuals test the covariance returned by the actual operator.
    for baseline in (34,66):
        err=test[f'bias_noise{baseline}_masyr'];cov=test[f'cov{baseline}_masyr2']
        maha=np.einsum('ni,nij,nj->n',err,np.linalg.inv(cov),err)
        rows.append(dict(scenario=scenario,method=f'null_covariance_{baseline}',n_test=ntest,
            mahalanobis_mean=float(np.mean(maha)),ellipse95_coverage=float(np.mean(maha<=5.991464547107979))))
    return result,rows

def summarize(p,out,thresholds,scenario):
    active=p['has_companion'];n=len(active);rows=[];periodrows=[]
    scale=K*p['distance_pc']/1000/p['vc_kms']
    true=np.column_stack((p['true_wide_vx'],p['true_wide_vy']))
    tnorm=np.linalg.norm(true,axis=1)
    null=~active;pass34=out['residual_score34']<=thresholds['34_residual']
    stored={}
    for method,(score,baseline) in METHODS.items():
        flag=out[score]>thresholds[method];keep=~flag
        observed=out[f'pm{baseline}_masyr']*scale[:,None]
        ib=out[f'inner_bias{baseline}_masyr']*scale[:,None]
        err=out[f'bias_noise{baseline}_masyr']*scale[:,None]
        norm=np.linalg.norm(observed,axis=1)
        sig=np.sqrt(np.trace(out[f'cov{baseline}_masyr2'],axis1=1,axis2=2)/2)*scale
        stored[f'flag_{method}']=flag
        for zone,zmask in [('ALL',np.ones(n,bool)),('V',p['zone_V']),('T',p['zone_T']),('D',p['zone_D'])]:
            zkeep=zmask&keep;zn=int(zmask.sum());nr=int(zkeep.sum())
            if nr==0:continue
            trueidx=np.median(tnorm[zkeep]);obsidx=np.median(norm[zkeep])
            row=dict(scenario=scenario,method=method,zone=zone,N_input=zn,N_retained=nr,
                retention=nr/zn,initial_companion_fraction=float(active[zmask].mean()),
                companion_detection=float(flag[zmask&active].mean()) if np.any(zmask&active) else np.nan,
                retained_companion_fraction=float(active[zkeep].mean()),
                median_norm_observed=float(obsidx),median_norm_true_retained=float(trueidx),
                squared_median_bias_index=float((obsidx/trueidx)**2),
                median_formal_sigma_normalized=float(np.median(sig[zkeep])),
                error_rms_normalized=float(np.sqrt(np.mean(np.sum(err[zkeep]**2,axis=1)))),
                inner_bias_rms_normalized=float(np.sqrt(np.mean(np.sum(ib[zkeep]**2,axis=1)))),
                tail_above_sqrt2=float(np.mean(norm[zkeep]>np.sqrt(2))),
                inherited34_companions_detected=float(flag[zmask&active&pass34].mean()) if np.any(zmask&active&pass34) else np.nan)
            rows.append(row)
        edges=[0,1,3,10,30,100,np.inf]
        for lo,hi in zip(edges[:-1],edges[1:]):
            mask=active&(p['period_year']>=lo)&(p['period_year']<hi);nm=int(mask.sum())
            kk=int(flag[mask].sum());ci=wilson(kk,nm)
            old=mask&pass34;nold=int(old.sum())
            periodrows.append(dict(scenario=scenario,method=method,period_lo=lo,period_hi=hi,N=nm,
                detected=kk,detection_fraction=kk/nm if nm else np.nan,detection_wilson_lo=ci[0],detection_wilson_hi=ci[1],
                N_survived34=nold,detected_among34_survivors=int(flag[old].sum()),
                fraction_among34_survivors=float(flag[old].mean()) if nold else np.nan))
    for baseline in (34,66):
        stored[f'observed_vx{baseline}']=out[f'pm{baseline}_masyr'][:,0]*scale
        stored[f'observed_vy{baseline}']=out[f'pm{baseline}_masyr'][:,1]*scale
        stored[f'cov_normalized{baseline}']=out[f'cov{baseline}_masyr2']*scale[:,None,None]**2
    return rows,periodrows,stored

def run(args):
    start=time.time();dest=HERE/'results';dest.mkdir(exist_ok=True)
    p,inner=population(args.n,72401)
    np.savez_compressed(dest/'latent_population.npz',**p)
    summaries=[];periods=[];calrows=[];thresholds={}
    for scenario in args.scenarios.split(','):
        print('calibrate',scenario,flush=True)
        th,cr=calibrate(scenario,args.ntrain,args.ntest);thresholds[scenario]=th;calrows+=cr
        # Save thresholds before generating/examining the companion outcomes.
        (dest/f'thresholds_{scenario}.json').write_text(json.dumps(th,indent=2))
        print('observe paired population',scenario,flush=True)
        out=observation(p,inner,scenario,82401)
        sr,pr,stored=summarize(p,out,th,scenario);summaries+=sr;periods+=pr
        np.savez_compressed(dest/f'observed_{scenario}.npz',**out,**stored,id=np.arange(args.n))
        pd.DataFrame(summaries).to_csv(dest/'comparison.csv',index=False)
        pd.DataFrame(periods).to_csv(dest/'period_detection.csv',index=False)
        pd.DataFrame(calrows).to_csv(dest/'null_validation.csv',index=False)
        print(pd.DataFrame(sr).query("zone=='ALL'")[['method','N_retained','companion_detection','retained_companion_fraction','squared_median_bias_index']].to_string(index=False),flush=True)
    manifest={'started_after_protocol':True,'protocol_sha256':hashlib.sha256((HERE/'PROTOCOL.md').read_bytes()).hexdigest(),
        'code_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [HERE/'epoch_operator.py',HERE/'run_comparison.py']},
        'input_sha256':hashlib.sha256(DATA.read_bytes()).hexdigest(),'arguments':vars(args),
        'python':platform.python_version(),'numpy':np.__version__,'thresholds':thresholds,'elapsed_seconds':time.time()-start,
        'limitations':'Controlled cadence/noise/companion scenarios; no actual DR4 data, Gaia selection reconstruction or gravity verdict.'}
    (dest/'manifest.json').write_text(json.dumps(manifest,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=30000)
    ap.add_argument('--ntrain',type=int,default=20000);ap.add_argument('--ntest',type=int,default=20000)
    ap.add_argument('--scenarios',default=','.join(SCENARIOS));run(ap.parse_args())
