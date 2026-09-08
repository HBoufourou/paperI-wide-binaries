"""Measured-only information variants; no simulator labels or energy fields."""
import numpy as np

MODES=('full66','no_diagnostics66','full34','covariates_only66','conditional66')
BINS=(2048,512,2048,8,2048)
GM_SUN=1.3271244e20;AU=149597870700.;A0=1.2e-10;K=4.740470463533349
VC=np.sqrt(GM_SUN/AU)/1000
BASE_FIELDS=('distance_pc','Mphot1','Mphot2','rproj_kau')
OBS_FIELDS={
    'full66':('pm34_masyr','pm66_masyr','pm_source34_masyr','pm_source66_masyr',
        'cov66_masyr2','cov_delta_source_masyr2','cov_delta_relative_masyr2',
        'residual_score_source34','residual_score_source66','acceleration_score_source34','acceleration_score_source66'),
    'no_diagnostics66':('pm66_masyr','cov66_masyr2'),
    'full34':('pm34_masyr','cov34_masyr2','residual_score_source34','acceleration_score_source34'),
    'covariates_only66':('cov66_masyr2',),
}
OBS_FIELDS['conditional66']=OBS_FIELDS['full66']


def finite(value,shape,name,positive=False,nonnegative=False):
    a=np.asarray(value,float)
    if a.shape!=shape or not np.isfinite(a).all():raise ValueError('Invalid '+name+' shape or finite support')
    if positive and np.any(a<=0):raise ValueError('Nonpositive '+name)
    if nonnegative and np.any(a<0):raise ValueError('Negative '+name)
    return a


def covariance(value,shape,name):
    a=finite(value,shape,name)
    if not np.allclose(a,a.swapaxes(-1,-2),rtol=1e-10,atol=1e-15):raise ValueError('Asymmetric '+name)
    if np.any(np.linalg.eigvalsh(a)<=0):raise ValueError('Nonpositive '+name)
    return a


def measured_features(population,observations,mode):
    if mode not in MODES:raise ValueError('Unknown mode')
    d=np.asarray(population['distance_pc'],float)
    if d.ndim!=1 or not len(d):raise ValueError('Nonempty distance vector required')
    n=len(d);d=finite(d,(n,),'distance',positive=True)
    m1=finite(population['Mphot1'],(n,),'Mphot1',positive=True)
    m2=finite(population['Mphot2'],(n,),'Mphot2',positive=True)
    r=finite(population['rproj_kau'],(n,2),'separation vector');s=np.linalg.norm(r,axis=1)
    if np.any(s<=0):raise ValueError('Nonpositive separation')
    m=m1+m2;rm=np.sqrt(GM_SUN*m/A0)/(1000*AU)
    factor=K*d/1000/(VC*np.sqrt(m/(1000*s)))
    baseline=34 if mode=='full34' else 66
    cov=covariance(observations[f'cov{baseline}_masyr2'],(n,2,2),'PM covariance')
    precision=np.sqrt(np.trace(cov,axis1=1,axis2=2)/2)*factor
    result={'separation':s/rm,'precision':precision}
    if not all(np.isfinite(x).all() for x in result.values()):raise ValueError('Nonfinite derived covariate')
    if mode=='covariates_only66':return result
    mu=finite(observations[f'pm{baseline}_masyr'],(n,2),'proper motion')
    vel=mu*factor[:,None];unit=r/s[:,None]
    result.update(parallel=np.einsum('ni,ni->n',vel,unit),perpendicular=vel[:,1]*unit[:,0]-vel[:,0]*unit[:,1])
    if not all(np.isfinite(x).all() for x in result.values()):raise ValueError('Nonfinite derived velocity')
    if mode=='no_diagnostics66':return result
    periods=(34,) if mode=='full34' else (34,66)
    scores=[finite(observations[f'{name}_source{b}'],(n,2),f'{name}{b}',nonnegative=True)
            for name in ('residual_score','acceleration_score') for b in periods]
    if mode!='full34':
        pm34=finite(observations['pm34_masyr'],(n,2),'PM34')
        a=finite(observations['pm_source34_masyr'],(n,2,2),'source PM34')
        b=finite(observations['pm_source66_masyr'],(n,2,2),'source PM66')
        source=covariance(observations['cov_delta_source_masyr2'],(n,2,2,2),'source PM difference covariance')
        relative=covariance(observations['cov_delta_relative_masyr2'],(n,2,2),'relative PM difference covariance')
        quadratic=lambda v,c:np.einsum('...i,...i->...',v,np.linalg.solve(c,v[...,None])[...,0])
        scores.extend([quadratic(b-a,source),quadratic(mu-pm34,relative)[:,None]])
    result['diagnostic']=np.max(np.column_stack(scores),axis=1)
    if not all(np.isfinite(x).all() for x in result.values()):raise ValueError('Nonfinite derived feature')
    return result


def encode(x,mode):
    sep=np.searchsorted([.3,1.,3.],x['separation'],side='right')
    noise=(x['precision']>=.05).astype(int)
    if mode=='covariates_only66':return (2*sep+noise).astype(np.uint16)
    edges=[-np.sqrt(2),-1.,-.5,0.,.5,1.,np.sqrt(2)]
    par=np.searchsorted(edges,x['parallel'],side='right');perp=np.searchsorted(edges,x['perpendicular'],side='right')
    v=(sep*8+par)*8+perp
    if mode=='no_diagnostics66':return (v*2+noise).astype(np.uint16)
    diag=np.searchsorted([3.,10.,100.],x['diagnostic'],side='right')
    return ((v*4+diag)*2+noise).astype(np.uint16)


def measured_codes(population,observations,mode):
    return encode(measured_features(population,observations,mode),mode)


def project_codes(full,mode):
    full=np.asarray(full)
    if mode in ('full66','conditional66'):return full
    if mode=='no_diagnostics66':return ((full//8)*2+full%2).astype(np.uint16)
    if mode=='covariates_only66':return ((full//512)*2+full%2).astype(np.uint16)
    raise ValueError('Full34 cannot be projected from full66')


def project_histogram(histogram,mode):
    a=np.asarray(histogram)
    if a.shape[-1]!=2048:raise ValueError('Expected full histogram')
    a=a.reshape(a.shape[:-1]+(4,8,8,4,2))
    if mode=='no_diagnostics66':return a.sum(axis=-2).reshape(a.shape[:-5]+(512,))
    if mode=='covariates_only66':return a.sum(axis=(-4,-3,-2)).reshape(a.shape[:-5]+(8,))
    raise ValueError('Expected a marginal mode')
