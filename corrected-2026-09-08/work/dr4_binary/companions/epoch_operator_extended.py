"""Independent paired along-scan observation scenarios; not a Gaia pipeline.

Public API: observe_pairs(wide_v_kms, distance_pc, pm_err_masyr, inner,
scenario, seed, ra_deg=None, dec_deg=None). Units: km/s, pc, mas/yr.
pm_err_masyr has shape(N,2 sources,2 axes). inner is None or a dictionary
with a_au,e,M0,mass_host_true,q,host,rotation(N,3,3),active arrays.
The relative convention is source2-source1; angular quantities are mas.
"""
from dataclasses import dataclass,asdict
import math
import numpy as np

T3=34/12
T4=66/12
K=4.740470463533349

@dataclass(frozen=True)
class Scenario:
    name:str
    rate:float=12.
    clustered:bool=False
    anisotropic_angles:bool=False
    temporal_rho:float=0.
    temporal_scale:float=.03
    pair_rho:float=0.
    pm_floor_fraction:float=0.

SCENARIOS={
    'uniform_white':Scenario('uniform_white'),
    'clustered_correlated':Scenario('clustered_correlated',clustered=True,temporal_rho=.25,pair_rho=.2),
    'sparse_floor':Scenario('sparse_floor',rate=8.,clustered=True,anisotropic_angles=True,temporal_rho=.25,pair_rho=.2,pm_floor_fraction=.2),
}

def cadence(scenario):
    s=SCENARIOS[scenario] if isinstance(scenario,str) else scenario
    def interval(lo,hi):
        n=int(round((hi-lo)*s.rate))
        if s.clustered:
            n=max(8,2*int(round(n/2)))
            c=np.linspace(lo+.05,hi-.05,n//2)
            t=np.column_stack((c,c+1.78/24/365.25)).ravel()
            k=np.repeat(np.arange(n//2),2)
        else:
            t=np.linspace(lo+.025,hi-.025,n)
            k=np.arange(n)
        angle=np.mod(k*np.pi*(3-np.sqrt(5)),2*np.pi)
        if s.anisotropic_angles:angle=np.where(k%2,1.05,0.)+.12*np.sin(k*1.37)
        return t,angle
    t3,a3=interval(0,T3);te,ae=interval(T3,T4)
    # Distinct angles after the shared prefix without changing either first epoch.
    ae=ae+.47
    return np.r_[t3,te],np.r_[a3,ae],len(t3)

def kepler_eccentric_anomaly(M,e):
    """Bracketed Newton solver, independent of the historical generator."""
    m=np.mod(M,2*np.pi);ee=np.broadcast_to(e,m.shape)
    lo=np.zeros_like(m);hi=np.full_like(m,2*np.pi)
    x=m+np.sign(np.sin(m))*.85*ee
    for _ in range(48):
        f=x-ee*np.sin(x)-m
        hi=np.where(f>0,x,hi);lo=np.where(f<0,x,lo)
        step=x-f/(1-ee*np.cos(x))
        step=np.where((step<=lo)|(step>=hi),.5*(lo+hi),step)
        step=np.where(np.abs(f)<1e-13,x,step)
        if np.max(np.abs(step-x))<2e-13:return step
        x=step
    if np.max(np.abs(x-ee*np.sin(x)-m))>2e-11:raise RuntimeError('Kepler solve failed')
    return x

def random_rotations(n,rng):
    q=rng.normal(size=(n,4));q/=np.linalg.norm(q,axis=1)[:,None]
    w,x,y,z=q.T
    return np.stack((1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w),
        2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w),
        2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)),axis=1).reshape(n,3,3)

def inner_positions(inner,t,distance_pc):
    """Per-source photocentre offsets (N,2 sources,n_epochs,2 axes), mas."""
    n=len(distance_pc);out=np.zeros((n,2,len(t),2))
    if inner is None:return out
    ids=np.flatnonzero(inner['active'])
    if not len(ids):return out
    a=inner['a_au'][ids];e=inner['e'][ids]
    period=np.sqrt(a**3/inner['mass_host_true'][ids])
    E=kepler_eccentric_anomaly(inner['M0'][ids,None]+2*np.pi*t[None,:]/period[:,None],e[:,None])
    plane=np.stack((a[:,None]*(np.cos(E)-e[:,None]),a[:,None]*np.sqrt(1-e[:,None]**2)*np.sin(E)),axis=2)
    xy=np.einsum('nij,ntj->nti',inner['rotation'][ids,:2,:2],plane)
    q=inner['q'][ids];coefficient=q**4/(1+q**4)-q/(1+q)
    xy*=coefficient[:,None,None]*1000/distance_pc[ids,None,None]
    out[ids,inner['host'][ids]]=xy
    return out

def astrometric_design(t,angle,ra_deg,dec_deg,acceleration=False):
    """Circular Earth approximation, ICRS tangent-plane parallax factors."""
    ra=np.deg2rad(ra_deg);de=np.deg2rad(dec_deg)
    earth=np.stack((np.cos(2*np.pi*t),np.sin(2*np.pi*t)*np.cos(np.deg2rad(23.43928)),
        np.sin(2*np.pi*t)*np.sin(np.deg2rad(23.43928))),axis=1)
    ura=np.stack((-np.sin(ra),np.cos(ra),np.zeros_like(ra)),axis=1)
    ude=np.stack((-np.cos(ra)*np.sin(de),-np.sin(ra)*np.sin(de),np.cos(de)),axis=1)
    pra=-ura@earth.T;pde=-ude@earth.T
    c=np.cos(angle);ss=np.sin(angle);tc=t-.5*(t.min()+t.max())
    x=np.empty((len(ra),len(t),7 if acceleration else 5))
    x[:,:,0]=c;x[:,:,1]=ss;x[:,:,2]=tc*c;x[:,:,3]=tc*ss
    x[:,:,4]=pra*c+pde*ss
    if acceleration:x[:,:,5]=.5*tc**2*c;x[:,:,6]=.5*tc**2*ss
    return x

def fit_geometry(t,angle,ra,de,covshape):
    invL=np.linalg.inv(np.linalg.cholesky(covshape))
    x7=astrometric_design(t,angle,ra,de,True)
    xw=np.einsum('ij,njk->nik',invL,x7)
    result={}
    for p in (5,7):
        x=xw[:,:,:p]
        normal=np.einsum('ntp,ntq->npq',x,x)
        inv=np.linalg.inv(normal)
        result[p]=(x,np.einsum('npq,ntq->npt',inv,x),inv)
    result['invL']=invL
    return result

def _scores(resid5,resid7,df):
    # Wilson-Hilferty is only a monotone score transform, never used as a p-value.
    z=((np.maximum(resid5,1e-100)/df)**(1/3)-(1-2/(9*df)))/np.sqrt(2/(9*df))
    tail=np.fromiter((.5*math.erfc(float(v)/np.sqrt(2)) for v in z.ravel()),float,z.size).reshape(z.shape)
    r=-np.log(np.maximum(tail,1e-300))
    a=np.maximum(resid5-resid7,0)/2
    return r.max(axis=1),a.max(axis=1),np.maximum(r,a).max(axis=1),r,a

def observe_pairs(wide_v_kms,distance_pc,pm_err_masyr,inner=None,
                  scenario='uniform_white',seed=1,ra_deg=None,dec_deg=None,batch_size=384):
    """Return paired PM observations, formal covariances and independent flags scores.

    wide_v_kms may be None: additive output bias_noise34/66_masyr does not depend
    on the straight-line wide velocity. No source/sample cuts are applied.
    Returned inner_bias is deterministic given orbit and cadence. Relative signs
    follow source2-source1. Fixed-design noise scale matches mean PM variance,
    not every measured marginal separately. Persistent slopes share both releases.
    """
    s=SCENARIOS[scenario] if isinstance(scenario,str) else scenario
    d=np.asarray(distance_pc,float);err=np.asarray(pm_err_masyr,float);n=len(d)
    if err.shape!=(n,2,2):raise ValueError('pm_err_masyr must be (N,2,2)')
    if np.any(d<=0) or np.any(err<=0) or not np.all(np.isfinite(err)):raise ValueError('Positive finite distances/errors required')
    ra=np.zeros(n) if ra_deg is None else np.asarray(ra_deg,float)
    de=np.full(n,45.) if dec_deg is None else np.asarray(dec_deg,float)
    t,angle,n3=cadence(s);dt=np.abs(t[:,None]-t[None,:])
    C=(1-s.temporal_rho)*np.eye(len(t))+s.temporal_rho*np.exp(-dt/s.temporal_scale)
    L=np.linalg.cholesky(C);rng=np.random.default_rng(seed)
    names=['pm','bias_noise','inner_bias','noise']
    out={f'{k}{b}_masyr':np.empty((n,2)) for k in names for b in (34,66)}
    for b in (34,66):
        out[f'cov{b}_masyr2']=np.empty((n,2,2))
        for k in ('residual_score','acceleration_score','combined_score','residual_chi2_max'):
            out[f'{k}{b}']=np.empty(n)
    # Additional observables and analytic shared-release covariance only.
    for b in (34,66):
        out[f'pm_source{b}_masyr']=np.empty((n,2,2))
        out[f'residual_score_source{b}']=np.empty((n,2))
        out[f'acceleration_score_source{b}']=np.empty((n,2))
    out['cov_delta_source_masyr2']=np.empty((n,2,2,2))
    out['cov_delta_relative_masyr2']=np.empty((n,2,2))
    out['cov_cross_relative_masyr2']=np.empty((n,2,2))
    out['sigma_al_mas']=np.empty((n,2))
    truepm=np.zeros((n,2)) if wide_v_kms is None else np.asarray(wide_v_kms)*1000/(K*d[:,None])
    for start in range(0,n,batch_size):
        sl=slice(start,min(start+batch_size,n));nb=len(d[sl])
        g3=fit_geometry(t[:n3],angle[:n3],ra[sl],de[sl],C[:n3,:n3])
        g4=fit_geometry(t,angle,ra[sl],de[sl],C)
        basevar=np.trace(g3[5][2][:,2:4,2:4],axis1=1,axis2=2)/2
        sigma=np.sqrt(np.mean(err[sl]**2,axis=2)*(1-s.pm_floor_fraction**2)/basevar[:,None])
        out['sigma_al_mas'][sl]=sigma
        common=rng.normal(size=(nb,len(t)))@L.T
        ind=rng.normal(size=(nb,2,len(t)))@L.T
        noise=(np.sqrt(s.pair_rho)*common[:,None,:]+np.sqrt(1-s.pair_rho)*ind)*sigma[:,:,None]
        floor=rng.normal(size=(nb,2,2))*err[sl]*s.pm_floor_fraction
        subinner=None if inner is None else {k:np.asarray(v)[sl] for k,v in inner.items()}
        xy=inner_positions(subinner,t,d[sl]);al=xy[:,:,:,0]*np.cos(angle)+xy[:,:,:,1]*np.sin(angle)
        for b,m,g in ((34,n3,g3),(66,len(t),g4)):
            y=al[:,:,:m]+noise[:,:,:m]
            yw=np.einsum('ij,nsj->nsi',g['invL'],y)
            iw=np.einsum('ij,nsj->nsi',g['invL'],al[:,:,:m])
            fit5=np.einsum('npt,nst->nsp',g[5][1],yw)
            fiti=np.einsum('npt,nst->nsp',g[5][1],iw)
            fit7=np.einsum('npt,nst->nsp',g[7][1],yw)
            res5=yw-np.einsum('ntp,nsp->nst',g[5][0],fit5)
            res7=yw-np.einsum('ntp,nsp->nst',g[7][0],fit7)
            chi5=np.sum(res5**2,axis=2)/sigma**2;chi7=np.sum(res7**2,axis=2)/sigma**2
            rs,ac,co,rs_source,ac_source=_scores(chi5,chi7,m-5)
            offset=fit5[:,1,2:4]-fit5[:,0,2:4]+floor[:,1]-floor[:,0]
            ib=fiti[:,1,2:4]-fiti[:,0,2:4]
            out[f'pm_source{b}_masyr'][sl]=fit5[:,:,2:4]+floor+truepm[sl,None,:]*np.array([-.5,.5])[None,:,None]
            out[f'residual_score_source{b}'][sl]=rs_source
            out[f'acceleration_score_source{b}'][sl]=ac_source
            out[f'pm{b}_masyr'][sl]=truepm[sl]+offset
            out[f'bias_noise{b}_masyr'][sl]=offset
            out[f'inner_bias{b}_masyr'][sl]=ib
            out[f'noise{b}_masyr'][sl]=offset-ib
            scale=sigma[:,0]**2+sigma[:,1]**2-2*s.pair_rho*sigma[:,0]*sigma[:,1]
            cov=scale[:,None,None]*g[5][2][:,2:4,2:4]
            cov[:,0,0]+=np.sum(err[sl,:,0]**2,axis=1)*s.pm_floor_fraction**2
            cov[:,1,1]+=np.sum(err[sl,:,1]**2,axis=1)*s.pm_floor_fraction**2
            out[f'cov{b}_masyr2'][sl]=cov
            out[f'residual_score{b}'][sl]=rs;out[f'acceleration_score{b}'][sl]=ac
            out[f'combined_score{b}'][sl]=co;out[f'residual_chi2_max{b}'][sl]=chi5.max(axis=1)
        B3=np.einsum('npt,tj->npj',g3[5][1][:,2:4,:],g3['invL'])
        B4=np.einsum('npt,tj->npj',g4[5][1][:,2:4,:],g4['invL'])
        cross=np.einsum('npi,ij,nqj->npq',B3,C[:n3,:],B4)
        delta=B4.copy();delta[:,:,:n3]-=B3
        dcov=np.einsum('npi,ij,nqj->npq',delta,C,delta)
        out['cov_delta_source_masyr2'][sl]=sigma[:,:,None,None]**2*dcov[:,None,:,:]
        out['cov_delta_relative_masyr2'][sl]=scale[:,None,None]*dcov
        cross*=scale[:,None,None]
        cross[:,0,0]+=np.sum(err[sl,:,0]**2,axis=1)*s.pm_floor_fraction**2
        cross[:,1,1]+=np.sum(err[sl,:,1]**2,axis=1)*s.pm_floor_fraction**2
        out['cov_cross_relative_masyr2'][sl]=cross
    out['times_year']=t;out['scan_angle_rad']=angle;out['n_epoch34']=np.array(n3)
    out['n_epoch66']=np.array(len(t))
    return out
