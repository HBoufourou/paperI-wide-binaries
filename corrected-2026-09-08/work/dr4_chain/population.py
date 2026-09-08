"""Energy/Lz-defined weighted parents and explicit inner-photocentre populations."""
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np
from orbits import units,tangent_basis,project,invariant_barrier
from joint_observer import random_rotations,kepler_eccentric_anomaly,inner_positions

HERE=Path(__file__).resolve().parent
SYSROOT=HERE.parent/'originality/dr4'
sys.path.insert(0,str(SYSROOT))
import run_comparison as empirical

QS=np.array([.1,.3,1.])

@lru_cache(maxsize=1)
def covariate_pool():
    return empirical.source_data()

def unit_directions(n,rng):
    x=rng.normal(size=(n,3));return x/np.linalg.norm(x,axis=1)[:,None]

def inner_orbits(n,rng,source_masses,rlo_au,distance,prior,*,host=None,qin=None):
    if host is None:host=rng.integers(0,2,n)
    if qin is None:qin=rng.uniform(.02,1,n)
    host=np.asarray(host);qin=np.asarray(qin)
    if host.shape!=(n,) or qin.shape!=(n,):raise ValueError('one frozen inner host/mass ratio per system required')
    mass=source_masses[np.arange(n),host]
    inner=dict(active=np.ones(n,bool),host=host,q=qin,mass_host_true=mass,
        a_au=np.empty(n),e=np.empty(n),M0=rng.uniform(0,2*np.pi,n),rotation=random_rotations(n,rng))
    if n==0:return inner,source_masses.copy(),np.zeros(0,bool),0
    todo=np.arange(n);attempts=0
    for _ in range(10000):
        if not len(todo):break
        count=len(todo);attempts+=count
        a=10**rng.normal(np.log10(prior),1.5,count)
        e=(-.4+np.sqrt(.16+2.4*rng.random(count)))/1.2
        valid=(a>=.01)&(a<=300)&(a*(1+e)<rlo_au[todo]/20)
        ids=todo[valid];inner['a_au'][ids]=a[valid];inner['e'][ids]=e[valid];todo=todo[~valid]
    if len(todo):raise RuntimeError('No inner prior support within hierarchy after declared rejection limit')
    e=inner['e'];a=inner['a_au'];anomaly=kepler_eccentric_anomaly(inner['M0'],e)
    plane=np.stack((a*(np.cos(anomaly)-e),a*np.sqrt(1-e*e)*np.sin(anomaly)),axis=1)
    xyz=np.einsum('nij,nj->ni',inner['rotation'][:,:,:2],plane)
    unresolved=np.linalg.norm(xyz[:,:2],axis=1)<distance
    photo=source_masses.copy()
    photo[np.arange(n),host]=mass*(1+qin**4)**.25/(1+qin)
    return inner,photo,unresolved,attempts

def proposal(n,rng,potential,cmax,component):
    """potential(x,q) returns negative relative interaction potential in batch."""
    d,pool=covariate_pool();index=rng.choice(pool,n,replace=True)
    phot_anchor=d['M1'][index]+d['M2'][index]
    q=QS[rng.integers(0,3,n)];fractions=np.column_stack((1/(1+q),q/(1+q)))
    if component=='pure':
        host_draw=np.zeros(n,int);qin_draw=np.zeros(n);flux_factor=np.ones(n)
    else:
        host_draw=rng.integers(0,2,n);qin_draw=rng.uniform(.02,1,n)
        g=(1+qin_draw**4)**.25/(1+qin_draw)
        flux_factor=1-fractions[np.arange(n),host_draw]*(1-g)
    total=phot_anchor/flux_factor
    u=units(total);radius=np.exp(rng.uniform(np.log(.03),np.log(30),n))
    x=radius[:,None]*unit_directions(n,rng);psi=potential(x,q)
    if psi.shape!=(n,) or not np.isfinite(psi).all() or np.any(psi>=0):raise ValueError('invalid supplied interaction potential')
    fraction_energy=rng.beta(2.5,1.5,n);binding=-psi*fraction_energy
    velocity=np.sqrt(2*(-psi)*(1-fraction_energy))[:,None]*unit_directions(n,rng)
    lz=np.cross(x,velocity)[:,2];rlo,valid_barrier=invariant_barrier(-binding,lz,cmax)
    valid=valid_barrier&(binding>=.1)&(binding<=5)&(rlo>.035)&(binding>cmax/30)
    sky=tangent_basis(d['ra1'][index],d['dec1'][index]);rproj=project(x,sky)[:,:2]*u['rM_au'][:,None]/1000
    s=np.linalg.norm(rproj,axis=1)
    # An unresolved companion has separation <distance AU and |photocentre
    # coefficient|<=.5. This necessary padded cut cannot remove a selected row.
    margin=np.zeros(n) if component=='pure' else .5*d['d'][index]/1000
    valid&=(s>=2-margin)&(s<=30+margin)
    ids=np.flatnonzero(valid)
    errors=np.stack((np.stack((d['mu1ra_err'][index],d['mu1dec_err'][index]),axis=1),
                     np.stack((d['mu2ra_err'][index],d['mu2dec_err'][index]),axis=1)),axis=1)
    p=dict(position=x[ids],velocity=velocity[ids],q=q[ids],energy=-binding[ids],lz=lz[ids],
        barrier_rlo=rlo[ids],mass_fraction=fractions[ids],Mtrue_total=total[ids],
        Mphot_total_anchor=phot_anchor[ids],true_mass_flux_factor=flux_factor[ids],
        rM_au=u['rM_au'][ids],rproj_bary_kau=rproj[ids],distance_pc=d['d'][index[ids]],
        ra_deg=d['ra1'][index[ids]],dec_deg=d['dec1'][index[ids]],pm_err_masyr=errors[ids],
        empirical_row_index=index[ids],Mtrue_source=fractions[ids]*total[ids,None],
        weight_beta1=radius[ids]**3*(-psi[ids])**3,
        weight_beta1p6=radius[ids]**3*(-psi[ids])**3*binding[ids]**.6,
        weight_beta1p3=radius[ids]**3*(-psi[ids])**3*binding[ids]**.3)
    nvalid=len(ids)
    if component=='pure':
        photo=p['Mtrue_source'];unresolved=np.ones(nvalid,bool);attempts=0
        inner=dict(active=np.zeros(nvalid,bool),host=np.zeros(nvalid,int),q=np.zeros(nvalid),
            mass_host_true=photo[:,0],a_au=np.ones(nvalid),e=np.zeros(nvalid),M0=np.zeros(nvalid),
            rotation=np.tile(np.eye(3),(nvalid,1,1)))
    else:
        prior=int(component.removeprefix('comp'))
        inner,photo,unresolved,attempts=inner_orbits(nvalid,rng,p['Mtrue_source'],p['barrier_rlo']*p['rM_au'],p['distance_pc'],prior,
            host=host_draw[ids],qin=qin_draw[ids])
    photocentre=inner_positions(inner,np.array([0.]),p['distance_pc'])
    offset=(photocentre[:,1,0]-photocentre[:,0,0])*p['distance_pc'][:,None]/1e6
    p.update(Mphot1=photo[:,0],Mphot2=photo[:,1],rproj_kau=p['rproj_bary_kau']+offset,
             initial_photocentre_offset_kau=offset)
    observed_s=np.linalg.norm(p['rproj_kau'],axis=1)
    selected=unresolved&(observed_s>=2)&(observed_s<=30)
    p.update({'inner_'+key:value for key,value in inner.items()})
    p={key:value[selected] for key,value in p.items()}
    return p,dict(proposals=n,invariant_and_projected_accepted=nvalid,resolved_rejected=int((~unresolved).sum()),
                  inner_prior_attempts=attempts,projected_final_rejected=int(np.sum(unresolved&~selected)),accepted=int(selected.sum()))

def generate_library(n,seed,potential,cmax,component):
    rng=np.random.default_rng(seed);chunks=[];count=0;account=[]
    while count<n:
        p,record=proposal(4096,rng,potential,cmax,component);account.append(record)
        if len(p['q']):chunks.append(p);count+=len(p['q'])
        if sum(r['proposals'] for r in account)>10000000:raise RuntimeError('Insufficient declared parent support')
    output={key:np.concatenate([p[key] for p in chunks],axis=0)[:n] for key in chunks[0]}
    metadata=dict(seed=seed,requested_rows=n,proposal_account=account,
        retained_after_final_batch=count,final_batch_discarded=count-n,
        mass_ratio_counts={str(q):int(np.sum(output['q']==q)) for q in QS},
        weights={})
    for key in ('weight_beta1','weight_beta1p6','weight_beta1p3'):
        w=output[key]/output[key].sum()
        metadata['weights'][key]=dict(ess=float(1/np.sum(w*w)),maximum_normalized_weight=float(w.max()),
            q_fractions={str(q):float(w[output['q']==q].sum()) for q in QS})
    return output,metadata
