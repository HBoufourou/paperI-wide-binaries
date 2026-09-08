"""Finite-mass QUMOND interaction energy, with a controlled spatial tail.

Single process / one BLAS thread. G=Mtotal=a0=1, uniform Newtonian E=1.
"""
import os
for _name in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"]:
    os.environ[_name]="1"
from pathlib import Path
import sys
import time
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"dr4_binary"/"physics"))
from binary_qumond import angular_rule,logarithmic_rule,binary_configuration,nu_parts

LEVELS={"base":dict(log_width=.6,radial_order=6,nmu=24,nphi=48),
        "reference":dict(log_width=.45,radial_order=8,nmu=40,nphi=80)}


def q_phantom(y):
    """Q(y^2)-y^2, Q(0)=0; derivative with respect to y is 2y(nu-1)."""
    y=np.asarray(y,float)
    if np.any(y<0) or not np.all(np.isfinite(y)):
        raise ValueError("nonnegative finite field magnitude required")
    out=np.empty_like(y);small=y<1e-3
    t=y[small]
    out[small]=4*t**1.5/3-.5*t*t+t**2.5/10-t**3.5/224+t**4.5/2304
    t=y[~small];s=np.sqrt(t*(t+4))
    out[~small]=s+2*t*t/(s+t)-4*np.arcsinh(np.sqrt(t)/2)
    return out


def mixed_phantom_difference(u1,u2,ue,weak_threshold=.01):
    """Inclusion/exclusion, using a mixed Hessian quadrature in weak fields."""
    norms=np.linalg.norm(u1,axis=1)+np.linalg.norm(u2,axis=1)
    weak=norms<weak_threshold
    result=np.empty(len(u1))
    strong=~weak
    if np.any(strong):
        a,b=u1[strong],u2[strong]
        result[strong]=(q_phantom(np.linalg.norm(ue+a+b,axis=1))-
                        q_phantom(np.linalg.norm(ue+a,axis=1))-
                        q_phantom(np.linalg.norm(ue+b,axis=1))+
                        q_phantom(np.linalg.norm(ue)))
    if np.any(weak):
        a,b=u1[weak],u2[weak]
        ab=np.sum(a*b,axis=1);total=np.zeros(len(a))
        nodes=[.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)]
        for s in nodes:
            for t in nodes:
                u=ue+s*a+t*b;y=np.linalg.norm(u,axis=1)
                eta,dn=nu_parts(y)
                total += .5*(eta*ab+dn/y*np.sum(u*a,axis=1)*np.sum(u*b,axis=1))
        result[weak]=total
    return result


def newton_interaction_energy(masses,b,separation):
    d=float(separation)
    landmarks=[b[0]/2,b[0],2*b[0],d/2,d,2*d]
    for offset in [b[1],10*b[1],100*b[1]]:landmarks.extend([d-offset,d+offset])
    r,w=logarithmic_rule(min(b)*1e-12,max(d,max(b))*1e6,landmarks,.12,16)
    average=2/(np.sqrt((r+d)**2+b[1]**2)+np.sqrt((r-d)**2+b[1]**2))
    shell=3*b[0]**2*r**3/(r*r+b[0]**2)**2.5
    return float(-np.prod(masses)*np.sum(w*shell*average))


def interaction_energy(q,separation,mu,*,epsilon=.00125,level="base",outer_factor=64.,
                       weak_threshold=.01):
    start=time.perf_counter()
    if q not in [.1,.3,1.] or separation<=0 or not -1<=mu<=1:
        raise ValueError("q in {.1,.3,1}, positive r, signed mu in [-1,1] required")
    masses,positions,b,ext=binary_configuration(q,separation,np.degrees(np.arccos(mu)),epsilon,1.)
    ue=-ext
    settings=LEVELS[level]
    e3=(positions[1]-positions[0])/separation
    side=ext-(ext@e3)*e3
    if np.linalg.norm(side)<1e-10:
        coordinate=np.eye(3)[np.argmin(np.abs(e3))];side=coordinate-(coordinate@e3)*e3
    e1=side/np.linalg.norm(side)
    basis=np.column_stack((e1,np.cross(e3,e1),e3))
    directions,wangle=angular_rule(settings["nmu"],settings["nphi"])
    directions=directions@basis.T
    integral=0.;point_count=0
    rmax=outer_factor*max(1.,separation)
    for centre in [0,1]:
        r,w=logarithmic_rule(1e-12*np.sqrt(masses[centre]),rmax,
                [b[centre]/2,b[centre],2*b[centre],separation/2,separation,2*separation,1.],
                settings["log_width"],settings["radial_order"])
        for first in range(0,len(r),32):
            rr=r[first:first+32]
            x=positions[centre]+(rr[:,None,None]*directions[None,:,:]).reshape(-1,3)
            offsets=[x-positions[i] for i in [0,1]]
            radius2=[np.sum(a*a,axis=1) for a in offsets]
            u=[masses[i]*offsets[i]/(radius2[i]+b[i]**2)[:,None]**1.5 for i in [0,1]]
            partition=radius2[1-centre]**3/(radius2[0]**3+radius2[1]**3)
            volume=(w[first:first+32,None]*rr[:,None]**3*wangle[None,:]).ravel()
            value=mixed_phantom_difference(u[0],u[1],ue,weak_threshold)
            integral+=float(np.sum(volume*partition*value))
            point_count+=len(x)
    phantom_volume=-integral/(8*np.pi)
    eta,dn=nu_parts(np.array([1.]));nu=float(1+eta[0]);ke=float(dn[0]/nu)
    phantom_tail=-float(np.prod(masses))*(nu*(1+ke/3)-1)/rmax
    un=newton_interaction_energy(masses,b,separation)
    energy=un+phantom_volume+phantom_tail
    reduced_mass=float(np.prod(masses))
    return dict(q=float(q),r=float(separation),mu=float(mu),epsilon=float(epsilon),
                level=level,settings=settings,outer_factor=outer_factor,weak_threshold=weak_threshold,
                U=energy,Psi=energy/reduced_mass,Newton_U=un,Newton_Psi=un/reduced_mass,
                phantom_volume=phantom_volume,phantom_spatial_tail=phantom_tail,
                reduced_mass=reduced_mass,quadrature_points=point_count,
                elapsed_seconds=time.perf_counter()-start)


if __name__=="__main__":
    import json
    print(json.dumps(interaction_energy(1.,1.,1/np.sqrt(2)),indent=2),flush=True)
