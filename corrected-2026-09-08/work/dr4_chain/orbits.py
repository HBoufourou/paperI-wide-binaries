"""Relative Hamiltonian trajectories for a supplied conservative batch potential."""
import numpy as np

GM_SUN=1.3271244e20
A0=1.2e-10
AU_M=149597870700.
YEAR_S=31557600.

def units(total_mass):
    total_mass=np.asarray(total_mass,float)
    if np.any(~np.isfinite(total_mass)) or np.any(total_mass<=0):raise ValueError('positive finite masses required')
    radius=np.sqrt(GM_SUN*total_mass/A0)
    return dict(rM_au=radius/AU_M,vM_kms=np.sqrt(A0*radius)/1000,tM_year=np.sqrt(radius/A0)/YEAR_S)

def tangent_basis(ra_deg,dec_deg):
    ra,de=np.deg2rad(ra_deg),np.deg2rad(dec_deg)
    east=np.stack((-np.sin(ra),np.cos(ra),np.zeros_like(ra)),axis=-1)
    north=np.stack((-np.cos(ra)*np.sin(de),-np.sin(ra)*np.sin(de),np.cos(de)),axis=-1)
    line=np.stack((np.cos(ra)*np.cos(de),np.sin(ra)*np.cos(de),np.sin(de)),axis=-1)
    return np.stack((east,north,line),axis=-2)

def project(vectors,basis):
    x=np.asarray(vectors,float)
    if x.ndim==2:return np.einsum('nij,nj->ni',basis,x)
    if x.ndim==3:return np.einsum('nij,ntj->nti',basis,x)
    raise ValueError('expected (N,3) or (N,T,3)')

def propagate(position,velocity,q,total_mass,times_year,acceleration,*,substeps=1):
    """Kick-drift-kick at every epoch; force callback must be a potential gradient.

    The t=0 state precedes the first epoch. The callback takes (N,3), (N,) q.
    Substeps are deterministic and shared by every system at each time interval.
    Returns positions and velocities in dimensionless binary units.
    """
    x=np.asarray(position,float).copy();v=np.asarray(velocity,float).copy()
    q=np.asarray(q,float);t=np.asarray(times_year,float)
    if x.ndim!=2 or x.shape[1]!=3 or v.shape!=x.shape or q.shape!=(len(x),):raise ValueError('invalid phase space shapes')
    if not np.isfinite(x).all() or not np.isfinite(v).all() or not np.isfinite(q).all():raise ValueError('finite initial states required')
    if t.ndim!=1 or not np.isfinite(t).all() or np.any(np.diff(t)<0) or np.any(t<0) or int(substeps)!=substeps or substeps<1:raise ValueError('ordered nonnegative times and positive integer substeps required')
    tau=units(total_mass)['tM_year']
    if tau.shape!=(len(x),):raise ValueError('one mass per system required')
    paths=np.empty((len(x),len(t),3));velocities=np.empty_like(paths)
    acc=acceleration(x,q);previous=0.
    for it,time in enumerate(t):
        h=(time-previous)/tau/substeps
        for _ in range(substeps):
            half=v+.5*h[:,None]*acc
            x=x+h[:,None]*half
            acc=acceleration(x,q)
            v=half+.5*h[:,None]*acc
        paths[:,it]=x;velocities[:,it]=v;previous=time
    return paths,velocities

def invariant_barrier(energy,lz,cmax):
    """Lower turning radius for Psi >= -cmax/r; use only inside certified domain."""
    energy,lz=np.asarray(energy,float),np.asarray(lz,float)
    discriminant=cmax*cmax+2*energy*lz*lz
    allowed=np.isfinite(discriminant)&(discriminant>=0)&(energy<0)
    radius=np.zeros_like(energy)
    radius[allowed]=lz[allowed]**2/(cmax+np.sqrt(discriminant[allowed]))
    return radius,allowed
