"""Independent QUMOND interaction action in prolate coordinates.

No production imports. G=Mtot=a0=1; finite Plummer sources, uniform Newtonian
external field. Newtonian energy is integrated independently in one dimension;
the phantom action is integrated to infinity, without a finite-radius tail
correction. All angular arguments are signed through theta in [0,180] degrees.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np


def gauss(edges, order):
    z,w=np.polynomial.legendre.leggauss(order)
    e=np.unique(np.asarray(edges,float)); lo=e[:-1,None]; hi=e[1:,None]
    return ((lo+hi+(hi-lo)*z)/2).ravel(), ((hi-lo)*w/2).ravel()


def newton_energy(m1,m2,b1,b2,d):
    edges=np.r_[0.,np.geomspace(b1/32,d*100,80),
                d-10*b2,d-b2,d,d+b2,d+10*b2]
    edges=np.unique(edges[edges>=0]); r,w=gauss(edges,64)
    t,tw=gauss([0.,1.],64); end=edges[-1]
    r=np.r_[r,end/(1-t)]; w=np.r_[w,tw*end/(1-t)**2]
    shell=3*b1*b1*r*r/(r*r+b1*b1)**2.5
    aa=np.sqrt((r+d)**2+b2*b2); bb=np.sqrt((r-d)**2+b2*b2)
    return float(-m1*m2*np.sum(w*shell*2/(aa+bb)))


def qphantom(y):
    """Q(y^2)-y^2, for the simple nu function; stable near 0 and infinity."""
    y=np.asarray(y,float); s=np.sqrt(y)*np.sqrt(y+4)
    den=s+y
    z=s+np.divide(2*y*y,den,out=np.zeros_like(y),where=den>0)-4*np.arcsinh(np.sqrt(y)/2)
    small=y<1e-3
    v=y[small]
    z[small]=v**1.5*(4/3+v/10-v*v/224+v**3/2304-5*v**4/90112)-v*v/2
    return z


def mixed_action(u1,u2,ue,threshold=.002):
    """Four-term Qph action; weak-field increments use its mixed Hessian.

    The 2x2 Gauss formula integrates the exact double-integral representation,
    avoiding subtraction of four almost equal numbers in the infinite tail.
    Its threshold is a declared numerical parameter, not a force-field model.
    """
    n1=np.linalg.norm(u1,axis=-1); n2=np.linalg.norm(u2,axis=-1)
    weak=n1+n2<threshold*np.linalg.norm(ue)
    out=np.empty(n1.shape)
    v1=u1[~weak]; v2=u2[~weak]
    out[~weak]=(qphantom(np.linalg.norm(ue+v1+v2,axis=-1))-
                qphantom(np.linalg.norm(ue+v1,axis=-1))-
                qphantom(np.linalg.norm(ue+v2,axis=-1))+
                float(qphantom(np.array([np.linalg.norm(ue)]))[0]))
    if np.any(weak):
        v1=u1[weak]; v2=u2[weak]; dot=np.sum(v1*v2,axis=-1)
        ans=np.zeros(len(v1)); ss=(.5-.5/math.sqrt(3),.5+.5/math.sqrt(3))
        for a in ss:
            for b in ss:
                v=ue+a*v1+b*v2; y=np.linalg.norm(v,axis=-1)
                extra=2/(y+np.sqrt(y)*np.sqrt(y+4))
                derivative=-1/(y**1.5*np.sqrt(y+4))
                ans+=.5*(extra*dot+derivative/y*np.sum(v1*v,axis=-1)*np.sum(v2*v,axis=-1))
        out[weak]=ans
    return out


def potential(q=.3,d=1.,theta_deg=45.,ext=1.,epsilon=.00125,
              order=8,nphi=None,radial_max=1000.,threshold=.002):
    if ext<=0:
        raise ValueError('This zero-at-infinity interaction action requires ext>0')
    nphi=nphi or 4*order
    m1=1/(1+q); m2=q/(1+q); b1=epsilon*np.sqrt(m1); b2=epsilon*np.sqrt(m2)
    a=d/2; floor=min(b1,b2)/a/32
    tx,wx=gauss(np.r_[0.,np.geomspace(floor,radial_max-1,44)],order)
    xi=1+tx; t,tw=gauss([0.,1.],order)
    xi=np.r_[xi,radial_max/(1-t)]; wx=np.r_[wx,tw*radial_max/(1-t)**2]
    pole=np.r_[0.,np.geomspace(floor,1.,32)]
    eta,ew=gauss(np.r_[-1+pole,1-pole],order)
    phi=(np.arange(nphi)+.5)*2*np.pi/nphi
    cp=np.cos(phi)[None,:]; sp=np.sin(phi)[None,:]; ee=eta[:,None]
    angle=math.radians(theta_deg); ue=ext*np.array([math.sin(angle),0.,-math.cos(angle)])
    integral=0.
    for x,w in zip(xi,wx):
        cyl=a*np.sqrt((x*x-1)*(1-ee*ee))
        z=np.broadcast_to(a*x*ee,(len(eta),nphi))
        r1=np.stack((cyl*cp,cyl*sp,z+a),axis=-1)
        r2=np.stack((cyl*cp,cyl*sp,z-a),axis=-1)
        u1=m1*r1/(np.sum(r1*r1,axis=-1)+b1*b1)[...,None]**1.5
        u2=m2*r2/(np.sum(r2*r2,axis=-1)+b2*b2)[...,None]**1.5
        action=mixed_action(u1,u2,ue,threshold)
        integral+=np.sum(action*w*ew[:,None]*(2*np.pi/nphi)*a**3*(x*x-ee*ee))
    un=newton_energy(m1,m2,b1,b2,d); up=-integral/(8*np.pi)
    return {'q':q,'d':d,'theta_deg':theta_deg,'ext_newtonian':ext,
            'epsilon':epsilon,'order':order,'nphi':nphi,'radial_max':radial_max,
            'threshold':threshold,'masses':[m1,m2],'softenings':[b1,b2],
            'newton_energy':un,'phantom_energy':float(up),'U':float(un+up),
            'Psi':float((un+up)/(m1*m2)),
            'nodes':len(xi)*len(eta)*nphi}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--q',type=float,default=.3)
    parser.add_argument('--d',type=float,default=1.)
    parser.add_argument('--theta',type=float,default=45.)
    parser.add_argument('--orders',type=int,nargs='+',default=[6,8])
    parser.add_argument('--threshold',type=float,default=.002)
    parser.add_argument('--radial-max',type=float,default=1000.)
    parser.add_argument('--out',type=Path,default=Path(__file__).with_name('energy_pilot.json'))
    args=parser.parse_args()
    out={'implementation':'independent prolate action; infinite rational tail; no production imports',
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'results':[]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    for order in args.orders:
        start=time.time()
        row=potential(args.q,args.d,args.theta,order=order,threshold=args.threshold,
                      radial_max=args.radial_max)
        row['seconds']=time.time()-start; out['results'].append(row)
        args.out.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(row),flush=True)


if __name__=='__main__': main()
