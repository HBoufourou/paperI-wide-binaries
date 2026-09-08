"""Conservative scalar interaction interpolant, finite Plummer sources, E_N=1.

The supported potential is the encoded piecewise polynomial, not a certified
continuous solution of QUMOND between independently checked locations.
"""
from pathlib import Path
import json
import math
import hashlib
import numpy as np

HERE=Path(__file__).resolve().parent
DEFAULT_PATH=HERE/'results'/'interaction_model.npz'


def spline_coefficients(nodes, values):
    """Not-a-knot cubic; axis0 is node axis; unscaled local power basis."""
    x=np.asarray(nodes,float); y=np.asarray(values,float); h=np.diff(x); n=len(x)
    if n<4 or len(y)!=n or np.any(h<=0):raise ValueError('need >=4 ordered nodes')
    matrix=np.zeros((n,n)); rhs=np.zeros_like(y)
    matrix[0,:3]=[-h[1],h[0]+h[1],-h[0]]
    matrix[-1,-3:]=[-h[-1],h[-2]+h[-1],-h[-2]]
    shape=(n-2,)+(1,)*(y.ndim-1)
    hp=h[:-1].reshape(shape); hn=h[1:].reshape(shape)
    rhs[1:-1]=6*((y[2:]-y[1:-1])/hn-(y[1:-1]-y[:-2])/hp)
    for i in range(1,n-1):matrix[i,i-1:i+2]=[h[i-1],2*(h[i-1]+h[i]),h[i]]
    m=np.linalg.solve(matrix,rhs.reshape(n,-1)).reshape(y.shape)
    hs=h.reshape((n-1,)+(1,)*(y.ndim-1))
    return np.stack([y[:-1],np.diff(y,axis=0)/hs-hs*(2*m[:-1]+m[1:])/6,
                     m[:-1]/2,np.diff(m,axis=0)/(6*hs)],axis=1)


def tensor_coefficients(x,y,values):
    cx=spline_coefficients(x,values)
    return spline_coefficients(y,cx.transpose(2,0,1)).transpose(2,0,3,1)


def bernstein_bounds(coefficients,x,y):
    """Outward-rounded interval power-to-Bernstein transform.

    Float coefficients and nodes define the real polynomial. Every arithmetic
    operation in its interval transform is rounded outwards with nextafter.
    This is a bound for that polynomial under ordinary IEEE binary64 arithmetic.
    """
    c=np.asarray(coefficients,float); hx=np.diff(x)[:,None]; hy=np.diff(y)[None,:]
    low=lambda a:np.nextafter(a,-np.inf)
    high=lambda a:np.nextafter(a,np.inf)
    def powers(h):
        lo=[np.ones_like(h)]; hi=[np.ones_like(h)]
        hl=low(h);hh=high(h)
        for i in range(3):lo.append(low(lo[-1]*hl));hi.append(high(hi[-1]*hh))
        return lo,hi
    xl,xh=powers(hx);yl,yh=powers(hy)
    b=np.zeros_like(c); bl=np.zeros_like(c); bh=np.zeros_like(c)
    for k in range(4):
        for l in range(4):
            lower=np.zeros(c.shape[:-2]);upper=np.zeros(c.shape[:-2]);nom=np.zeros_like(lower)
            for i in range(k+1):
                for j in range(l+1):
                    tx=math.comb(k,i)/math.comb(3,i);ty=math.comb(l,j)/math.comb(3,j)
                    fl=low(low(low(xl[i]*yl[j])*low(tx))*low(ty))
                    fh=high(high(high(xh[i]*yh[j])*high(tx))*high(ty))
                    cc=c[...,i,j]
                    tl=low(cc*np.where(cc>=0,fl,fh));th=high(cc*np.where(cc>=0,fh,fl))
                    lower=low(lower+tl);upper=high(upper+th)
                    nom+=cc*hx**i*hy**j*tx*ty
            b[...,k,l]=nom;bl[...,k,l]=lower;bh[...,k,l]=upper
    return b,bl,bh


def write_model(path,log_r,mu_nodes,q_nodes,avalues,metadata):
    coefficients=np.empty((2,len(q_nodes),len(log_r)-1,len(mu_nodes)-1,4,4))
    bern=np.empty_like(coefficients);lower=np.empty_like(coefficients);upper=np.empty_like(coefficients)
    for g in range(2):
        for iq in range(len(q_nodes)):
            coefficients[g,iq]=tensor_coefficients(log_r,mu_nodes,avalues[g,iq])
            bern[g,iq],lower[g,iq],upper[g,iq]=bernstein_bounds(coefficients[g,iq],log_r,mu_nodes)
    cmax=np.max(upper,axis=(2,3,4,5));cmin=np.min(lower,axis=(2,3,4,5))
    outer_amax=np.max(upper[:,:,-1,:,3,:],axis=(2,3))
    if np.any(cmin<=0):raise ValueError('Positive A is not certified by Bernstein lower bounds')
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,log_r=log_r,mu_nodes=mu_nodes,q_nodes=q_nodes,
        gravity_names=np.array(['newton','qumond']),A_values=avalues,
        coefficients=coefficients,bernstein_coefficients=bern,
        bernstein_lower=lower,bernstein_upper=upper,cmax=cmax,cmin=cmin,outer_amax=outer_amax,
        metadata=json.dumps(metadata,sort_keys=True))


class InteractionModel:
    domain={'r_min':.03,'r_max':30.,'q':[.1,.3,1.], 'external_newtonian':1.,
            'epsilon':.00125,'units':'G=Mtotal=a0=1'}

    @classmethod
    def load(cls,path=DEFAULT_PATH,require_validation=True):
        path=Path(path)
        if require_validation:
            validation_path=path.parent.parent/'validation.json'
            if not validation_path.exists():raise RuntimeError('Final physics validation absent')
            validation=json.loads(validation_path.read_text(encoding='utf-8'))
            if not validation.get('passed',False):raise RuntimeError('Final physics validation failed')
            sha=hashlib.sha256(path.read_bytes()).hexdigest()
            if validation.get('model_sha256')!=sha:raise RuntimeError('Model differs from validated file')
            if validation.get('evaluator_sha256')!=hashlib.sha256(Path(__file__).read_bytes()).hexdigest():
                raise RuntimeError('Evaluator differs from validated source')
        obj=cls()
        with np.load(path,allow_pickle=False) as data:
            obj.x=data['log_r'];obj.y=data['mu_nodes'];obj.q_nodes=data['q_nodes']
            obj.coefficients=data['coefficients'];obj.bounds=data['cmax'];obj.lower_bounds=data['cmin']
            obj.outer_bounds=data['outer_amax'];obj.metadata=json.loads(str(data['metadata']))
        obj.path=path;return obj

    def _indices(self,q,gravity):
        if gravity not in ('newton','qumond'):raise ValueError('gravity must be newton or qumond')
        q=np.asarray(q,float);dist=np.abs(q[...,None]-self.q_nodes)
        iq=np.argmin(dist,axis=-1)
        if not np.all(np.min(dist,axis=-1)<=1e-12):raise ValueError('Only discrete q=.1,.3,1 supported')
        return int(gravity=='qumond'),iq

    def cmax(self,q,gravity='qumond'):
        g,iq=self._indices(q,gravity);return self.bounds[g,iq]

    def cmin(self,q,gravity='qumond'):
        g,iq=self._indices(q,gravity);return self.lower_bounds[g,iq]

    def outer_amax(self,q,gravity='qumond'):
        g,iq=self._indices(q,gravity);return self.outer_bounds[g,iq]

    def _evaluate(self,positions,q,gravity):
        p=np.asarray(positions,float)
        if p.shape[-1:]!=(3,) or not np.all(np.isfinite(p)):raise ValueError('Finite positions (...,3) required')
        shape=np.broadcast_shapes(p.shape[:-1],np.shape(q));p=np.broadcast_to(p,shape+(3,))
        g,iq=self._indices(np.broadcast_to(q,shape),gravity)
        r=np.linalg.norm(p,axis=-1)
        if np.any(r<.03*(1-1e-12)) or np.any(r>30*(1+1e-12)):raise ValueError('Radius outside [.03,30]')
        rr=np.clip(r,.03,30);x=np.clip(np.log(rr),self.x[0],self.x[-1]);mu=np.clip(p[...,2]/r,-1,1)
        i=np.clip(np.searchsorted(self.x,x,side='right')-1,0,len(self.x)-2)
        j=np.clip(np.searchsorted(self.y,mu,side='right')-1,0,len(self.y)-2)
        dx=x-self.x[i];dy=mu-self.y[j];c=self.coefficients[g,iq,i,j]
        a=np.zeros(shape);ax=np.zeros(shape);ay=np.zeros(shape)
        xp=[np.ones(shape),dx,dx*dx,dx*dx*dx];yp=[np.ones(shape),dy,dy*dy,dy*dy*dy]
        for k in range(4):
            for l in range(4):
                a+=c[...,k,l]*xp[k]*yp[l]
                if k:ax+=k*c[...,k,l]*xp[k-1]*yp[l]
                if l:ay+=l*c[...,k,l]*xp[k]*yp[l-1]
        return a,ax,ay,rr,mu,p/r[...,None]

    def potential(self,positions,q,gravity='qumond'):
        a,_,_,r,_,_=self._evaluate(positions,q,gravity);return -a/r

    psi=potential

    def acceleration(self,positions,q,gravity='qumond'):
        a,ax,ay,r,mu,rhat=self._evaluate(positions,q,gravity)
        out=-(a-ax+mu*ay)[...,None]*rhat/(r*r)[...,None]
        out[...,2]+=ay/(r*r)
        return out
