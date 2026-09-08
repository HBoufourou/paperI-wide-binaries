"""Independent reconstruction of the encoded conservative interaction model.

No production imports. Cubics are solved as a full constraint system, rather
than a second-derivative spline system. Bernstein interval bounds are checked
against exact rational arithmetic on the stored IEEE nodes and coefficients.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[_key]='1'
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cubic_constraints(nodes,values):
    """Solve normalized local cubic coefficients using all interpolation constraints."""
    nodes=np.asarray(nodes,float); values=np.asarray(values,float); h=np.diff(nodes)
    n=len(nodes); m=n-1
    if n<4 or np.any(h<=0) or len(values)!=n:raise ValueError('Invalid cubic input')
    mat=np.zeros((4*m,4*m)); rhs=np.zeros((4*m,)+values.shape[1:]); row=0
    for j in range(m):
        mat[row,4*j]=1; rhs[row]=values[j];row+=1
        mat[row,4*j:4*j+4]=1;rhs[row]=values[j+1];row+=1
    for j in range(m-1):
        mat[row,4*j:4*j+4]=np.array([0,1,2,3])/h[j]
        mat[row,4*(j+1)+1]=-1/h[j+1];row+=1
        mat[row,4*j:4*j+4]=np.array([0,0,2,6])/h[j]**2
        mat[row,4*(j+1)+2]=-2/h[j+1]**2;row+=1
    mat[row,3]=1/h[0]**3;mat[row,7]=-1/h[1]**3;row+=1
    mat[row,4*(m-2)+3]=1/h[-2]**3;mat[row,4*(m-1)+3]=-1/h[-1]**3;row+=1
    assert row==4*m
    coeff=np.linalg.solve(mat,rhs.reshape(4*m,-1)).reshape((m,4)+values.shape[1:])
    scale=h[:,None]**np.arange(4)[None,:]
    return coeff/scale.reshape(scale.shape+(1,)*(values.ndim-1))


def tensor_constraints(x,y,values):
    cc=cubic_constraints(x,values)
    return cubic_constraints(y,np.moveaxis(cc,-1,0)).transpose(2,0,3,1)


class IndependentModel:
    def __init__(self,path):
        self.path=Path(path)
        with np.load(path,allow_pickle=False) as a:self.raw={k:a[k] for k in a.files}
        self.x=self.raw['log_r'];self.y=self.raw['mu_nodes'];self.q=self.raw['q_nodes']
        self.c=self.raw['coefficients']

    def evaluate(self,position,q,gravity='qumond'):
        p=np.asarray(position,float);qs=np.broadcast_to(q,p.shape[:-1])
        r=np.linalg.norm(p,axis=-1);mu=p[...,2]/r;x=np.log(r)
        if np.any(x<self.x[0]-1e-12) or np.any(x>self.x[-1]+1e-12):raise ValueError('Outside model radius')
        qi=np.argmin(abs(qs[...,None]-self.q),axis=-1)
        if np.any(abs(qs-self.q[qi])>1e-12):raise ValueError('Unsupported mass ratio')
        if gravity not in ('newton','qumond'):raise ValueError('Unsupported gravity')
        ii=np.searchsorted(self.x,x,side='right')-1; jj=np.searchsorted(self.y,mu,side='right')-1
        ii=np.clip(ii,0,len(self.x)-2);jj=np.clip(jj,0,len(self.y)-2)
        cc=self.c[int(gravity=='qumond'),qi,ii,jj];dx=x-self.x[ii];dy=mu-self.y[jj]
        # Horner along each axis, differentiated analytically.
        py=np.zeros(cc.shape[:-2]+(4,));dpy=np.zeros_like(py)
        for k in range(3,-1,-1):
            dpy=dpy*dy[...,None]+py;py=py*dy[...,None]+cc[...,k]
        aa=np.zeros_like(r);ax=np.zeros_like(r);ay=np.zeros_like(r)
        for k in range(3,-1,-1):
            ax=ax*dx+aa;aa=aa*dx+py[...,k];ay=ay*dx+dpy[...,k]
        # Build ∇Psi from spherical coordinate gradients, separately from the
        # algebraically combined expression used by the production evaluator.
        unit=p/r[...,None];grad_r=(aa-ax)/r**2;grad_mu=-ay/r
        dmu=-mu[...,None]*unit/r[...,None];dmu[...,2]+=1/r
        force=-(grad_r[...,None]*unit+grad_mu[...,None]*dmu)
        return -aa/r,force

    def potential(self,p,q,gravity='qumond'):return self.evaluate(p,q,gravity)[0]


def exact_bernstein(coeff,x,y,lower,upper):
    """Every output float interval must contain the exact real coefficient."""
    failures=[];count=0;maxwidth=0.
    ratio=[[F(math.comb(k,i),math.comb(3,i)) for i in range(k+1)] for k in range(4)]
    for ix in range(len(x)-1):
        hx=F(float(x[ix+1]))-F(float(x[ix]));xp=[hx**i for i in range(4)]
        for iy in range(len(y)-1):
            hy=F(float(y[iy+1]))-F(float(y[iy]));yp=[hy**j for j in range(4)]
            cc=[[F(float(coeff[ix,iy,i,j]))*xp[i]*yp[j] for j in range(4)] for i in range(4)]
            for k in range(4):
                for l in range(4):
                    exact=sum((cc[i][j]*ratio[k][i]*ratio[l][j] for i in range(k+1) for j in range(l+1)),F(0))
                    lo=float(lower[ix,iy,k,l]);hi=float(upper[ix,iy,k,l])
                    if not F(lo)<=exact<=F(hi):failures.append([ix,iy,k,l,lo,hi,float(exact)])
                    maxwidth=max(maxwidth,hi-lo);count+=1
    return count,failures,maxwidth


def fixtures():
    x=np.array([-.7,-.1,.2,1.1,1.5]);y=np.array([-1.,-.7,.1,.4,1.])
    vals=2+x[:,None]**3*(1+y[None,:]+y[None,:]**2)+x[:,None]*y[None,:]**3
    c=tensor_constraints(x,y,vals)
    rng=np.random.default_rng(268900001);err=0.
    for j in range(4):
        for k in range(4):
            dx=rng.uniform(0,x[j+1]-x[j],20);dy=rng.uniform(0,y[k+1]-y[k],20)
            p=np.polynomial.polynomial.polyval2d(dx,dy,c[j,k])
            xx=x[j]+dx;yy=y[k]+dy
            exact=2+xx**3*(1+yy+yy**2)+xx*yy**3
            err=max(err,float(np.max(abs(p-exact))))
    mock=IndependentModel.__new__(IndependentModel);mock.x=x;mock.y=y;mock.q=np.array([.1,.3,1.])
    mock.c=np.tile(c,(2,3,1,1,1,1))
    pos=rng.normal(size=(30,3));pos/=np.linalg.norm(pos,axis=1)[:,None]
    pos*=np.exp(rng.uniform(-.5,1.3,len(pos)))[:,None]
    def exact_potential(p):
        r=np.linalg.norm(p,axis=1);l=np.log(r);mu=p[:,2]/r
        return -(2+l**3*(1+mu+mu**2)+l*mu**3)/r
    pot,force=mock.evaluate(pos,.3)
    pe=float(np.max(abs(pot-exact_potential(pos))))
    h=1e-5;grad=np.empty_like(pos)
    for axis in range(3):
        d=np.eye(3)[axis]*h
        grad[:,axis]=(-exact_potential(pos+2*d)+8*exact_potential(pos+d)-8*exact_potential(pos-d)+exact_potential(pos-2*d))/(12*h)
    ge=float(np.max(abs(force+grad)))
    return {'bicubic_known_polynomial_absolute_error':err,
            'horner_evaluation_known_potential_absolute_error':pe,
            'gradient_closed_potential_difference_absolute_error':ge,
            'passed':bool(err<1e-11 and pe<1e-11 and ge<1e-7)}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--model',type=Path,default=Path(__file__).parents[1]/'physics/results/interaction_model.npz')
    parser.add_argument('--out',type=Path,default=Path(__file__).with_name('interpolant_verification.json'))
    parser.add_argument('--probes',type=Path)
    parser.add_argument('--fixture-only',action='store_true')
    args=parser.parse_args();start=time.time();out={'source_sha256':sha(__file__),'fixtures':fixtures(),'complete':False}
    if args.fixture_only:
        out['complete']=True;out['passed']=out['fixtures']['passed']
    elif not args.model.exists():
        out['status']='waiting_for_model';out['required_path']=str(args.model)
    else:
        model=IndependentModel(args.model);a=model.raw;x=model.x;y=model.y
        out['model_sha256']=sha(args.model);out['model_path']=str(args.model)
        out['models']=[];out['scope']='Encoded scalar polynomial reconstruction and bounds; no new physical accuracy guarantee.'
        for g,name in enumerate(('newton','qumond')):
            for qi,q in enumerate(model.q):
                cc=tensor_constraints(x,y,a['A_values'][g,qi]);target=model.c[g,qi]
                coefferr=float(np.max(abs(cc-target)))
                count,fails,width=exact_bernstein(target,x,y,a['bernstein_lower'][g,qi],a['bernstein_upper'][g,qi])
                bounds=(float(a['cmax'][g,qi])>=float(a['bernstein_upper'][g,qi].max()) and
                        float(a['cmin'][g,qi])<=float(a['bernstein_lower'][g,qi].min()) and
                        float(a['cmin'][g,qi])>0 and
                        float(a['outer_amax'][g,qi])>=float(a['bernstein_upper'][g,qi,-1,:,3,:].max()))
                row={'gravity':name,'q':float(q),'coefficient_max_absolute_difference':coefferr,
                     'exact_bernstein_intervals_checked':count,'exact_interval_failures':fails,
                     'max_interval_width':width,'global_and_outer_bounds_passed':bool(bounds),
                     'passed':bool(coefferr<1e-9 and not fails and bounds)}
                out['models'].append(row);print(json.dumps(row),flush=True)
        out['passed']=out['fixtures']['passed'] and all(r['passed'] for r in out['models'])
        if args.probes:
            with np.load(args.probes,allow_pickle=False) as z:probes={k:z[k] for k in z.files}
            metadata=json.loads(str(probes['metadata']));probe_rows=[]
            assert metadata['model_sha256']==out['model_sha256'],'Probe model hash mismatch'
            for case in metadata['cases']:
                k=case['key'];pos=probes[k+'_position'];qs=probes[k+'_q']
                shape=np.broadcast_shapes(pos.shape[:-1],qs.shape)
                pos=np.broadcast_to(pos,shape+(3,));qs=np.broadcast_to(qs,shape)
                pp,aa=model.evaluate(pos,qs,case['gravity'])
                targetp=probes[k+'_potential'];targeta=probes[k+'_acceleration']
                pe=float(np.max(abs(pp-targetp)));ae=float(np.max(abs(aa-targeta)))
                good=(pp.shape==targetp.shape and aa.shape==targeta.shape and
                      np.allclose(pp,targetp,rtol=3e-13,atol=3e-13) and np.allclose(aa,targeta,rtol=3e-12,atol=3e-12))
                probe_rows.append(dict(case,potential_max_absolute_error=pe,acceleration_max_absolute_error=ae,passed=bool(good)))
            out['api_probes']={'file_sha256':sha(args.probes),'metadata':metadata,'cases':probe_rows,
                              'passed':all(r['passed'] for r in probe_rows) and all(r['rejected'] for r in metadata['invalid_inputs'])}
            out['passed']=out['passed'] and out['api_probes']['passed']
        out['complete']=True
    out['seconds']=time.time()-start;args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps(out),flush=True)

if __name__=='__main__':main()
