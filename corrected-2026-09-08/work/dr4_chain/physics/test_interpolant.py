"""Algebraic checks independent of physical energy-grid results."""
from pathlib import Path
import tempfile,json
import numpy as np
from interaction_model import tensor_coefficients,bernstein_bounds,write_model,InteractionModel

def run():
    x=np.linspace(np.log(.03),np.log(30),17);y=np.linspace(-1,1,9)
    xx,yy=np.meshgrid(x,y,indexing='ij')
    a=2+.01*xx**3+.04*xx*yy+.02*yy**3
    c=tensor_coefficients(x,y,a);b,bl,bh=bernstein_bounds(c,x,y)
    errors=[];contain=True
    for i in range(len(x)-1):
        for j in range(len(y)-1):
            for t,s in [(.13,.87),(.5,.5),(.91,.07),(0,0),(1,1)]:
                dx=t*(x[i+1]-x[i]);dy=s*(y[j+1]-y[j]);xp=x[i]+dx;yp=y[j]+dy
                val=sum(c[i,j,k,l]*dx**k*dy**l for k in range(4) for l in range(4))
                truth=2+.01*xp**3+.04*xp*yp+.02*yp**3
                errors.append(abs(val-truth));contain &= bl[i,j].min()<=val<=bh[i,j].max()
    values=np.stack([np.stack([np.ones_like(a)]*3),np.stack([1.6*(1-.14*(1-yy*yy))]*3)])
    with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as td:
        p=Path(td)/'toy.npz';write_model(p,x,y,np.array([.1,.3,1]),values,{'test':True})
        model=InteractionModel.load(p,require_validation=False)
        pos=np.array([[.03,0,0],[.4,.3,.2],[0,0,30],[-.3,.1,-.5]])
        r=np.linalg.norm(pos,axis=1);mu=pos[:,2]/r
        truth=-1.6*(1-.14*(1-mu*mu))/r
        efe_error=float(np.max(np.abs(model.psi(pos,.3)-truth)))
        newton_error=float(np.max(np.abs(model.acceleration(pos,.1,'newton')+pos/r[:,None]**3)))
        refused=[]
        for p0,q0 in [([.029,0,0],.3),([30.01,0,0],1),([1,0,0],.2)]:
            try:model.psi(p0,q0);refused.append(False)
            except ValueError:refused.append(True)
    result=dict(cubic_max_error=max(errors),bernstein_contains_samples=bool(contain),
                efe_exact_max_error=efe_error,newton_exact_max_error=newton_error,domain_refusal=refused)
    result['passed']=bool(max(errors)<1e-12 and contain and efe_error<1e-12 and newton_error<1e-11 and all(refused))
    print(json.dumps(result,indent=2));return result

if __name__=='__main__':
    result=run();(Path(__file__).resolve().parent/'interpolant_unit_checks.json').write_text(json.dumps(result,indent=2))
    assert result['passed']
