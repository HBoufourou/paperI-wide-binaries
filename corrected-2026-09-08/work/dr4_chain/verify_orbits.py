"""Analytical Hamiltonian fixtures for the trajectory integrator, not QUMOND evidence."""
from pathlib import Path
import hashlib,json
import numpy as np
from orbits import propagate,units,invariant_barrier

HERE=Path(__file__).resolve().parent

def potential(x,anisotropic):
    r=np.linalg.norm(x,axis=1);mu=x[:,2]/r
    return -(1.4+.2*mu*mu if anisotropic else np.ones(len(x)))/r

def acceleration(x,q,anisotropic):
    r=np.linalg.norm(x,axis=1);h=x/r[:,None];mu=h[:,2]
    a=1.4+.2*mu*mu if anisotropic else np.ones(len(x))
    am=.4*mu if anisotropic else np.zeros(len(x))
    return -a[:,None]*h/r[:,None]**2+am[:,None]*(np.array([0,0,1])-mu[:,None]*h)/r[:,None]**2

def main():
    rng=np.random.default_rng(268000001);n=24
    theta=rng.uniform(0,2*np.pi,n);radius=rng.uniform(.8,1.5,n)
    x=np.column_stack((radius*np.cos(theta),radius*np.sin(theta),rng.uniform(-.2,.2,n)))
    r=np.linalg.norm(x,axis=1);v=np.column_stack((-np.sin(theta),np.cos(theta),rng.uniform(-.15,.15,n)))/np.sqrt(r[:,None])
    q=np.full(n,.3);mass=np.ones(n);times=np.linspace(0,12,1201)*float(units(1.)['tM_year'])
    records=[]
    def check(name,value,threshold):
        row=dict(name=name,value=float(value),threshold=threshold,passed=bool(value<threshold));records.append(row)
        if not row['passed']:raise AssertionError(row)
    for anisotropic in (False,True):
        label='homogeneous_anisotropic' if anisotropic else 'Newton_point'
        acc=lambda p,qq:acceleration(p,qq,anisotropic)
        p1,v1=propagate(x,v,q,mass,times,acc,substeps=1)
        p2,v2=propagate(x,v,q,mass,times,acc,substeps=2)
        initial=.5*np.sum(v*v,axis=1)+potential(x,anisotropic)
        energy1=.5*np.sum(v1*v1,axis=2)+potential(p1.reshape(-1,3),anisotropic).reshape(n,-1)
        energy2=.5*np.sum(v2*v2,axis=2)+potential(p2.reshape(-1,3),anisotropic).reshape(n,-1)
        err1=float(np.max(abs((energy1-initial[:,None])/initial[:,None])))
        err2=float(np.max(abs((energy2-initial[:,None])/initial[:,None])))
        check(label+'/relative_energy',err2,1e-4)
        check(label+'/refinement_ratio',err2/err1,.35)
        lz=np.cross(p2,v2)[:,:,2];lz0=np.cross(x,v)[:,2]
        check(label+'/Lz_conservation',np.max(abs(lz-lz0[:,None])),1e-10)
        reverse_x,reverse_v=propagate(p1[:,-1],-v1[:,-1],q,mass,times,acc,substeps=1)
        check(label+'/reversibility_position',np.max(abs(reverse_x[:,-1]-x)),1e-9)
        check(label+'/reversibility_velocity',np.max(abs(reverse_v[:,-1]+v)),1e-9)
    energy=np.array([-.2,-.5,-1.]);lz=np.array([.4,.7,.8]);cmax=2.
    lower,valid=invariant_barrier(energy,lz,cmax)
    check('barrier_root',np.max(abs(lz**2/(2*lower**2)-cmax/lower-energy)),1e-10)
    assert valid.all()
    result=dict(passed=True,checks=len(records),records=records,
        scope='Analytical Newton point and homogeneous anisotropic Hamiltonians; does not validate two-mass QUMOND interpolant.',
        physical_model_validated=False,sky_ready=False,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'orbits.py',Path(__file__)]})
    (HERE/'ORBIT_IMPLEMENTATION_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=True,checks=len(records))))

if __name__=='__main__':main()
