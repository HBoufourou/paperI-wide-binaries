"""Independent potential differentiation against frozen earlier force references."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
from independent_energy import potential, qphantom, mixed_action, newton_energy


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def analytic():
    y=np.geomspace(1e-10,1e8,100); h=y*1e-3
    fd=(-qphantom(y+2*h)+8*qphantom(y+h)-8*qphantom(y-h)+qphantom(y-2*h))/(12*h)
    expected=4*y/(y+np.sqrt(y)*np.sqrt(y+4))
    derivative=float(np.max(np.abs(fd/expected-1)))
    rng=np.random.default_rng(275000001)
    ue=np.array([0.,0.,-1.]); a=rng.normal(size=(100,3)); b=rng.normal(size=(100,3))
    a*=1e-3/np.linalg.norm(a,axis=1)[:,None]; b*=1e-3/np.linalg.norm(b,axis=1)[:,None]
    direct=mixed_action(a,b,ue,0.)
    stable=mixed_action(a,b,ue,.1)
    mix=float(np.max(np.abs(direct-stable))/np.max(np.abs(stable)))
    # Leading mixed-Hessian action under a truly weak perturbation.
    eps=1e-5; vv=mixed_action(a*eps,b*eps,ue,.002)/eps**2
    nu=(1+np.sqrt(5.))/2; nup=-1/np.sqrt(5.)
    leading=2*((nu-1)*np.sum(a*b,axis=1)+nup*a[:,2]*b[:,2])
    tail=float(np.max(np.abs(vv-leading))/np.max(np.abs(leading)))
    return {'qphantom_derivative_relative_error':derivative,
            'mixed_integral_vs_direct_at_002_absolute_normalized':mix,
            'external_tail_hessian_limit_relative_error':tail,
            'passed':bool(derivative<1e-7 and mix<1e-7 and tail<1e-7)}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=Path(__file__).with_name('energy_validation.json'))
    parser.add_argument('--force-root',type=Path,default=Path(__file__).parents[2]/'dr4_binary/physics/results/jobs')
    args=parser.parse_args(); start=time.time()
    source=Path(__file__).with_name('independent_energy.py')
    out={'implementation_sha256':sha(source),'driver_sha256':sha(Path(__file__)),
         'analytic':analytic(),'integrals':[],'gradients':[],
         'scope':'Independent prolate potential: finite Plummer q=1,.3; E_N=1. No population/sky conclusion.',
         'predeclared_tolerances':{'last_order_force_relative':.0001,'last_step_force_relative':.0001,
                                 'force_reference_relative':.0001},'complete':False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    def save(): args.out.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    def run(q,theta,order,d=1.,offset_theta=0.):
        t=time.time(); v=potential(q,d,theta+offset_theta,order=order)
        v['seconds']=time.time()-t; out['integrals'].append(v); save()
        print(json.dumps({'integral':len(out['integrals']),'q':q,'theta':theta+offset_theta,
                          'd':d,'order':order,'U':v['U'],'seconds':v['seconds']}),flush=True)
        return v['U']
    for q,theta in [(1.,45.),(.3,45.),(.3,135.)]:
        centre=run(q,theta,12)
        ref=None
        if theta==45.:
            name='reference_E1_q'+('1' if q==1 else '0.3')+'_d1_t45_reference.json'
            path=args.force_root/name
            payload=json.loads(path.read_text())
            # Reference force field is inspected explicitly below, never imported.
            force=np.asarray(payload['force'])[0] if 'force' in payload else np.asarray(payload['forces'])[0]
            ref={'path':str(path),'sha256':sha(path),'force':force.tolist()}
        for order,step in ([(8,.004),(12,.004),(12,.002)] if theta==45. else [(12,.004)]):
            um=run(q,theta,order,d=1-step); up=run(q,theta,order,d=1+step)
            tm=run(q,theta,order,offset_theta=-math.degrees(step))
            tp=run(q,theta,order,offset_theta=math.degrees(step))
            dr=(up-um)/(2*step); dt=(tp-tm)/(2*step)
            t=math.radians(theta); rh=np.array([math.sin(t),0.,math.cos(t)])
            th=np.array([math.cos(t),0.,-math.sin(t)])
            f=dr*rh+dt*th
            row={'q':q,'theta':theta,'order':order,'step_r':step,'step_theta_radians':step,
                 'U_centre':centre,'radial_derivative':dr,'angular_derivative':dt,
                 'force1_lab':f.tolist(),'reference':ref}
            if ref:
                ff=np.asarray(ref['force'])
                row['force_reference_relative_error']=float(np.linalg.norm(f-ff)/np.linalg.norm(ff))
                row['force_reference_component_difference']=(f-ff).tolist()
            out['gradients'].append(row); save()
            print(json.dumps({'gradient':row}),flush=True)
    errors=[]; out['convergence']=[]
    for q in [1.,.3]:
        rows=[v for v in out['gradients'] if v['q']==q and v['theta']==45.]
        f=[np.asarray(r['force1_lab']) for r in rows]
        ordererr=float(np.linalg.norm(f[1]-f[0])/np.linalg.norm(f[1]))
        steperr=float(np.linalg.norm(f[2]-f[1])/np.linalg.norm(f[2]))
        out['convergence'].append({'q':q,'order_relative_change':ordererr,
                                   'step_relative_change':steperr})
        errors.extend([ordererr,steperr,rows[-1]['force_reference_relative_error']])
    out['numerical_checks_passed']=bool(out['analytic']['passed'] and max(errors)<.0001)
    out['seconds']=time.time()-start; out['complete']=True; save()
    print(json.dumps({'complete':True,'passed':out['numerical_checks_passed'],
                      'seconds':out['seconds'],'max_checked_error':max(errors)}),flush=True)


if __name__=='__main__': main()
