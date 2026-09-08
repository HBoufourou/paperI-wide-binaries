"""Independent recount of the final pointwise physical-validation evidence."""
import argparse,json
from pathlib import Path
import numpy as np
from verify_interpolant import IndependentModel,sha


def rel(a,b):return float(np.linalg.norm(np.asarray(a)-b)/np.linalg.norm(b))


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).parents[1])
    p.add_argument('--out',type=Path,default=Path(__file__).with_name('physics_report_verification.json'))
    args=p.parse_args();root=args.root.resolve();path=root/'physics/validation.json'
    out={'source_sha256':sha(__file__),'independent_model_source_sha256':sha(Path(__file__).with_name('verify_interpolant.py')),
         'complete':False,'scope':'Recounts saved pointwise evidence; direct force integrations are not rerun.'}
    if not path.exists():out['status']='waiting_for_final_physics_report'
    else:
        report=json.loads(path.read_text());model=IndependentModel(root/'physics/results/interaction_model.npz')
        checks=[];metrics=[];hashes={str(path.relative_to(root)):sha(path)}
        def check(name,condition,error=None):checks.append({'name':name,'passed':bool(condition),'absolute_error':error})
        def close(name,a,b,atol=3e-9,rtol=1e-5):
            error=float(np.max(abs(np.asarray(a)-b)));check(name,np.allclose(a,b,atol=atol,rtol=rtol),error)
        check('model_hash',sha(model.path)==report['model_sha256'])
        for row in report['rows']:
            q,r,mu,g=row['q'],row['r'],row['mu'],row['gravity'];pos=r*np.array([np.sqrt(1-mu*mu),0.,mu])
            fpath=root/'physics/results/forces'/f'final_q{q:g}_r{r:g}_mu{mu:g}.json'
            raw=json.loads(fpath.read_text());hashes[str(fpath.relative_to(root))]=sha(fpath)
            if g=='qumond':reference=np.asarray(raw['relative_acceleration'])
            else:
                f=np.asarray(raw['newton_force']);mass=np.asarray(raw['masses']);reference=f[1]/mass[1]-f[0]/mass[0]
            acc=model.evaluate(pos,q,g)[1];error=rel(acc,reference)
            label=f'{g}/q{q}/r{r}/mu{mu}'
            close(label+'/reference_acceleration',row['reference_acceleration'],reference,atol=1e-12,rtol=1e-12)
            close(label+'/acceleration',row['interpolated_acceleration'],acc,atol=1e-10,rtol=1e-12)
            close(label+'/relative_error',row['relative_force_error'],error,atol=1e-11,rtol=1e-9)
            h=1e-5*r
            fd=np.array([-(model.potential(pos+h*axis,q,g)-model.potential(pos-h*axis,q,g))/(2*h) for axis in np.eye(3)])
            grad=rel(fd,acc);h=1e-4*r
            jac=np.column_stack([(model.evaluate(pos+h*axis,q,g)[1]-model.evaluate(pos-h*axis,q,g)[1])/(2*h) for axis in np.eye(3)])
            curl=float(np.linalg.norm(jac-jac.T)/np.linalg.norm(jac))
            close(label+'/finite_difference',row['finite_difference_gradient_error'],grad)
            close(label+'/jacobian',row['jacobian_asymmetry'],curl)
            metrics.append({'q':q,'r':r,'mu':mu,'gravity':g,'force_error':error,'gradient_error':grad,'jacobian_asymmetry':curl,
                            'criteria_passed':error<.01 and grad<1e-5 and curl<1e-4})
        energy=[]
        for row in report['energy_controls']:
            q,r=row['q'],row['r'];ee=[]
            for name in ('base','reference','outer'):
                ep=root/'physics/results/energy_controls'/f'q{q:g}_r{r:g}_{name}.json'
                rr=json.loads(ep.read_text());hashes[str(ep.relative_to(root))]=sha(ep);ee.append(rr['Psi'])
            mesh=abs(ee[0]/ee[1]-1);outer=abs(ee[0]/ee[2]-1)
            close(f'energy/q{q}/r{r}/mesh',row['mesh_relative'],mesh,atol=1e-14,rtol=1e-12)
            close(f'energy/q{q}/r{r}/outer',row['outer_relative'],outer,atol=1e-14,rtol=1e-12)
            energy.append({'q':q,'r':r,'mesh_relative':mesh,'outer_relative':outer,'criteria_passed':mesh<.001 and outer<.001})
        old=[]
        for row in report['old_reference_comparisons']:
            q=row['q'];fp=root.parent/'dr4_binary/physics/results/jobs'/f'reference_E1_q{q:g}_d1_t45_reference.json'
            raw=json.loads(fp.read_text());hashes[str(fp.relative_to(root.parent))]=sha(fp)
            acc=model.evaluate(np.array([2**-.5,0.,2**-.5]),q)[1];e=rel(acc,np.array(raw['relative_acceleration']))
            close(f'old_reference/q{q}',row['relative_error'],e,atol=1e-12,rtol=1e-8);old.append({'q':q,'relative_error':e})
        expected_points={(q,r,mu,g) for q in (.1,.3,1.) for r,mu in ((.045,.61),(.12,-.22),(.41,.87),(.9,-.71),(2.8,.19),(11.5,-.46),(24,.79)) for g in ('newton','qumond')}
        check('complete_final_point_set',len(metrics)==42 and {(r['q'],r['r'],r['mu'],r['gravity']) for r in metrics}==expected_points)
        science=all(r['criteria_passed'] for r in metrics+energy)
        check('reported_final_status',science==report['passed'])
        close('reported_max_force_error',report['max_force_error'],max(r['force_error'] for r in metrics),atol=1e-12,rtol=1e-9)
        out.update(complete=True,report_reconstruction_passed=all(c['passed'] for c in checks),physics_criteria_passed=science,
                   checks=checks,checks_count=len(checks),force_metrics=metrics,energy_metrics=energy,old_reference_metrics=old,
                   input_sha256=hashes,model_sha256=sha(model.path))
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('checks','force_metrics','energy_metrics','old_reference_metrics','input_sha256')}),flush=True)

if __name__=='__main__':main()
