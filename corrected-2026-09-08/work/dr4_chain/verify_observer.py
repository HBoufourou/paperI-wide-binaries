"""Bounded joint-observer checks against a frozen legacy and analytic signals."""
from pathlib import Path
import hashlib,importlib.util,json,sys
import numpy as np
import joint_observer as new

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('legacy_chain_check',HERE.parent/'dr4_binary/companions/epoch_operator_extended.py')
old=importlib.util.module_from_spec(spec);sys.modules[spec.name]=old;spec.loader.exec_module(old)

def main():
    rng=np.random.default_rng(269000001);n=48
    distance=rng.uniform(40,180,n);errors=rng.uniform(.04,.12,(n,2,2))
    ra=rng.uniform(0,360,n);de=rng.uniform(-75,75,n)
    fraction=np.column_stack((rng.uniform(.5,.91,n),np.empty(n)));fraction[:,1]=1-fraction[:,0]
    inner=dict(active=rng.random(n)<.5,host=rng.integers(0,2,n),q=rng.uniform(.1,.9,n),
        mass_host_true=rng.uniform(.5,2,n),a_au=rng.uniform(.3,10,n),e=rng.uniform(0,.7,n),
        M0=rng.uniform(0,2*np.pi,n),rotation=old.random_rotations(n,rng))
    records=[]
    def close(name,a,b,atol=1e-11,rtol=1e-11):
        a,b=np.asarray(a),np.asarray(b)
        passed=bool(np.allclose(a,b,atol=atol,rtol=rtol))
        records.append(dict(name=name,passed=passed,max_abs=float(np.max(np.abs(a-b))),atol=atol,rtol=rtol))
        if not passed:raise AssertionError(records[-1])
    for si,scenario in enumerate(new.SCENARIOS):
        t,angle,n3=new.cadence(scenario)
        zero=np.zeros((n,len(t),2))
        args=(None,distance,errors,inner,scenario,269000002+si,ra,de)
        ref=old.observe_pairs(*args)
        null=new.observe_pairs(*args,outer_relative_xy_au=zero,mass_fraction=fraction)
        for key in ref:close(scenario+'/legacy_zero/'+key,null[key],ref[key])
        velocity=rng.uniform(-1,1,(n,2))
        xy=velocity[:,None,:]/new.K*t[None,:,None]
        straight=new.observe_pairs(*args,outer_relative_xy_au=xy,mass_fraction=fraction)
        legacy=old.observe_pairs(velocity,distance,errors,inner,scenario,269000002+si,ra,de)
        for b in (34,66):
            expected=velocity*1000/(new.K*distance[:,None])
            close(scenario+f'/straight_relative{b}',straight[f'pm{b}_masyr'],legacy[f'pm{b}_masyr'])
            close(scenario+f'/straight_outer{b}',straight[f'outer{b}_masyr'],expected)
            for score in ('residual_score','acceleration_score','combined_score'):
                close(scenario+f'/straight_invariance_{score}{b}',straight[f'{score}{b}'],null[f'{score}{b}'],1e-8,1e-8)
            coeff=np.column_stack((-fraction[:,1],fraction[:,0]))
            close(scenario+f'/source_slopes{b}',straight[f'pm_source{b}_masyr']-null[f'pm_source{b}_masyr'],coeff[:,:,None]*expected[:,None,:])
        acceleration=rng.uniform(-.015,.015,(n,2))
        curved_xy=xy+.5*acceleration[:,None,:]*t[None,:,None]**2
        curved=new.observe_pairs(*args,outer_relative_xy_au=curved_xy,mass_fraction=fraction)
        for b in (34,66):
            close(scenario+f'/decomposition{b}',curved[f'pm{b}_masyr'],curved[f'outer{b}_masyr']+curved[f'inner_bias{b}_masyr']+curved[f'noise{b}_masyr'],1e-10,1e-10)
            close(scenario+f'/source_difference{b}',curved[f'pm_source{b}_masyr'][:,1]-curved[f'pm_source{b}_masyr'][:,0],curved[f'pm{b}_masyr'])
            close(scenario+f'/noise_unchanged{b}',curved[f'noise{b}_masyr'],null[f'noise{b}_masyr'],1e-10,1e-10)
            close(scenario+f'/covariance_unchanged{b}',curved[f'cov{b}_masyr2'],null[f'cov{b}_masyr2'],0,0)
        for key in ('cov_delta_relative_masyr2','cov_delta_source_masyr2'):
            close(scenario+'/'+key,curved[key],null[key],0,0)
            assert np.linalg.eigvalsh(curved[key]).min()>0
        for label,kwargs in [('shape',dict(outer_relative_xy_au=zero[:,:-1])),
                             ('nonfinite',dict(outer_relative_xy_au=zero+np.nan)),
                             ('mass',dict(outer_relative_xy_au=zero,mass_fraction=fraction*2)),
                             ('duplicate_slope',dict(outer_relative_xy_au=zero))]:
            invalid_args=args if label!='duplicate_slope' else (velocity,)+args[1:]
            try:new.observe_pairs(*invalid_args,**kwargs)
            except ValueError:records.append(dict(name=scenario+'/reject_'+label,passed=True))
            else:raise AssertionError('Invalid inputs accepted: '+label)
    result=dict(passed=True,checks=len(records),records=records,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'joint_observer.py',HERE/'OBSERVER_PROTOCOL.md',Path(__file__)]},
        fixture_systems=n,scenarios=list(new.SCENARIOS),
        scope='Analytic linear signal, combined fits, additive decomposition and legacy covariance; no gravity decision or population validation.')
    (HERE/'OBSERVER_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=True,checks=len(records))))

if __name__=='__main__':main()
