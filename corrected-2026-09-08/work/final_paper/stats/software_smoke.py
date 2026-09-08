"""Complete file pipeline on isolated artificial banks; no physical evidence."""
from pathlib import Path
import csv,hashlib,json,shutil,subprocess,sys,time
import numpy as np
import engine
from test_engine import fixture

HERE=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    start=time.perf_counter();base=HERE/'software_fixture'
    if base.exists():raise ValueError('Preserve previous fixture evidence; refuse overwrite')
    source=base/'source/stats';source.mkdir(parents=True)
    for name in ('engine.py','features.py'):shutil.copyfile(HERE/name,source/name)
    shutil.copyfile(HERE.parent/'PROTOCOL.md',source.parent/'PROTOCOL.md')
    inputs=base/'inputs';inputs.mkdir();records=[];rng=np.random.default_rng(310000010)
    for i,key in enumerate(sorted(engine.required_records())):
        gravity,role,realization,component=key;p,o=fixture(32)
        p['distance_pc']=rng.uniform(40,170,32);p['rproj_kau']=rng.uniform(1,12,(32,2));p['energy']=-rng.uniform(.3,2,32)
        p['weight_beta1']=rng.uniform(.2,2,32)
        for beta in (1.3,1.6):p[engine.weight_key(beta)]=p['weight_beta1']*(-p['energy'])**(beta-1)
        pp=inputs/f'bank{i:02d}_population.npz';np.savez_compressed(pp,**p)
        obs={};hashes={}
        for si,scenario in enumerate(engine.SCENARIOS):
            a={k:v.copy() for k,v in o.items()}
            for b in (34,66):
                a[f'pm{b}_masyr']+=rng.normal(0,.1,(32,2))
                a[f'pm_source{b}_masyr']+=rng.normal(0,.03,(32,2,2))
                a[f'residual_score_source{b}']+=rng.uniform(0,15,(32,2))
            op=inputs/f'bank{i:02d}_{scenario}.npz';np.savez_compressed(op,**a)
            obs[scenario]=op.relative_to(base).as_posix();hashes[scenario]=sha(op)
        records.append(dict(gravity=gravity,role=role,realization=realization,component=component,
            prior_au=None if component=='pure' else int(component[4:]),population_path=pp.relative_to(base).as_posix(),observations=obs,
            sha256=dict(population=sha(pp),observations=hashes),metadata=dict(software_only=True,input_seed=310000010)))
    manifest=base/'manifest.json';manifest.write_text(json.dumps(dict(purpose='software_fixture',production_authorized=False,
        scenarios=list(engine.SCENARIOS),records=records),indent=2)+'\n',encoding='utf-8')
    out=base/'outputs'
    with (base/'execution.log').open('w',encoding='utf-8') as log:
        for stage in ('declare','freeze','calibrate','test'):
            args=[sys.executable,str(source/'engine.py'),stage]
            if stage!='declare':args+=['--manifest',str(manifest),'--out',str(out)]
            if stage=='freeze':args+=['--software-fixture','--fixture-repeats','2']
            result=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,text=True)
            if result.returncode:raise RuntimeError('Software fixture stage failed; inspect preserved execution.log: '+stage)
    rows=list(csv.DictReader((out/'metrics.csv').open()))
    assert len(rows)==5400 and all(int(r['indeterminate'])==2 and int(r['both_retained'])==2 for r in rows)
    assert len(list(csv.DictReader((out/'paired.csv').open())))==7020
    cal=json.loads((out/'calibration/manifest.json').read_text());test=json.loads((out/'test/manifest.json').read_text())
    assert cal['catalogue_views']==1404 and test['catalogue_views']==1080
    e=json.loads((out/'EXECUTION_PROTOCOL.json').read_text())
    streams=[x['seed'] for x in e['calibration_cases']+e['test_cases']]
    assert len(set(streams))==414 and min(streams)>=310100001 and max(streams)<311000000
    proof=dict(passed=True,software_only=True,physical_evidence=False,fixture_input_seed=310000010,case_streams=streams,
        checks=['Full copied-source pipeline','44 artificial banks','234 calibration and180 test cases',
            '1404 calibration and1080 test views','5400 metric rows','7020 paired rows','Every decision indeterminate at insufficient rank resolution',
            '414 fixture-only case seeds, no scientific sampling'],
        code_sha256={n:sha(source/n) for n in ('engine.py','features.py')},protocol_sha256=sha(source.parent/'PROTOCOL.md'),
        manifest_sha256=sha(manifest),execution_protocol_sha256=sha(out/'EXECUTION_PROTOCOL.json'),test_manifest_sha256=sha(out/'test/manifest.json'),
        elapsed_seconds=time.perf_counter()-start)
    (base/'VERIFICATION.json').write_text(json.dumps(proof,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in proof.items() if k!='case_streams'}))


if __name__=='__main__':main()
