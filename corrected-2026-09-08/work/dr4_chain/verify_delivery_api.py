"""Bounded saved-catalogue replay through the public API; not independent inference."""
from pathlib import Path
import argparse, hashlib, json, sys
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'decision'))
sys.path.insert(0,str(HERE/'physics'))
import engine
from interaction_model import InteractionModel


def main(out):
    manifest=engine.read_json(HERE/'bank_manifest.json')
    result=HERE/'decision/results'
    spec=engine.read_json(result/'EXECUTION_PROTOCOL.json')
    assert spec['engine_sha256']==engine.sha(HERE/'decision/engine.py')
    assert spec['protocol_sha256']==engine.sha(HERE/'decision/PROTOCOL.md')
    assert spec['bank_manifest_sha256']==engine.sha(HERE/'bank_manifest.json')
    tm=engine.read_json(result/'template_manifest.json')
    assert tm['template_sha256']==engine.sha(result/'templates.npz')
    cm=engine.read_json(result/'calibration/manifest.json')
    assert cm['reference_sha256']==engine.sha(result/'calibration/reference.npz')
    test=engine.read_json(result/'test/manifest.json')
    assert test['calibration_manifest_sha256']==engine.sha(result/'calibration/manifest.json')
    template=engine.load_npz(result/'templates.npz')
    calibration=engine.load_npz(result/'calibration/reference.npz')['sorted']
    records={(r['gravity'],r['role'],r['component']):r for r in manifest['records']}
    checks=[]
    def check(name,ok,**detail):
        checks.append(dict(name=name,passed=bool(ok),**detail))
        if not ok:raise AssertionError(name)
    model=InteractionModel.load()
    for gravity in ('newton','qumond'):
        x=np.array([[1.,.3,.5],[2.,-.4,-.7]])
        check('validated_model/'+gravity,np.isfinite(model.potential(x,.3,gravity)).all() and np.isfinite(model.acceleration(x,.3,gravity)).all())
    selected=[r for r in test['records'] if r['beta']==1. and r['prior']==160 and r['f'] in (0.,.3)]
    assert len(selected)==4
    for rec in selected:
        path=result/'test'/rec['path'];check('test_hash/'+rec['path'],engine.sha(path)==rec['sha256'])
        saved=engine.load_npz(path);lineage=rec['lineage'];ix=saved['parent_index'][0]
        populations=[];observed={s:[] for s in engine.SCENARIOS}
        for component in lineage['components']:
            record=records[(rec['gravity'],'test',component)]
            pp=HERE/record['population_path'];check('population_hash/'+rec['path']+'/'+component,engine.sha(pp)==record['sha256']['population'])
            p=engine.load_npz(pp)
            populations.append({k:p[k] for k in ('distance_pc','Mphot1','Mphot2','rproj_kau')})
            for scenario in engine.SCENARIOS:
                op=HERE/record['observations'][scenario]
                check('observation_hash/'+rec['path']+'/'+component+'/'+scenario,engine.sha(op)==record['sha256']['observations'][scenario])
                o=engine.load_npz(op)
                observed[scenario].append({k:v for k,v in o.items() if v.ndim>0 and len(v)==len(p['distance_pc'])})
        pcat={k:np.concatenate([p[k] for p in populations])[ix] for k in populations[0]}
        for si,scenario in enumerate(engine.SCENARIOS):
            ocat={k:np.concatenate([o[k] for o in observed[scenario]])[ix] for k in observed[scenario][0]}
            answer=engine.infer_observed(pcat,ocat,scenario,template,calibration,True)
            expected=int(saved['decision'][si,0]);label=engine.GRAVITIES[expected] if expected>=0 else 'indeterminate'
            expected_reason=('both_retained','both_rejected','one_retained')[int(saved['reason'][si,0])]
            check('public_api/'+rec['path']+'/'+scenario,answer['status']=='conditional' and answer['decision']==label and answer['reason']==expected_reason and np.allclose(answer['p_values'],saved['p_joint'][si,0],rtol=0,atol=1e-14),decision=label,p_values=answer['p_values'])
            counts=np.bincount(engine.encode_features(engine.measured_features(pcat,ocat)),minlength=engine.NBIN)
            check('saved_counts/'+rec['path']+'/'+scenario,np.array_equal(counts,saved['counts'][si,0]))
            refusal=engine.infer_observed(pcat,ocat,scenario,template,calibration,False)
            check('scope_refusal/'+rec['path']+'/'+scenario,refusal['status']=='unsupported' and refusal['decision']=='indeterminate')
    value=dict(passed=True,checks=checks,checks_count=len(checks),saved_parent_catalogues=4,paired_views=12,
        scope='Bounded replay through the relocated production public API; independent decision reconstruction is a separate review.',
        full_banks_regenerated=False,full_calibration_rerun=False,sky_ready=False,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    engine.write_json(out,value)
    print(json.dumps({k:v for k,v in value.items() if k!='checks'}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=HERE/'PUBLIC_API_VALIDATION.json')
    main(parser.parse_args().out)
