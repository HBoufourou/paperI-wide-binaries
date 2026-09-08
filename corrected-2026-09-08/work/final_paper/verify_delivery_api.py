"""Bounded replay of the released measured-only catalogue interface.

This does not regenerate physical banks or duplicate the independent full
statistical reconstruction. Run from a relocated subset or complete archive.
"""
from pathlib import Path
import argparse,json,hashlib,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'stats'))
import engine,features

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def run(out):
    start=time.perf_counter();results=HERE/'stats/results';checks=[]
    def check(name,ok):
        checks.append(dict(name=name,passed=bool(ok)))
        if not ok:raise AssertionError(name)
    declaration=engine.check_declaration();execution=read(results/'EXECUTION_PROTOCOL.json')
    check('declaration',sha(HERE/'stats/STATISTICAL_DECLARATION.json')==execution['statistical_declaration_sha256'])
    check('bank_manifest',sha(HERE/'bank_manifest.json')==execution['bank_manifest_sha256'])
    template_meta=read(results/'template_manifest.json')
    check('templates',sha(results/'templates.npz')==template_meta['template_sha256'])
    cm=read(results/'calibration/manifest.json');tm=read(results/'test/manifest.json')
    check('calibration_reference',sha(results/'calibration/reference.npz')==cm['reference_sha256'])
    check('calibration_manifest',sha(results/'calibration/manifest.json')==tm['calibration_manifest_sha256'])
    records=read(HERE/'bank_manifest.json')['records'];case=tm['records'][3]
    check('test_case',sha(results/'test'/case['path'])==case['sha256'])
    saved=engine.load_npz(results/'test'/case['path']);template=engine.load_npz(results/'templates.npz')
    reference=engine.load_npz(results/'calibration/reference.npz')['sorted']
    pools=[];observed=[]
    for component in case['lineage']['components']:
        row=next(r for r in records if (r['gravity'],r['role'],r['realization'],r['component'])==(case['gravity'],'test',case['realization'],component))
        pp=HERE/row['population_path'];check(component+'/population',sha(pp)==row['sha256']['population'])
        pools.append(engine.load_npz(pp));views={}
        for scenario in engine.SCENARIOS:
            op=HERE/row['observations'][scenario];check(component+'/'+scenario,sha(op)==row['sha256']['observations'][scenario]);views[scenario]=engine.load_npz(op)
        observed.append(views)
    parent={k:np.concatenate([p[k] for p in pools]) for k in features.BASE_FIELDS}
    n=0
    for ri in (0,1,499,999):
        indices=saved['parent_index'][ri];p={k:v[indices] for k,v in parent.items()}
        for si,scenario in enumerate(engine.SCENARIOS):
            for mi,mode in enumerate(engine.MODES):
                o={k:np.concatenate([v[scenario][k] for v in observed])[indices] for k in features.OBS_FIELDS[mode]}
                for rule_index,rule in enumerate(engine.RULES):
                    reply=engine.predict_catalogue(p,o,scenario,template,reference,mode,rule,declared_scope_confirmed=True)
                    label=f'{ri}/{scenario}/{mode}/{rule}'
                    expected=int(saved['decision'][rule_index,mi,si,ri]);reason=int(saved['reason'][rule_index,mi,si,ri])
                    check(label+'/supported',reply['status']=='conditional' and reply['sky_ready'] is False)
                    check(label+'/p_values',np.array_equal(reply['p_values'],saved['p_joint'][rule_index,mi,si,ri]))
                    check(label+'/decision',reply['decision']==('indeterminate' if expected<0 else engine.GRAVITIES[expected]))
                    check(label+'/reason',reply['reason']==('both_retained','both_rejected','one_retained')[reason]);n+=1
    for args in [({},{}),({'distance_pc':np.ones(3)},{}),({k:v[:1000] for k,v in parent.items()},{})]:
        reply=engine.predict_catalogue(*args,'uniform_white',template,reference)
        check('unsupported/'+str(len(checks)),reply['status']=='unsupported' and reply['decision']=='indeterminate' and reply['sky_ready'] is False)
    record=dict(passed=True,checks_count=len(checks),checks=checks,catalogues=4,physical_case_index=3,
        scenario_views=12,method_rule_calls=n,all_p_values_and_decisions_exact=True,source_sha256=sha(__file__),
        runtime_seconds=time.perf_counter()-start,physical_banks_regenerated=False,full_statistics_reconstructed=False,
        scope='Bounded public-interface replay from saved measured fields; exact calibration values and classifications, including reduced-mode field isolation. No external Gaia validation.')
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k!='checks'}))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=HERE/'review/delivery_api_verification.json');run(parser.parse_args().out)
