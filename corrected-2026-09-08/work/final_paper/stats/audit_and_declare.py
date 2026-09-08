"""Freeze source/case declaration with an immutable namespace-audit history."""
from pathlib import Path
import hashlib,json
import engine

HERE=Path(__file__).resolve().parent


def main():
    declaration=engine.declare();records=declaration['calibration_cases']+declaration['test_cases']
    science=[r['seed'] for r in records];physical=declaration['physical_streams']
    fixture=declaration['fixture_unit_seeds']+[declaration['fixture_input_seed']]+[310100001+1000*i for i in range(414)]
    old=engine.read(HERE.parents[1]/'dr4_chain/decision/SEED_AUDIT.json')
    historical=[r['seed'] for k in ('decision_streams','upstream_streams','fixture_streams') for r in old[k]]+old['corrected_software_streams']
    groups=dict(scientific_catalogues=science,fresh_physical=physical,software_fixtures=fixture,historical_chain=historical)
    checks=[]
    def check(name,value):checks.append(dict(name=name,passed=bool(value)))
    for name,values in groups.items():check(name+'/within_unique',len(values)==len(set(values)))
    names=list(groups)
    for i,a in enumerate(names):
        for b in names[i+1:]:check(a+'/'+b+'/disjoint',not set(groups[a])&set(groups[b]))
    check('case_counts',len(declaration['calibration_cases'])==234 and len(declaration['test_cases'])==180)
    check('main_science_namespace',min(science)>=320000000 and max(science)<327000000)
    manifest_path=HERE.parent/'bank_manifest.json';manifest=engine.read(manifest_path)
    executed=[]
    for r in manifest['records']:
        if r['origin']=='fresh_physical_generation':
            spec=r['metadata']['source_specification'];executed.extend([spec['population_seed'],spec['observer_seed']])
    check('64_executed_physical_streams_match',len(executed)==64 and set(executed)==set(physical))
    proof=engine.read(HERE/'software_fixture/VERIFICATION.json')
    check('complete_software_pass',proof['passed'] and proof['software_only'] and not proof['physical_evidence'])
    check('software_same_source',proof['code_sha256']==engine.source_hashes())
    check('software_same_root_protocol',proof['protocol_sha256']==engine.sha(HERE.parent/'PROTOCOL.md'))
    inspection=engine.read(HERE/'INPUT_INSPECTION.json')
    check('preimplementation_historical_literal_scan_clear',not inspection['namespace_scan']['candidates'] and not inspection['namespace_scan']['syntax_errors'])
    result=dict(passed=all(r['passed'] for r in checks),checks=checks,streams=groups,
        statistical_declaration_sha256=engine.sha(HERE/'STATISTICAL_DECLARATION.json'),bank_manifest_sha256=engine.sha(manifest_path),
        old_audit_sha256=engine.sha(HERE.parents[1]/'dr4_chain/decision/SEED_AUDIT.json'),
        preimplementation_inspection_sha256=engine.sha(HERE/'INPUT_INSPECTION.json'),software_proof_sha256=engine.sha(HERE/'software_fixture/VERIFICATION.json'),
        source_sha256=engine.source_hashes(),audit_source_sha256=engine.sha(Path(__file__)),
        limitation='Exact declared streams and inherited source snapshot, not a universal proof for arbitrary future arithmetic. Independent verification intentionally replays the same indices; it does not claim fresh experimental data.')
    payload=json.dumps(result,indent=2)+'\n';digest=hashlib.sha256(payload.encode()).hexdigest()
    archive=HERE/'audit_history'/('SEED_AUDIT_'+digest[:16]+'.json');archive.parent.mkdir(exist_ok=True)
    if archive.exists():assert archive.read_text(encoding='utf-8')==payload
    else:archive.write_text(payload,encoding='utf-8')
    (HERE/'SEED_AUDIT.json').write_text(payload,encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],checks=len(checks),group_counts={k:len(v) for k,v in groups.items()},audit_sha256=digest)))
    if not result['passed']:raise AssertionError([r for r in checks if not r['passed']])


if __name__=='__main__':main()
