"""Verify saved evidence integrity and summarize the executed decision gates."""
from pathlib import Path
import csv, hashlib, json

HERE=Path(__file__).resolve().parent


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    checks=[];evidence={}
    def check(name,ok):
        checks.append(dict(name=name,passed=bool(ok)))
        if not ok:raise AssertionError(name)
    proofs={
        'physics/validation.json':('passed',),
        'OBSERVER_VALIDATION.json':('passed',),
        'POPULATION_IMPLEMENTATION_VALIDATION.json':('passed',),
        'ORBIT_IMPLEMENTATION_VALIDATION.json':('passed',),
        'physical_validation.json':('passed',),
        'review/ENERGY_REPORT.json':('passed','completed'),
        'review/interpolant_refined_verification.json':('passed','complete'),
        'review/grid_refined_verification.json':('passed',),
        'review/physics_report_verification.json':('report_reconstruction_passed','physics_criteria_passed','complete'),
        'review/population_verification.json':('passed','complete'),
        'review/decision_verification.json':('passed','complete'),
        'decision/SEED_AUDIT.json':('passed',),
    }
    for name,fields in proofs.items():
        value=read(HERE/name)
        for field in fields:check(name+'/'+field,value.get(field) is True)
        evidence[name]=dict(sha256=sha(HERE/name),checks_count=value.get('checks_count',value.get('checks') if isinstance(value.get('checks'),int) else len(value.get('checks',[]))))
    physics=read(HERE/'physics/validation.json')
    check('scalar_model',physics['model_sha256']==sha(HERE/'physics/results/interaction_model.npz'))
    check('scalar_evaluator',physics['evaluator_sha256']==sha(HERE/'physics/interaction_model.py'))
    physical=read(HERE/'physical_validation.json')
    check('physical_validated_model',physical['scalar_model_sha256']==physics['model_sha256'])
    check('physical_physics_proof',physical['physics_validation_sha256']==sha(HERE/'physics/validation.json'))
    for name,digest in physical['implementation_proof_sha256'].items():check('implementation_proof/'+name,sha(HERE/name)==digest)
    check('physical_scope_not_sky_ready',physical['sky_ready'] is False and physical['astrophysical_population_validated'] is False)
    for name,digest in physical['source_sha256'].items():
        path=HERE/'physics'/name if name=='interaction_model.py' else HERE/name
        check('physical_source/'+name,sha(path)==digest)
    for record in physical['files']:check('trajectory_file/'+record['path'],sha(HERE/record['path'])==record['sha256'])
    generation=read(HERE/'GENERATION_PROTOCOL.json')
    for name,digest in generation['source_sha256'].items():check('generation_source/'+name,sha(HERE/name)==digest)
    check('generation_physical_proof',generation['physical_validation_sha256']==sha(HERE/'physical_validation.json'))
    bank=read(HERE/'bank_manifest.json')
    check('22_authorized_banks',bank['production_authorized'] is True and len(bank['records'])==22)
    check('bank_physical_proof',bank['physical_validation']['sha256']==sha(HERE/bank['physical_validation']['path']))
    for record in bank['records']:
        name=record['population_path'];check('bank/'+name,sha(HERE/name)==record['sha256']['population'])
        for scenario,name in record['observations'].items():check('bank/'+name,sha(HERE/name)==record['sha256']['observations'][scenario])
    result=HERE/'decision/results';execution=read(result/'EXECUTION_PROTOCOL.json')
    check('production_not_fixture',execution['fixture'] is False and execution['purpose']=='physics_production')
    check('execution_engine',execution['engine_sha256']==sha(HERE/'decision/engine.py'))
    check('execution_protocol',execution['protocol_sha256']==sha(HERE/'decision/PROTOCOL.md'))
    check('execution_banks',execution['bank_manifest_sha256']==sha(HERE/'bank_manifest.json'))
    template=read(result/'template_manifest.json')
    check('template',template['template_sha256']==sha(result/'templates.npz'))
    check('template_execution',template['execution_protocol_sha256']==sha(result/'EXECUTION_PROTOCOL.json'))
    calibration=read(result/'calibration/manifest.json');test=read(result/'test/manifest.json')
    check('calibration_dimensions',len(calibration['records'])==52 and calibration['views']==156000)
    check('test_dimensions',len(test['records'])==72 and test['views']==216000)
    for stage,manifest in [('calibration',calibration),('test',test)]:
        check(stage+'/execution',manifest['execution_protocol_sha256']==sha(result/'EXECUTION_PROTOCOL.json'))
        for record in manifest['records']:check(stage+'/'+record['path'],sha(result/stage/record['path'])==record['sha256'])
    check('calibration_reference',calibration['reference_sha256']==sha(result/'calibration/reference.npz'))
    check('test_calibration',test['calibration_manifest_sha256']==sha(result/'calibration/manifest.json'))
    check('metrics',test['metrics_sha256']==sha(result/'metrics.csv'))
    check('paired',test['paired_sha256']==sha(result/'paired_cadences.csv'))
    review=read(HERE/'review/decision_verification.json')
    check('decision_review_source',review['verifier_sha256']==sha(HERE/'review/verify_decisions.py'))
    check('reviewed_engine',review['production_engine_sha256']==sha(HERE/'decision/engine.py'))
    check('reviewed_protocol',review['production_protocol_sha256']==sha(HERE/'decision/PROTOCOL.md'))
    for name,digest in review['input_sha256'].items():check('decision_review_input/'+name,sha(HERE/name)==digest)
    for name,source in [('review/population_verification.json','review/verify_populations.py'),
                        ('review/interpolant_refined_verification.json','review/verify_interpolant.py'),
                        ('review/grid_refined_verification.json','review/verify_energy_grid.py'),
                        ('review/physics_report_verification.json','review/verify_physics_report.py')]:
        proof=read(HERE/name);check('review_source/'+name,proof['source_sha256']==sha(HERE/source))
        if 'model_sha256' in proof:check('reviewed_model/'+name,proof['model_sha256']==physics['model_sha256'])
    with (result/'metrics.csv').open(encoding='utf-8',newline='') as f:metrics=list(csv.DictReader(f))
    check('216_test_cells',len(metrics)==216)
    groups={}
    predicates={
        'all':lambda r:True,
        'fitted_family':lambda r:r['family']=='endpoint' and r['beta_family']=='endpoint',
        'outside_inner_prior':lambda r:r['family']=='outside_inner_prior' and r['beta_family']=='endpoint',
        'outside_beta_grid':lambda r:r['beta_family']=='outside_beta_grid',
    }
    for name,predicate in predicates.items():
        rows=[r for r in metrics if predicate(r)]
        groups[name]=dict(cells=len(rows),count_pass=sum(r['count_target_passed']=='True' for r in rows),
            confidence_pass=sum(r['confidence_target_passed']=='True' for r in rows),
            minimum_correct_rate=min(float(r['correct_rate']) for r in rows),maximum_wrong_rate=max(float(r['wrong_rate']) for r in rows),
            maximum_both_rejected_rate=max(int(r['both_rejected'])/int(r['N_catalogues']) for r in rows),
            minimum_true_retained_rate=min(float(r['true_retained_rate']) for r in rows))
    value=dict(evidence_integrity_passed=True,checks=checks,checks_count=len(checks),evidence=evidence,
        complete_conditional_orbital_population=True,complete_conditional_forward_decision_chain=True,
        populations=440000,observation_views=1320000,calibration_catalogue_views=156000,test_catalogue_views=216000,
        catalogue_size=execution['N'],decision_gate=groups,
        all_declared_cells_meet_count_target=groups['all']['count_pass']==groups['all']['cells'],
        all_declared_cells_meet_confidence_target=groups['all']['confidence_pass']==groups['all']['cells'],
        nominal_one_percent_test_level_established_between_libraries=False,
        stationary_parent_only=True,full_gaia_population_validated=False,continuum_qumond_uniform_error_certified=False,
        actual_DR4_observations=False,sky_ready=False,
        scope='Integrity of executed conditional evidence; this assembler does not independently repeat physical or statistical calculations.',
        source_sha256=sha(Path(__file__)))
    (HERE/'STATUS.json').write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in value.items() if k not in ('checks','evidence')}),flush=True)


if __name__=='__main__':main()
