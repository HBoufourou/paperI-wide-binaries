"""Prospective namespace audit; no physical generation or decision sampling."""
from pathlib import Path
import ast, hashlib, json
import engine

HERE = Path(__file__).resolve().parent
CHAIN = HERE.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def arithmetic(node, names):
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return node.value
    if isinstance(node, ast.Name):
        return names[node.id]
    if isinstance(node, ast.BinOp):
        a, b = arithmetic(node.left, names), arithmetic(node.right, names)
        if isinstance(node.op, ast.Add): return a+b
        if isinstance(node.op, ast.Mult): return a*b
        if isinstance(node.op, ast.Sub): return a-b
    raise ValueError('Not a declared arithmetic seed expression: '+ast.dump(node))


def main():
    checks = []
    def check(name, value):
        checks.append(dict(name=name, passed=bool(value)))

    actual = list(engine.calibration_cases())+list(engine.test_cases())
    predicted = []
    mixtures = [(f, eta) for f in (0., .075, .15, .225, .3)
                for eta in ((0.,) if f == 0 else (0., .5, 1.))]
    for hi, gravity in enumerate(('Newton', 'QUMOND')):
        for bi, beta in enumerate((1., 1.6)):
            for mi, (f, eta) in enumerate(mixtures):
                predicted.append(dict(gravity=gravity, role='calibration', beta=beta,
                    f=f, eta=eta, prior=None, seed=270000001+1000000*hi+10000*bi+100*mi))
        for bi, beta in enumerate((1., 1.6, 1.3)):
            for pi, prior in enumerate((40, 80, 60, 160)):
                for fi, f in enumerate((0., .075, .225, .3)):
                    if bi == 2 and (prior not in (60, 160) or f not in (0., .3)): continue
                    predicted.append(dict(gravity=gravity, role='test', beta=beta,
                        f=f, eta=0., prior=prior,
                        seed=274000001+3000000*hi+10000*bi+1000*pi+100*fi))
    keys = ('gravity','role','beta','f','eta','prior','seed')
    canon = lambda rows: sorted(tuple(str(row[k]) for k in keys) for row in rows)
    check('independent_enumeration_matches_engine', canon(actual)==canon(predicted))
    science = {r['seed'] for r in actual}
    check('124_unique_decision_streams',len(actual)==len(science)==124)
    check('52_calibration_72_test',sum(r['role']=='calibration' for r in actual)==52 and sum(r['role']=='test' for r in actual)==72)

    source = CHAIN/'generate_banks.py'
    tree = ast.parse(source.read_text(encoding='utf-8'))
    expr = {node.targets[0].id: node.value for node in ast.walk(tree)
            if isinstance(node, ast.Assign) and len(node.targets)==1
            and isinstance(node.targets[0],ast.Name) and node.targets[0].id in ('ps','os')}
    upstream = []
    for gi, gravity in enumerate(('Newton','QUMOND')):
        for ri, role in enumerate(('fitting','calibration','test')):
            for ci, component in enumerate(('pure','comp40','comp80','comp60','comp160')):
                if ci>=3 and role!='test': continue
                for kind, base, variable in [('population',260000001,'ps'),('observer',263000001,'os')]:
                    seed = base+gi*1000000+ri*10000+ci*100
                    check('/'.join((kind,gravity,role,component,'source_expression')),
                          arithmetic(expr[variable],dict(gi=gi,ri=ri,ci=ci))==seed)
                    upstream.append(dict(kind=kind,gravity=gravity,role=role,component=component,seed=seed))
    upstream_seeds = {r['seed'] for r in upstream}
    check('44_unique_upstream_streams',len(upstream)==len(upstream_seeds)==44)
    check('no_decision_upstream_collision',not science.intersection(upstream_seeds))

    fixtures = [
        ('verify_orbits.py','analytic_trajectories',268000001),
        *[('verify_population.py','population_'+c,268000011+j) for j,c in enumerate(('pure','comp40','comp160'))],
        ('verify_population.py','empty_helper',268000020),
        *[('validate_physical_chain.py','initial_'+g,268000101+1000*gi) for gi,g in enumerate(('newton','qumond'))],
        *[('validate_physical_chain.py','observer_'+g,268000201+1000*gi) for gi,g in enumerate(('newton','qumond'))],
        ('verify_observer.py','analytic_inputs',269000001),
        *[('verify_observer.py','paired_noise_'+s,269000002+si) for si,s in enumerate(('uniform_white','clustered_correlated','sparse_floor'))],
        ('review/verify_interpolant.py','polynomial_fixture',268900001),
        ('review/export_model_probes.py','interpolation_probes',268900101),
        ('review/check_energy.py','executed_analytic_action_fixture',275000001),
    ]
    fixture_seeds = {r[2] for r in fixtures}
    check('16_unique_upstream_fixture_streams',len(fixtures)==len(fixture_seeds)==16)
    check('no_decision_fixture_collision',not science.intersection(fixture_seeds))
    check('no_upstream_fixture_collision',not upstream_seeds.intersection(fixture_seeds))
    check('275000001_reserved_only_for_executed_analytic_fixture',275000001 in fixture_seeds and 275000001 not in science|upstream_seeds)

    sources = sorted([*CHAIN.glob('*.py'),*(CHAIN/'review').glob('*.py'),*(CHAIN/'physics').glob('*.py')])
    candidates = []
    for path in sources:
        lines = path.read_text(encoding='utf-8').splitlines()
        for node in ast.walk(ast.parse('\n'.join(lines))):
            if isinstance(node,ast.Constant) and type(node.value) is int and 250000000<=node.value<300000000:
                candidates.append(dict(path=path.relative_to(CHAIN).as_posix(),line=node.lineno,
                    value=node.value,expression_line=lines[node.lineno-1].strip()))
    expected_literals = {260000001,263000001,268000000,268000001,268000011,268000020,
        268000101,268000201,269000000,269000001,269000002,268900001,268900101,275000001}
    # A verifier added after scientific freeze deliberately replays the existing
    # streams. Check its complete case declaration before classifying its seed
    # constants as a replay; no random catalogue is sampled by this audit.
    replay_path=CHAIN/'review/verify_decisions.py'
    replay_records=[]
    if replay_path.exists():
        replay_tree=ast.parse(replay_path.read_text(encoding='utf-8'))
        function=next(n for n in replay_tree.body if isinstance(n,ast.FunctionDef) and n.name=='declared_cases')
        namespace=dict(GRAVITIES=('Newton','QUMOND'),BETAS=(1.,1.6),PRIORS=(40,80,60,160),FTEST=(0.,.075,.225,.3),
            MIXES=tuple(mixtures),THETA=tuple((b,f,e) for b in (1.,1.6) for f,e in mixtures))
        exec(compile(ast.Module(body=[function],type_ignores=[]),str(replay_path),'exec'),namespace)
        replay_cal,replay_test=namespace['declared_cases']();replay_records=replay_cal+replay_test
        check('independent_verifier_replays_exact_124_cases',canon(replay_records)==canon(actual))
    unknown = [r for r in candidates if r['value'] not in expected_literals and not
        (r['path']=='review/verify_decisions.py' and r['value'] in (270000001,274000001) and replay_records)]
    check('all_large_integer_candidates_accounted_for',not unknown)
    check('analytic_fixture_literal_present',any(r['path']=='review/check_energy.py' and r['value']==275000001 for r in candidates))

    software_seeds = {279000001,279000002,279000003,279000010}
    software_seeds.update(279100001+1000*i for i in range(52))
    software_seeds.update(279300001+1000*i for i in range(72))
    check('corrected_software_streams_distinct',len(software_seeds)==128 and not software_seeds.intersection(science|upstream_seeds|fixture_seeds))
    archive = HERE/'software_fixture'/'executed_source'
    proof = json.loads((HERE/'software_fixture'/'VERIFICATION.json').read_text())
    # Exact executed source files remain available even though future scheduling changed.
    check('preserved_executed_engine_hash',sha(archive/'engine.py')==proof['engine_sha256'])
    check('preserved_executed_protocol_hash',sha(archive/'PROTOCOL.md')==proof['protocol_sha256'])
    result = dict(passed=all(r['passed'] for r in checks),checks=checks,
        scientific_execution_performed=False,decision_streams=actual,upstream_streams=upstream,
        fixture_streams=[dict(path=p,purpose=n,seed=s) for p,n,s in fixtures],
        corrected_software_streams=sorted(software_seeds),large_integer_candidates=candidates,independent_verifier_replay_cases=replay_records,
        unknown_candidates=unknown,source_sha256={p.relative_to(CHAIN).as_posix():sha(p) for p in sources},
        decision_source_sha256={p.name:sha(p) for p in [HERE/'engine.py',HERE/'PROTOCOL.md',HERE/'AMENDMENT_01_TEST_SEEDS.md',Path(__file__)]},
        intentional_reuse=['Same selected parents across cadence observations and comparisons.',
            'Same population seed for predeclared 20k/40k prefix refinement.',
            'Same null noise when verifying analytic versus legacy observers.',
            'Source-seeded deterministic reruns of an existing declared check.',
            'Independent verification replays all124 existing decision streams to compare indices exactly; it does not claim new data.'],
        historical_exception='The first artificial plumbing fixture inherited future decision seeds; its complete unchanged archive is disclosed in SOFTWARE_SEED_AUDIT.md. It is not an independent physical library or scientific test.',
        limitations='Static source snapshot plus explicit current expression enumeration. Not an assertion about arbitrary uninspected future code or RNG substreams. Different seeds and finite-library role separation do not eliminate empirical covariate reuse.')
    out=HERE/'SEED_AUDIT.json';out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],checks=len(checks),decision_streams=len(science),upstream_streams=len(upstream_seeds),fixture_streams=len(fixture_seeds),unknown_candidates=unknown)))
    if not result['passed']:raise AssertionError([r for r in checks if not r['passed']])


if __name__=='__main__': main()
