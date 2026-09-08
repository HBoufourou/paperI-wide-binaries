"""Read-only inventory of external dependencies; writes only its inventory JSON."""
from pathlib import Path
import hashlib, importlib.metadata, json, platform

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    external=[
      ('work/dr4_binary/physics/binary_qumond.py',
       ['work/dr4_chain/physics/interaction_energy.py','work/dr4_chain/physics/validate_model.py'],
       'Shared quadrature/source-field helpers and independent force-reference evaluator; NumPy and standard library only.'),
      *[(f'work/dr4_binary/physics/results/jobs/reference_E1_q{q}_d1_t45_reference.json',
       ['work/dr4_chain/physics/validate_model.py']+(['work/dr4_chain/review/check_energy.py'] if q!='0.1' else []),
       'Frozen earlier direct two-source force reference; numeric payload is read without importing an earlier experiment.') for q in ('0.1','0.3','1')],
      ('work/dr4_binary/companions/epoch_operator_extended.py',
       ['work/dr4_chain/build_observer.py','work/dr4_chain/verify_observer.py'],
       'Exact source for the derived observer and frozen comparator; the builder checks its SHA256. Standalone NumPy and standard library.'),
      ('work/originality/dr4/run_comparison.py',
       ['work/dr4_chain/population.py','work/dr4_chain/generate_banks.py'],
       'Only source_data() supplies the inherited, already selected empirical covariates. Import also requires pandas and epoch_operator.py; the old experiment main is not executed.'),
      ('work/originality/dr4/epoch_operator.py',
       ['work/originality/dr4/run_comparison.py'],
       'Required by the old selection-provider module import. Its old orbit/observation experiment is not executed.'),
      ('work/reanalysis/data/observation_arrays.npz',
       ['work/originality/dr4/run_comparison.py','work/dr4_chain/generate_banks.py','work/dr4_chain/review/verify_populations.py'],
       'Frozen reduced empirical observed catalogue used as covariates, not as the new synthetic decision truth.')
    ]
    rows=[]
    for relative,consumers,why in external:
        path=ROOT/relative
        if not path.is_file():raise FileNotFoundError(path)
        rows.append(dict(path=relative,bytes=path.stat().st_size,sha256=sha(path),consumers=consumers,purpose=why))
    core=[]
    for path in sorted(HERE.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts:continue
        if 'software_fixture_initial_seed_audit' in path.parts or 'software_fixture' in path.parts:continue
        if path.suffix=='.py' or (path.suffix=='.md' and any(t in path.name for t in ('PROTOCOL','AMENDMENT','CONTRACT','DERIVATION','API'))):
            core.append(dict(path=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=sha(path)))
    bank=HERE/'bank_manifest.json';bp=json.loads(bank.read_text()) if bank.exists() else {}
    status={}
    for name in ('physics/validation.json','physical_validation.json','POPULATION_IMPLEMENTATION_VALIDATION.json'):
        path=HERE/name
        if path.exists():
            value=json.loads(path.read_text());status[name]=dict(present=True,sha256=sha(path),passed=value.get('passed'),status=value.get('status'))
        else:status[name]=dict(present=False)
    result=dict(schema_version=1,scope='Minimal external-file closure for the entire new dr4_chain subtree, including independent controls; this is an inventory, not a frozen completion certificate.',
        path_base='Package/project root containing work/',new_subtree=dict(path='work/dr4_chain',include='Preserve the entire new subtree, including protocols, sources, saved libraries, manifests, controls and disclosed software history.',exclude=['**/__pycache__/**','**/*.pyc']),
        external_files=rows,external_file_count=len(rows),external_total_bytes=sum(r['bytes'] for r in rows),
        runtime=dict(executed_python=platform.python_version(),numpy=importlib.metadata.version('numpy'),pandas=importlib.metadata.version('pandas'),
            required=['Python standard library','NumPy','pandas (inherited source_data module import)'],
            independent_review=['Python standard library','NumPy'],network_required=False,optional_plotting_dependencies=[]),
        source_snapshot=core,
        execution_status_snapshot=dict(proofs=status,bank_manifest_present=bank.exists(),bank_record_count=len(bp.get('records',[])),
            production_authorized=bp.get('production_authorized',False),decision_frozen=(HERE/'decision/results/EXECUTION_PROTOCOL.json').exists(),
            decision_calibration_complete=(HERE/'decision/results/calibration/manifest.json').exists(),decision_test_complete=(HERE/'decision/results/test/manifest.json').exists(),
            warning='Status is a read-only point-in-time snapshot while upstream work may continue; final manifests and reviews control completion.'),
        unnecessary_old_trees=['Earlier physical/paired/decision libraries','Earlier companion regression models and calibration/test catalogues','Earlier whole ZIP archives','Earlier axisymmetric solver tree','Earlier manuscript and plotting outputs'],
        raw_catalogue_note='The frozen observation_arrays.npz is sufficient for this chain. Regenerating that earlier reduction from the official CSV is a separate upstream task with its own dependencies; it is not silently required here.',
        path_portability=dict(runtime='New scientific bank population/observation paths are relative to work/dr4_chain; code locates dependencies from __file__. Run documented commands from the project root.',
            metadata='Some archived proof/model/manifest path strings are execution metadata. Decision bank_manifest_path in EXECUTION_PROTOCOL.json is descriptive; runtime uses the explicit --manifest path and its bytes hash.',
            resume='Existing GENERATION_PROTOCOL.json and decision execution protocol are intentionally immutable. Do not edit paths or hashes to bypass a mismatch. New generation in a clean source copy produces its own provenance; preserve original evidence separately.',
            hashing='Elapsed times and compressed-container details can differ between independent regeneration runs. Compare scientific arrays/numerics with declared tolerances and retain each run\'s own hashes; never claim identical old-file hashes without verifying bytes.'),
        inventory_source_sha256=sha(Path(__file__)))
    out=HERE/'PACKAGE_DEPENDENCIES.json';out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(external_files=len(rows),external_bytes=result['external_total_bytes'],source_snapshot_files=len(core),status=result['execution_status_snapshot'])))


if __name__=='__main__':main()
