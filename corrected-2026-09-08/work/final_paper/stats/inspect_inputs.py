"""Read existing inputs for a prospective ablation design; no catalogue sampling."""
from pathlib import Path
import ast, hashlib, json
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CHAIN=ROOT/'work/dr4_chain'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest_path=CHAIN/'bank_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    rows=[];total_bytes=0
    fields={'pm34_masyr':(2,), 'cov34_masyr2':(2,2),
            'residual_score_source34':(2,), 'acceleration_score_source34':(2,)}
    for record in manifest['records']:
        pop=CHAIN/record['population_path'];total_bytes+=pop.stat().st_size
        with np.load(pop,allow_pickle=False) as z:n=len(z['distance_pc'])
        for scenario,relative in record['observations'].items():
            path=CHAIN/relative;total_bytes+=path.stat().st_size
            with np.load(path,allow_pickle=False) as z:
                available=sorted(z.files)
                checks={k:bool(k in z and z[k].shape==(n,)+shape and np.isfinite(z[k]).all()) for k,shape in fields.items()}
                cov=z['cov34_masyr2'];eigen=float(np.linalg.eigvalsh(cov).min())
                symmetry=float(np.max(abs(cov-cov.swapaxes(-1,-2))))
                checks['cov34_positive_definite']=eigen>0
                checks['cov34_symmetric']=np.allclose(cov,cov.swapaxes(-1,-2),rtol=1e-10,atol=1e-15)
                checks['scores34_nonnegative']=bool(np.all(z['residual_score_source34']>=0) and np.all(z['acceleration_score_source34']>=0))
            rows.append(dict(gravity=record['gravity'],role=record['role'],component=record['component'],scenario=scenario,
                N=n,path=relative,expected_sha256=record['sha256']['observations'][scenario],checks=checks,
                covariance34_min_eigenvalue=eigen,covariance34_max_asymmetry=symmetry,available_fields=available))
    timings=[r['metadata']['elapsed_seconds'] for r in manifest['records']]
    old_results=CHAIN/'decision/results'
    storage={name:sum(p.stat().st_size for p in (old_results/name).rglob('*') if p.is_file()) for name in ('calibration','test')}
    candidates=[];syntax_errors=[];scanned=0
    for path in (ROOT/'work').rglob('*.py'):
        if '__pycache__' in path.parts or 'vendor' in path.parts or HERE in path.parents:continue
        try:tree=ast.parse(path.read_text(encoding='utf-8-sig'))
        except (SyntaxError,UnicodeError) as error:
            syntax_errors.append(dict(path=path.relative_to(ROOT).as_posix(),error=type(error).__name__));continue
        scanned+=1
        for node in ast.walk(tree):
            if isinstance(node,ast.Constant) and type(node.value) is int and 310000000<=node.value<400000000:
                candidates.append(dict(path=path.relative_to(ROOT).as_posix(),line=node.lineno,value=node.value))
    result=dict(purpose='Read-only preparation for a prospective new protocol. No new samples, scores or decisions.',
        bank_manifest_sha256=sha(manifest_path),authorized_old_manifest=manifest['production_authorized'],banks=len(manifest['records']),views=len(rows),
        input_checks_passed=all(all(r['checks'].values()) for r in rows),checks=sum(len(r['checks']) for r in rows),rows=rows,
        existing_bank_bytes=total_bytes,existing_decision_storage_bytes=storage,
        recorded_generation_seconds=dict(sum=float(sum(timings)),minimum=float(min(timings)),maximum=float(max(timings))),
        prior_observed_decision_seconds=dict(freeze=11.847059000050649,calibration=89.93158200010657,test=146.20199740002863,
            origin='Recorded original command outputs; these are observations of the old run, not new timing experiments.'),
        namespace_scan=dict(python_files=scanned,integer_range='310000000<=value<400000000',candidates=candidates,syntax_errors=syntax_errors,
            limitation='Read-only literal scan; not a proof about arbitrary future code or values assembled from small integers. New expressions need an exact pre-freeze enumeration.'),
        source_sha256={name:sha(CHAIN/name) for name in ('decision/engine.py','decision/PROTOCOL.md','joint_observer.py','generate_banks.py')},
        inspector_sha256=sha(Path(__file__)))
    (HERE/'INPUT_INSPECTION.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','source_sha256')}))
    if not result['input_checks_passed']:raise AssertionError('Missing or invalid proposed34-month inputs')


if __name__=='__main__':main()
