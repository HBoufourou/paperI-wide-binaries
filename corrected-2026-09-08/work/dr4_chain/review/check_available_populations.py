"""Early row checks on completed bank records; never a final readiness report."""
import argparse,json
from pathlib import Path
from verify_populations import IndependentModel,Checks,load,sha,check_population


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).parents[1]);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();root=args.root.resolve();path=root/'bank_manifest.json'
    text=path.read_text();manifest=json.loads(text);checks=Checks();rows=[]
    model=IndependentModel(root/'physics/results/interaction_model.npz')
    model.empirical=load(root.parent/'reanalysis/data/observation_arrays.npz')
    for record in manifest['records']:
        pp=root/record['population_path'];label='/'.join(record[k] for k in ('gravity','role','component'))
        checks.yes(label+'/population_sha256',sha(pp)==record['sha256']['population'])
        rows.append(check_population(load(pp),model,record['gravity'].lower(),label,checks))
    out={'complete':False,'status':'intermediate_row_verification_only',
         'source_sha256':sha(__file__),'row_checker_sha256':sha(Path(__file__).with_name('verify_populations.py')),
         'model_sha256':sha(model.path),'records_checked':len(rows),'rows_checked':sum(r['rows'] for r in rows),
         'row_checks_passed':all(r['passed'] for r in checks.rows),'checks':checks.rows,
         'scope':'Saved-row snapshot only. Final22-bank/role/seed/observation checks are not asserted here.'}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='checks'}),flush=True)

if __name__=='__main__':main()
