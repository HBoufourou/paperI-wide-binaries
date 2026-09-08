"""End-to-end SOFTWARE fixture: no physical population or gravity result."""
from pathlib import Path
import contextlib,io
import numpy as np
import engine as e

def run():
    base=e.HERE/'software_fixture';inputs=base/'inputs';out=base/'outputs'
    if base.exists():raise ValueError('Refuse to overwrite fixture execution')
    inputs.mkdir(parents=True);rng=np.random.default_rng(279000010);records=[]
    for gravity,role,component in sorted(e.required_records()):
        n=48;tag=f'{gravity}_{role}_{component}';folder=inputs/tag;folder.mkdir()
        angle=rng.uniform(0,2*np.pi,n);s=rng.uniform(2,30,n)
        pop=dict(rproj_kau=s[:,None]*np.column_stack((np.cos(angle),np.sin(angle))),distance_pc=rng.uniform(50,150,n),
            Mphot1=rng.uniform(.5,1,n),Mphot2=rng.uniform(.1,.5,n),weight_beta1=np.ones(n),weight_beta1p6=rng.uniform(.5,1.5,n),weight_beta1p3=rng.uniform(.8,1.2,n))
        pp=folder/'population.npz';np.savez_compressed(pp,**pop);obsfiles={};hashes={}
        for si,scenario in enumerate(e.SCENARIOS):
            amplitude=.1 if gravity=='Newton' else .3
            pm=rng.normal(0,amplitude,(n,2));delta=rng.normal(0,.01,(n,2));source=np.stack((-pm/2,pm/2),axis=1)
            obs=dict(pm34_masyr=pm-delta,pm66_masyr=pm,pm_source34_masyr=source-np.stack((-delta/2,delta/2),axis=1),pm_source66_masyr=source,
                cov66_masyr2=np.tile(np.eye(2)*.0004,(n,1,1)),cov_delta_relative_masyr2=np.tile(np.eye(2)*.0002,(n,1,1)),
                cov_delta_source_masyr2=np.tile(np.eye(2)*.0001,(n,2,1,1)))
            for key in e.SOURCE_SCORES:obs[key]=rng.uniform(.1,5 if component=='pure' else 50,(n,2))
            op=folder/(scenario+'.npz');np.savez_compressed(op,**obs)
            obsfiles[scenario]=op.relative_to(inputs).as_posix();hashes[scenario]=e.sha(op)
        records.append(dict(gravity=gravity,role=role,component=component,prior_au=None if component=='pure' else int(component[4:]),
            population_path=pp.relative_to(inputs).as_posix(),observations=obsfiles,sha256=dict(population=e.sha(pp),observations=hashes),
            metadata=dict(software_fixture=True,fixture_seed=279000010)))
    mp=inputs/'manifest.json';e.write_json(mp,dict(schema_version=1,purpose='software_fixture',production_authorized=False,
        scenarios=list(e.SCENARIOS),records=records,physical_scope='None: artificial arrays exercise file/provenance/calibration/test plumbing only'))
    refused=False
    try:e.load_banks(mp)
    except ValueError:refused=True
    assert refused
    buffer=io.StringIO()
    with contextlib.redirect_stdout(buffer):
        e.freeze_input(mp,out,True,fixture_repeats=8)
        e.run_calibration(mp,out,True)
        e.run_test(mp,out,True)
    (base/'execution.log').write_text(buffer.getvalue(),encoding='utf-8')
    cal=e.read_json(out/'calibration/manifest.json');test=e.read_json(out/'test/manifest.json')
    assert cal['views']==1248 and test['views']==1728 and test['scientific_status']=='software_only'
    decisions=[]
    for rec in test['records']:
        data=e.load_npz(out/'test'/rec['path']);decisions.append(data['decision'])
        assert np.all(data['decision']==-1) #8 calibration draws cannot reach a Bonferroni p<=.01.
    e.write_json(base/'VERIFICATION.json',dict(status='passed_software_only',gravity_evidence=False,production_mode_refused_fixture=True,
        fixture_seed=279000010,fixture_repeats=8,calibration_views=cal['views'],test_views=test['views'],all_decisions_indeterminate=True,
        fixture_calibration_seed_rule='279100001+1000*case_index',fixture_test_seed_rule='279300001+1000*case_index',
        fixture_file_sha256=e.sha(__file__),engine_sha256=e.sha(e.HERE/'engine.py'),protocol_sha256=e.sha(e.HERE/'PROTOCOL.md'),
        input_manifest_sha256=e.sha(mp),test_manifest_sha256=e.sha(out/'test/manifest.json'),
        scope='Full file/provenance/template/profile/rank-calibration/test pipeline on artificial arrays. Eight calibration repeats intentionally cannot pass alpha=.01; no physical model or efficacy is tested.'))
    print('END-TO-END SOFTWARE FIXTURE PASSED; physical evidence = false')

if __name__=='__main__':run()
