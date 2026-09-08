"""Report saved conditional decisions without rerunning generation or scoring."""
from pathlib import Path
import argparse,csv,hashlib,json
import numpy as np

def load(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def table(rows,columns):
    result=['| '+' | '.join(columns)+' |','| '+' | '.join('---' for _ in columns)+' |']
    for row in rows:result.append('| '+' | '.join(str(row[k]) for k in columns)+' |')
    return '\n'.join(result)

def run(out):
    out=Path(out);spec=load(out/'EXECUTION_PROTOCOL.json');test=load(out/'test/manifest.json');cal=load(out/'calibration/manifest.json')
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    if digest(out/'metrics.csv')!=test['metrics_sha256']:raise ValueError('Saved metrics hash mismatch')
    if digest(out/'EXECUTION_PROTOCOL.json')!=test['execution_protocol_sha256']:raise ValueError('Execution hash mismatch')
    if digest(out/'calibration/manifest.json')!=test['calibration_manifest_sha256']:raise ValueError('Calibration manifest hash mismatch')
    with (out/'metrics.csv').open(encoding='utf-8',newline='') as f:metrics=list(csv.DictReader(f))
    fixture=spec['fixture'];summary=[];extremes=[]
    for gravity in ('Newton','QUMOND'):
        for scenario in ('uniform_white','clustered_correlated','sparse_floor'):
            for family in ('endpoint','outside_inner_prior'):
                for beta_family in ('endpoint','outside_beta_grid'):
                    selected=[r for r in metrics if r['gravity']==gravity and r['scenario']==scenario and r['family']==family and r['beta_family']==beta_family]
                    if not selected:continue
                    a=min(selected,key=lambda r:float(r['correct_rate']));b=max(selected,key=lambda r:float(r['wrong_rate']))
                    c=min(selected,key=lambda r:float(r['true_retained_rate']))
                    for why,row in [('minimum_correct',a),('maximum_wrong',b),('minimum_true_retained',c)]:
                        extremes.append(dict(extreme=why,**row))
                    summary.append(dict(gravity=gravity,scenario=scenario,family=family,beta_family=beta_family,cells=len(selected),
                        count_pass=sum(r['count_target_passed']=='True' for r in selected),confidence_pass=sum(r['confidence_target_passed']=='True' for r in selected),
                        worst_correct=float(a['correct_rate']),worst_correct_beta=a['beta'],worst_correct_prior=a['prior'],worst_correct_f=a['f'],
                        worst_wrong=float(b['wrong_rate']),worst_wrong_beta=b['beta'],worst_wrong_prior=b['prior'],worst_wrong_f=b['f'],
                        worst_correct_wilson_lo=float(a['correct_wilson_lo']),worst_correct_wilson_hi=float(a['correct_wilson_hi']),
                        worst_wrong_wilson_lo=float(b['wrong_wilson_lo']),worst_wrong_wilson_hi=float(b['wrong_wilson_hi']),
                        worst_true_retained=float(c['true_retained_rate']),worst_true_retained_beta=c['beta'],worst_true_retained_prior=c['prior'],worst_true_retained_f=c['f'],
                        maximum_both_rejected_rate=max(float(r['both_rejected'])/int(r['N_catalogues']) for r in selected),
                        maximum_both_retained_rate=max(float(r['both_retained'])/int(r['N_catalogues']) for r in selected)))
    with (out/'group_summary.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    with (out/'worst_cells.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(extremes[0]));w.writeheader();w.writerows(extremes)
    display=[]
    for r in summary:display.append(dict(Hypothesis=r['gravity'],Cadence=r['scenario'],Population=r['family']+'/'+r['beta_family'],Cells=r['cells'],
        Count_pass=r['count_pass'],Confidence_pass=r['confidence_pass'],Minimum_correct=f"{100*r['worst_correct']:.1f}%",Maximum_wrong=f"{100*r['worst_wrong']:.1f}%",Minimum_true_retained=f"{100*r['worst_true_retained']:.1f}%"))
    title='SOFTWARE FIXTURE ONLY — no physical decision evidence' if fixture else 'Conditional coupled-library decisions'
    primary=[r for r in metrics if r['family']=='endpoint' and r['beta_family']=='endpoint']
    endpoint_beta=[r for r in metrics if r['beta_family']=='endpoint']
    failed_confidence=[r for r in metrics if r['confidence_target_passed']!='True']
    retention=min(primary,key=lambda r:float(r['true_retained_rate']))
    failure_display=[dict(Gravity=r['gravity'],Beta=r['beta'],Prior_AU=r['prior'],f=r['f'],Cadence=r['scenario'],
        Correct=r['correct']+'/'+r['N_catalogues'],Wrong=r['wrong']+'/'+r['N_catalogues'],
        Correct_Wilson_lower=f"{100*float(r['correct_wilson_lo']):.3f}%",Wrong_Wilson_upper=f"{100*float(r['wrong_wilson_hi']):.3f}%") for r in failed_confidence]
    text=[f'# {title}','',
        f"The recorded execution contains {cal['views']:,} calibration views and {test['views']:,} test views, with {spec['repeats']:,} catalogue draws per cell and {spec['N']:,} parents per catalogue. The same parent index arrays couple the three cadence views and both hypothesis scores.",'',
        f"The observed count targets pass in {sum(r['count_target_passed']=='True' for r in metrics)}/{len(metrics)} cells. The stronger marginal Wilson-bound targets pass in {len(metrics)-len(failed_confidence)}/{len(metrics)}. For the fitted beta endpoints, including the separately labelled60/160 AU prior tests, the corresponding confidence count is {sum(r['confidence_target_passed']=='True' for r in endpoint_beta)}/{len(endpoint_beta)}. These counts summarize this conditional library experiment, not an unconditional future-survey guarantee.",'',
        'Each row below retains the worst correct, wrong and true-hypothesis-retention rates over its own cells. These extrema may occur at different nuisance settings; group_summary.csv identifies each setting, and worst_cells.csv retains the full corresponding cell and Wilson intervals. metrics.csv contains every individual cell. No average substitutes for a failing cell.','',
        table(display,list(display[0])),'',
        'Count-pass means an observed correct rate of at least95% and wrong rate of at most1%. Confidence-pass additionally requires the Wilson95% lower endpoint for correct decisions to reach95% and the Wilson95% upper endpoint for wrong decisions to be at most1%. These are different statements. The intervals are marginal Monte Carlo summaries, not simultaneous guarantees over all cells.','',
        'The individual cells failing the stronger confidence target are listed here; all individual results remain in metrics.csv.','',
        table(failure_display,list(failure_display[0])) if failure_display else 'No cell fails the stronger confidence target.','',
        'An indeterminate result has a recorded reason: both gravity hypotheses remain compatible, both fail the combined adequacy/preference test, or the measured/provenance contract is unsupported. The engine never forces a decision by choosing the larger likelihood. Near-best nuisance counts and the full profiled deviance surfaces remain saved so loss of separation through nuisance flexibility is visible. A deviance difference of2 used for those descriptive near-best counts is not a confidence threshold.','',
        'The calibrated discrete family uses beta1/1.6 and mixtures of pure,comp40,comp80; real comp60/comp160 and beta1.3 results are separate tests of transfer. The entries with f=0 under different inner-prior names are additional random repetitions of the same pure population, not different physical prior cases. Each fitted-component law is normalized after the upstream selection, so f is the fraction among those selected parents.','',
        'In the saved test records, eta=0 is an unused interface placeholder whenever a physical prior is supplied. It is not the true mixture parameter of every test population: prior40 corresponds to the fitted eta=0 endpoint, prior80 to eta=1, and the actual60/160 priors have no true eta in that fitted family. Interpret the prior column as the simulated physical prior; at f=0 the companion prior and eta are immaterial.','',
        'Rank p-values have their usual conditional calibration property when the new catalogue and references share a calibrating empirical law. Separate finite test and calibration pools need not have that identical law, even when their underlying generators match. The reported test rates therefore measure library transfer rather than inherit an exact finite-sample guarantee. The experiment does not quantify every uncertainty associated with regenerating the entire fitting/calibration pipeline.','',
        f"Across the saved test views, {sum(int(r['both_retained']) for r in metrics):,} indeterminate outputs retain both hypotheses and {sum(int(r['both_rejected']) for r in metrics):,} reject both. Even within the fitted family, the minimum observed true-model retention is {100*float(retention['true_retained_rate']):.3f}% (marginal Wilson95% interval [{100*float(retention['true_retained_wilson_lo']):.3f}%, {100*float(retention['true_retained_wilson_hi']):.3f}%]), at {retention['gravity']}, beta={retention['beta']}, prior={retention['prior']} AU, f={retention['f']}, {retention['scenario']}. Classification success therefore must not be described as demonstrated99% true-model retention in every test law. This worst inspected interval is descriptive, not a simultaneous calibration test.",'',
        'The new external trajectory, inner photocentre motion and shared astrometric noise are coupled before fitting by the upstream observer. This engine only consumes its declared measured fields. Photometric masses, distances, source/noise distributions, the specific stationary parent and its instantaneous selection remain model assumptions. A numerically successful declared-domain result is not a calibration on actual Gaia measurements.','',
        'The prospective photometric-anchor amendment fixes each empirical total photometric mass before deriving true mass and orbital scales. Measured separation includes the initial noiseless photocentre offset. The covariate library nevertheless inherits earlier chance-alignment, RUWE, precision, RV and transverse-velocity selection; this conditional source distribution must not be equated with an unselected future Gaia parent.','',
        '2048 joint bins represent separation, both signed projected velocity components, an astrometric diagnostic and formal precision. Overflow categories and total-concentration1 smoothing keep the score computable but do not prove physical support. EXECUTION_PROTOCOL.json and template_manifest.json retain all component weight/empty-bin/smoothing diagnostics. High-weight or sparse cases must not be described as an effective-size certificate for a larger future catalogue.','',
        'No threshold, bin edge, mixture grid or smoothing constant is chosen after these test outcomes. The source and input hashes are frozen in the execution protocol, with a separate physical-scope authorization upstream. The final public-sky readiness status remains false.','']
    if fixture:text.extend(['Every array in this software run is artificial. Its eight calibration draws cannot yield a Bonferroni-combined p-value at or below0.01, so indeterminate outputs test that the pipeline respects its calibration resolution. Do not interpret any rate or table in this fixture as gravity power.',''])
    (out/'REPORT.md').write_text('\n'.join(text),encoding='utf-8')
    print('REPORTED',len(metrics),'cells; software_only=',fixture)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);run(p.parse_args().out)
