"""Describe saved results, all failures and paired changes; never rescore data."""
from pathlib import Path
import argparse,csv,hashlib,itertools,json,re
import numpy as np

HERE=Path(__file__).resolve().parent
MODES=('full66','no_diagnostics66','full34','covariates_only66','conditional66')
RULES=('baseline','envelope');SCENARIOS=('uniform_white','clustered_correlated','sparse_floor')


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def csvread(p):
    with Path(p).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
def csvwrite(p,rows,fieldnames=None):
    names=list(rows[0]) if rows else fieldnames
    if names is None:raise ValueError('Empty table requires explicit columns')
    with Path(p).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=names);w.writeheader();w.writerows(rows)
def rate(row,key):return float(row[key+'_rate'])
def yes(row,key):return row[key]=='True'
def table(rows):
    if not rows:return 'No entries.'
    keys=list(rows[0]);lines=['| '+' | '.join(keys)+' |','| '+' | '.join('---' for _ in keys)+' |']
    return '\n'.join(lines+['| '+' | '.join(str(row[k]) for k in keys)+' |' for row in rows])
def identity(row):return {k:row[k] for k in ('mode','rule','realization','gravity','beta','f','prior','scenario','seed')}
def family(row):
    if row['beta_family']=='outside_beta_grid':return 'outside_beta_grid'
    return 'outside_inner_prior' if row['companion_family']=='outside_inner_prior' else 'fitted_family'


def summarize(rows,keys):
    groups={}
    for r in rows:groups.setdefault(tuple(r[k] for k in keys),[]).append(r)
    output=[];extremes=[]
    for key,items in groups.items():
        record=dict(zip(keys,key));record['cells']=len(items)
        for label in ('classification_count_pass','retention_count_pass','overall_count_pass',
                      'classification_confidence_pass','retention_confidence_pass','overall_confidence_pass'):
            record[label]=sum(yes(r,label) for r in items)
        for label,chooser in [('correct',min),('wrong',max),('true_retained',min)]:
            worst=chooser(items,key=lambda r:rate(r,label));record['worst_'+label]=rate(worst,label)
            record['worst_'+label+'_identity']=json.dumps(identity(worst),sort_keys=True)
            record['worst_'+label+'_wilson_lo']=float(worst[label+'_wilson_lo'])
            record['worst_'+label+'_wilson_hi']=float(worst[label+'_wilson_hi'])
            extremes.append(dict(extreme=label,**worst))
        record['maximum_both_retained_rate']=max(int(r['both_retained'])/int(r['N_catalogues']) for r in items)
        record['maximum_both_rejected_rate']=max(int(r['both_rejected'])/int(r['N_catalogues']) for r in items)
        record['maximum_raw_fitting_empty_fraction_mean']=max(float(r['raw_fitting_empty_fraction_mean']) for r in items)
        record['maximum_near_best_true_mean']=max(float(r['near_best_true_mean']) for r in items)
        output.append(record)
    return output,extremes


def envelope_effects(out,test):
    rows=[]
    for case in test['records']:
        path=out/'test'/case['path']
        if sha(path)!=case['sha256']:raise ValueError('Test file hash mismatch for descriptive paired analysis')
        with np.load(path,allow_pickle=False) as z:d=z['decision'];rej=z['rejected'];reason=z['reason']
        h=case['hypothesis_index']
        for mi,mode in enumerate(MODES):
            for si,scenario in enumerate(SCENARIOS):
                a=d[1,mi,si];b=d[0,mi,si];n=len(a)
                row={k:case[k] for k in ('realization','gravity','beta','f','prior','seed','beta_family','companion_family')}
                row.update(mode=mode,scenario=scenario,N_catalogues=n,difference_direction='envelope_minus_baseline')
                for label,aa,bb in [('correct',a==h,b==h),('wrong',(a>=0)&(a!=h),(b>=0)&(b!=h)),
                                   ('indeterminate',a<0,b<0),('true_retained',~rej[1,mi,si,:,h],~rej[0,mi,si,:,h])]:
                    row[label+'_gained']=int(np.sum(aa&~bb));row[label+'_lost']=int(np.sum(bb&~aa))
                    row[label+'_difference']=float(aa.mean()-bb.mean())
                row['baseline_double_rejected_to_correct']=int(np.sum((reason[0,mi,si]==1)&(a==h)))
                row['baseline_double_rejected_to_wrong']=int(np.sum((reason[0,mi,si]==1)&(a>=0)&(a!=h)))
                row['baseline_correct_to_both_retained']=int(np.sum((b==h)&(reason[1,mi,si]==0)))
                if row['true_retained_lost']!=0:raise ValueError('Envelope violates nested true-model retention')
                rows.append(row)
    return rows


def information_summary(pairs):
    groups={}
    for row in pairs:
        if row['comparison_type']!='information':continue
        key=tuple(row[k] for k in ('first_rule','first_mode','second_rule','second_mode'))
        groups.setdefault(key,[]).append(row)
    output=[]
    for key,rows in groups.items():
        item=dict(zip(('first_rule','first_mode','second_rule','second_mode'),key));item['cells']=len(rows)
        for label in ('correct','wrong','indeterminate'):
            delta=[float(r[label+'_difference']) for r in rows]
            item[label+'_increased_cells']=sum(v>1e-14 for v in delta)
            item[label+'_unchanged_cells']=sum(abs(v)<=1e-14 for v in delta)
            item[label+'_decreased_cells']=sum(v<-1e-14 for v in delta)
            item[label+'_minimum_difference']=min(delta);item[label+'_maximum_difference']=max(delta)
            for suffix in ('first_only','second_only'):
                item[label+'_'+suffix+'_total']=sum(int(r[label+'_'+suffix]) for r in rows)
        output.append(item)
    return output


def run(out):
    out=Path(out);execution=read(out/'EXECUTION_PROTOCOL.json');test=read(out/'test/manifest.json');cal=read(out/'calibration/manifest.json')
    for name,key in [('metrics.csv','metrics_sha256'),('paired.csv','paired_sha256')]:
        if sha(out/name)!=test[key]:raise ValueError('Saved summary hash mismatch')
    if test['execution_protocol_sha256']!=sha(out/'EXECUTION_PROTOCOL.json') or test['calibration_manifest_sha256']!=sha(out/'calibration/manifest.json'):raise ValueError('Stage provenance mismatch')
    metrics=csvread(out/'metrics.csv');pairs=csvread(out/'paired.csv')
    if len(metrics)!=5400 or len(pairs)!=7020:raise ValueError('Expected complete5400 metrics and7020 pair rows')
    for row in metrics:row['evaluation_family']=family(row)
    groups,extremes=summarize(metrics,('mode','rule'))
    strata,_=summarize(metrics,('mode','rule','evaluation_family'))
    realizations,_=summarize(metrics,('mode','rule','realization'))
    failures=[r for r in metrics if not yes(r,'overall_confidence_pass')]
    count_failures=[r for r in metrics if not yes(r,'overall_count_pass')]
    effects=envelope_effects(out,test);effect_summary=[]
    information=information_summary(pairs)
    for mode in MODES:
        cells=[r for r in effects if r['mode']==mode];record=dict(mode=mode,cells=len(cells))
        for label in ('correct','wrong','indeterminate','true_retained'):
            delta=[r[label+'_difference'] for r in cells]
            record[label+'_increased_cells']=sum(x>1e-14 for x in delta)
            record[label+'_unchanged_cells']=sum(abs(x)<=1e-14 for x in delta)
            record[label+'_decreased_cells']=sum(x<-1e-14 for x in delta)
            record[label+'_minimum_difference']=min(delta);record[label+'_maximum_difference']=max(delta)
            record[label+'_paired_gained']=sum(r[label+'_gained'] for r in cells)
            record[label+'_paired_lost']=sum(r[label+'_lost'] for r in cells)
        record['double_rejected_to_wrong_total']=sum(r['baseline_double_rejected_to_wrong'] for r in cells)
        effect_summary.append(record)
    for name,rows in [('method_summary.csv',groups),('family_summary.csv',strata),('realization_summary.csv',realizations),
        ('worst_cells.csv',extremes),('confidence_failures.csv',failures),('count_failures.csv',count_failures),
        ('envelope_effects.csv',effects),('envelope_summary.csv',effect_summary),('information_summary.csv',information)]:csvwrite(out/name,rows,list(metrics[0]) if name in ('confidence_failures.csv','count_failures.csv') else None)
    library_rows=[]
    for r in execution['library_support']:
        g,role,realization,component=r['key']
        library_rows.append(dict(gravity=g,role=role,realization=realization,component=component,**{k:v for k,v in r.items() if k!='key'}))
    support=read(out/'template_manifest.json')
    if support['template_sha256']!=sha(out/'templates.npz') or support['execution_protocol_sha256']!=sha(out/'EXECUTION_PROTOCOL.json'):raise ValueError('Support provenance mismatch')
    csvwrite(out/'library_support.csv',library_rows);csvwrite(out/'fitting_component_support.csv',support['component_support'])
    support_summary=dict(minimum_effective_rows=min(r['effective_rows'] for r in library_rows),
        maximum_parent_weight=max(r['max_weight'] for r in library_rows),library_beta_rows=len(library_rows),
        fitted_component_rows=len(support['component_support']))
    durations={}
    for name in ('freeze','calibration','test'):
        path=out.parent/(name+'.log')
        if path.exists():
            matches=re.findall(r'COMPLETE\s+(\w+)\s+seconds\s+([\d.]+)',path.read_text(encoding='utf-8-sig'))
            if matches:durations[name]=float(matches[-1][1])
    artifact=dict(complete=True,scientific_status='software_only' if execution['fixture'] else 'conditional_development_experiment',
        methods=groups,by_family=strata,by_realization=realizations,envelope=effect_summary,information=information,support=support_summary,
        metric_rows=len(metrics),paired_rows=len(pairs),physical_cells=540,confidence_failure_rows=len(failures),count_failure_rows=len(count_failures),
        durations_seconds=durations,summary_source_sha256=sha(Path(__file__)),
        execution_protocol_sha256=sha(out/'EXECUTION_PROTOCOL.json'),test_manifest_sha256=sha(out/'test/manifest.json'),
        input_sha256={name:sha(out/name) for name in ('metrics.csv','paired.csv')},
        scope='Descriptive reconstruction from saved immutable results. Cellwise marginal intervals; no universal or Gaia guarantee.')
    (out/'report_artifacts.json').write_text(json.dumps(artifact,indent=2)+'\n',encoding='utf-8')
    display=[dict(Mode=r['mode'],Rule=r['rule'],Cells=r['cells'],Classification_count=r['classification_count_pass'],Retention_count=r['retention_count_pass'],
        Classification_CI=r['classification_confidence_pass'],Retention_CI=r['retention_confidence_pass'],Both_CI=r['overall_confidence_pass'],
        Worst_correct=f"{100*r['worst_correct']:.1f}%",Worst_wrong=f"{100*r['worst_wrong']:.1f}%",Worst_true_retained=f"{100*r['worst_true_retained']:.1f}%") for r in groups]
    paired_display=[dict(Mode=r['mode'],Correct_increased=r['correct_increased_cells'],Correct_decreased=r['correct_decreased_cells'],
        Wrong_increased=r['wrong_increased_cells'],Wrong_decreased=r['wrong_decreased_cells'],Retention_increased=r['true_retained_increased_cells'],
        Correct_delta_range=f"[{100*r['correct_minimum_difference']:.2f}, {100*r['correct_maximum_difference']:.2f}] pp",
        Wrong_delta_range=f"[{100*r['wrong_minimum_difference']:.2f}, {100*r['wrong_maximum_difference']:.2f}] pp") for r in effect_summary]
    info_display=[dict(Rule=r['first_rule'],Comparison='full66 minus '+r['second_mode'],
        Correct_increased=r['correct_increased_cells'],Correct_decreased=r['correct_decreased_cells'],
        Correct_delta_range=f"[{100*r['correct_minimum_difference']:.2f}, {100*r['correct_maximum_difference']:.2f}] pp",
        Wrong_delta_range=f"[{100*r['wrong_minimum_difference']:.2f}, {100*r['wrong_maximum_difference']:.2f}] pp") for r in information]
    text=['# Conditional information-block and calibration-library experiment','',
        f"The completed execution used {cal['catalogue_views']:,} paired calibration catalogue views and {test['catalogue_views']:,} test views, each evaluated in five fixed information variants. Each test cell contains {execution['repeats']:,} catalogue draws of {execution['N']:,} parents. Baseline and envelope rules act on those same scores and catalogues. The10 method/rule rows below each summarize540 cells; they are not10 independent physical experiments.",'',
        'All previous releases remain unchanged. This protocol was fixed after inspection of the earlier development experiment. It expands the beta grid prospectively, retains the original fitting banks, adds two calibration-library realizations and reserves two fresh test-library realizations. It is not a blind human study, a preregistration or a validation of the actual Gaia sky.','',
        '## Separate classification and retention criteria','',table(display),'',
        'Classification count-pass requires observed correct≥95% and wrong≤1%. Retention count-pass separately requires true-model retention≥99%. The stronger CI gates require Wilson95% lower(correct)≥95%, upper(wrong)≤1%, and lower(true-retained)≥99%, with classification and retention kept separate. Both-CI requires both gates. These are marginal descriptive criteria; they are not simultaneous confidence statements over all cells or selected extrema.','',
        f"There are {len(count_failures):,} method/rule/cell rows failing at least one count criterion and {len(failures):,} failing at least one stronger confidence criterion. Every failure is retained in count_failures.csv or confidence_failures.csv. metrics.csv contains all5,400 rows; family_summary.csv separates fitted-family, outside-inner-prior and outside-beta-grid tests; realization_summary.csv preserves the two test realizations. worst_cells.csv supplies the exact parameters and intervals behind each reported extreme. A good average never replaces a failed cell.",'',
        '## Paired effect of the calibration-library envelope','',table(paired_display),'',
        'An increased wrong rate is adverse, even though its signed change is positive. All rates above compare exactly paired catalogues; envelope_effects.csv retains all2,700 per-cell/mode contrasts, recovered/lost correct and wrong decisions, true-model retention changes and transitions from double rejection. paired.csv also retains all prescribed information-block comparisons under both calibration rules. The tables do not treat the three cadences as independent replicates.','',
        'The envelope is the maximum across libraries of each already Bonferroni-combined hypothesis p-value. It never lowers an individual p-value, so lost true-model retention must be zero; this identity was checked in the descriptive reconstruction. That does not imply wrong classifications can never increase: an earlier double rejection may become a single wrong retention. The observed paired transitions, rather than an assumed superiority, determine the interpretation.','',
        '## What the information comparisons identify','',table(info_display),'',
        'Each contrast retains all540 paired cells; information_summary.csv includes positive, zero and negative differences plus first-only/second-only discordant totals for correct, wrong and indeterminate outcomes. A positive wrong-rate difference is adverse. These summaries do not replace the exact-cell paired.csv or transform reused cadences into independent observations.','',
        'full66 versus no_diagnostics66 removes exactly the diagnostic block while retaining the same projected velocity and uncertainty/separation cells. full34 instead uses its own measured PM/covariance and only34-month source residual/acceleration diagnostics; its interface never reads66-month or cross-release information. That contrast concerns the whole later information package, not observation duration alone, and full34 is not a recreation of the public Gaia DR3 products.','',
        'covariates_only66 retains normalized separation and formal66-month precision. Separation is a physical outcome and selection variable. The control therefore measures information supplied by the stipulated parent and selection; successful classification by that control is not itself label leakage. conditional66 subtracts the marginal-z deviance for each nuisance point before profiling, using the full mixed distribution. It removes the direct marginal likelihood term while retaining differences in available stratum counts and conditional populations. Its calibration is not exactly conditional on each possible observed n_z, and it does not prove that parent assumptions have been eliminated.','',
        'The fitted grid is beta1/1.3/1.6 and the comp40/comp80 mixture family. The new1.15/1.45 cases are outside that discrete grid, and actual60/160 AU inner priors remain physical transfer tests. Pure systems are counted once per beta/hypothesis/realization and have no meaningful inner-prior or mixture-eta truth. Eta is null in test records; templates use eta only to mix40/80.','',
        'The f=0.15 points and eta=0.5 mixtures enter fitting/profiling and calibration, but have no dedicated test-library cells. Calibration coverage of those nuisance points is not an independent transfer test of them.','',
        '## Finite-library and physical limits','',
        'Rank calibration applies under its calibrating empirical law. The library envelope covers a finite union of those laws, not every fresh empirical library or continuous nuisance value. The two new test realizations measure transfer and provide a limited sensitivity check. They do not measure fitting-library variability: the original six fitting banks stay fixed. Repeated catalogues reuse finite banks; counts, effective row numbers and smoothing do not certify a larger sample size or a future-survey guarantee.','',
        f"Across the44 libraries and five stored beta weight sets, the minimum effective row count is {support_summary['minimum_effective_rows']:.1f} and the largest normalized parent weight is {support_summary['maximum_parent_weight']:.6g}. These are global weight-concentration diagnostics, not binwise support guarantees or counts of independent repeated catalogues. library_support.csv and fitting_component_support.csv retain the source records; every metric row separately records the mean fraction entering cells empty in the raw fitting union. Near-best nuisance counts use a descriptive2-deviance window and are not parameter confidence intervals.",'',
        'The same scalar two-mass Newton/QUMOND model, fixed synthetic external field, restricted stationary parent, empirical covariates, photocentre prescription and observation operator are inherited. True total mass is tied to a fixed photometric anchor before orbital scaling. The inherited empirical pool already reflects chance-alignment, RUWE, precision, RV and transverse-velocity selection. Source-level error/photometry relations, distance/parallax uncertainties and real-survey selection remain conditional approximations. The unresolved companion enters the outer interaction through a monopole approximation; this is not a complete QUMOND three-body solution.','',
        'All five variants use measured fields only. Gravity/companion labels and energy weights define supervised simulated laws and evaluation truth, never the individual measured feature vector. The numerical results require an independent implementation check; that review is a separate artifact and is not external human peer review. The public-sky readiness status stays false.','',
        '## Reproducibility','',
        'STATISTICAL_DECLARATION.json fixes sources, protocol,414 catalogue streams and all cases before execution. EXECUTION_PROTOCOL.json adds immutable input hashes and cached measured codes. Every calibration/test case saves one parent-index array shared by all methods/cadences, fullfloat64 nuisance surfaces, scores and decisions. Counts are reconstructible from those indices and codes and have canonical hashes; they are not duplicated on disk. The independent checker can reconstruct every case.','',
        f"Recorded stage durations (seconds, descriptive single-thread measurements): {json.dumps(durations,sort_keys=True)}. Code hashes, data hashes and any reproduction-run timing differences must be retained rather than edited into an old manifest.",'']
    if execution['fixture']:text.insert(2,'SOFTWARE FIXTURE ONLY. Every input is artificial; these numbers provide no physical gravity evidence.\n')
    (out/'REPORT.md').write_text('\n'.join(text),encoding='utf-8')
    print(json.dumps(dict(complete=True,metric_rows=len(metrics),group_rows=len(groups),envelope_rows=len(effects),confidence_failure_rows=len(failures))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=HERE/'results');run(p.parse_args().out)
