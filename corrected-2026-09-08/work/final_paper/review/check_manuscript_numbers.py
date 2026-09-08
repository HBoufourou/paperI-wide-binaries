"""Bounded editorial number check; no production imports and no rescoring."""
from pathlib import Path
import csv,hashlib,json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'review'
MAN=ROOT/'MANUSCRIPT_REVIEW_VERSION.md'
RESULTS=ROOT/'stats/results'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def rows(path):
    with Path(path).open(encoding='utf-8',newline='') as stream:return list(csv.DictReader(stream))


def main():
    text=MAN.read_text(encoding='utf-8');manhash=sha(MAN);checks=[]
    inputs={str(MAN.relative_to(ROOT)):manhash}
    def load(name):
        p=RESULTS/name;inputs[str(p.relative_to(ROOT))]=sha(p)
        return rows(p) if p.suffix=='.csv' else read(p)
    def proof(relative):
        p=ROOT/relative;inputs[relative]=sha(p);return read(p)
    def check(name,observed,expected,source,tolerance=0):
        okay=abs(observed-expected)<=tolerance if isinstance(observed,(int,float)) and isinstance(expected,(int,float)) else observed==expected
        checks.append(dict(name=name,passed=bool(okay),observed=observed,expected=expected,tolerance=tolerance,source=source))
    def claim(anchor,values,source):
        # The anchor ties the numerical comparison to the exact manuscript claim.
        check('manuscript claim exists: '+anchor[:80],''.join(anchor.split()) in ''.join(text.split()),True,'MANUSCRIPT_REVIEW_VERSION.md; whitespace normalized')
        for label,observed,expected,tol in values:check(label,observed,expected,source,tol)

    metrics=load('metrics.csv');pairs=load('paired.csv');groups=load('method_summary.csv')
    info=load('information_summary.csv');env=load('envelope_summary.csv');effects=load('envelope_effects.csv')
    real=load('realization_summary.csv');support=load('library_support.csv');cal=load('calibration/manifest.json');test=load('test/manifest.json');exe=load('EXECUTION_PROTOCOL.json')
    lookup={(r['mode'],r['rule']):r for r in groups}
    mode_map={'Full66':'full66','NoDiagnostics66':'no_diagnostics66','Full34':'full34','CovariatesOnly66':'covariates_only66','Conditional66':'conditional66'}
    table_count={4:0,6:0}
    for line in text.splitlines():
        if not line.startswith('| '):continue
        cells=[v.strip() for v in line.strip('|').split('|')]
        if ' / ' not in cells[0]:continue
        mode,rule=cells[0].split(' / ')
        if mode not in mode_map:continue
        key=(mode_map[mode],'baseline' if rule=='base' else rule);r=lookup[key]
        if len(cells)==4:
            table_count[4]+=1
            for value,column in zip(cells[1:],('worst_correct','worst_wrong','worst_true_retained')):
                check('Table1 '+str(key)+' '+column,float(value),100*float(r[column]),'method_summary.csv',.0500000001)
        elif len(cells)==6:
            table_count[6]+=1
            for value,column in zip(cells[1:],('classification_count_pass','classification_confidence_pass','retention_count_pass','retention_confidence_pass','overall_confidence_pass')):
                check('Table2 '+str(key)+' '+column,int(value),int(r[column]),'method_summary.csv')
    check('Table1 covers every mode/rule',table_count[4],10,'manuscript table parsing')
    check('Table2 covers every mode/rule',table_count[6],10,'manuscript table parsing')

    full=lookup['full66','envelope'];base=lookup['full66','baseline'];nod=lookup['no_diagnostics66','envelope']
    claim('Two fresh test-library realizations provide 540 scenario-specific cells',[
        ('Abstract/two test realizations',len({r['realization'] for r in test['records']}),2,0),
        ('Abstract/540 cells',len(test['records'])*3,540,0),('Abstract/1000 catalogue repeats',exe['repeats'],1000,0),('Abstract/1000 parents',exe['N'],1000,0)],'EXECUTION_PROTOCOL/test manifest')
    check('Abstract/five variants',len({r['mode'] for r in groups}),5,'method_summary.csv')
    check('Abstract/three calibration libraries',len({r['realization'] for r in cal['records']}),3,'calibration/manifest.json')
    claim('targets hold in531 and 501 cells',[
        ('Abstract/531 Wilson classification',int(full['classification_confidence_pass']),531,0),
        ('Abstract/501 Wilson retention',int(full['retention_confidence_pass']),501,0),
        ('Abstract/98.1 minimum retention',float(full['worst_true_retained'])*100,98.1,.00000001),
        ('Abstract/540 no-diagnostic Wilson classification',int(nod['classification_confidence_pass']),540,0)],'method_summary.csv')
    claim('targets in all 540 cells', [('Abstract/540 Full66 count classification',int(full['classification_count_pass']),540,0)],'method_summary.csv')
    claim('recovers126 correct catalogue-view decisions but also creates15 wrong exclusive decisions',[
        ('Abstract/126 correct recovered',int(next(r for r in env if r['mode']=='full66')['correct_paired_gained']),126,0),
        ('Abstract/15 wrong created',int(next(r for r in env if r['mode']=='full66')['wrong_paired_gained']),15,0)],'envelope_summary.csv')

    physical=proof('review/physical_banks_verification.json');api=proof('review/delivery_api_verification.json')
    claim('All32 fresh libraries met the effective-size rule at 20,000 accepted parents',[
        ('5.1/new libraries',physical['fresh_banks'],32,0),('5.1/all library row counts',all(p['rows']==20000 for p in physical['populations']),True,0),
        ('5.1/fresh parents',physical['fresh_rows'],640000,0),('5.1/fresh scenario views',physical['fresh_rows']*3,1920000,0),
        ('5.1/historical inputs',physical['historical_banks'],12,0),('5.1/total libraries',len(physical['populations']),44,0),
        ('5.1/total parents',physical['total_accepted_rows'],880000,0),
        ('5.1/minimum ESS',min(float(r['effective_rows']) for r in support),16320.7,.05),
        ('5.1/max normalized weight',max(float(r['max_weight']) for r in support),.000193509,.0000000005)],'physical_banks_verification.json/library_support.csv')
    claim('passed6,143 grouped checks',[
        ('5.1/physical checks',physical['checks_count'],6143,0),('5.1/all observation files',physical['all_observation_files'],132,0),('5.1/physical complete/pass',physical['complete'] and physical['passed'],True,0)],'physical_banks_verification.json')
    claim('120 calls and 497 checks passed',[
        ('5.1/API catalogues',api['catalogues'],4,0),('5.1/API scenario views',api['scenario_views'],12,0),('5.1/API calls',api['method_rule_calls'],120,0),
        ('5.1/API checks',api['checks_count'],497,0),('5.1/API exact probabilities/decisions',api['all_p_values_and_decisions_exact'],True,0)],'delivery_api_verification.json')
    claim('observed99 per cent retention target passes in537 cells',[
        ('5.2/537 retention count',int(full['retention_count_pass']),537,0)],'method_summary.csv')

    missed=[r for r in metrics if r['mode']=='full66' and r['rule']=='envelope' and r['retention_count_pass']=='False']
    claim('three Full66-envelope cells missing the observed retention target',[
        ('5.2/three retention misses',len(missed),3,0),
        ('5.2/misses Newton beta1 realization2',all(r['gravity']=='Newton' and float(r['beta'])==1 and r['realization']=='2' for r in missed),True,0)],'metrics.csv')
    def cell(mode,gravity,beta,f,prior,scenario,realization):
        selected=[r for r in metrics if r['mode']==mode and r['rule']=='envelope' and r['gravity']==gravity and float(r['beta'])==beta and float(r['f'])==f and r['prior']==prior and r['scenario']==scenario and r['realization']==str(realization)]
        if len(selected)!=1:raise ValueError('Cell lookup is not unique')
        return selected[0]
    pure=cell('full66','Newton',1,0,'','sparse_floor',2)
    sparse=cell('full66','Newton',1,.3,'160','sparse_floor',2)
    uniform=cell('full66','Newton',1,.3,'160','uniform_white',2)
    for label,r,expected in [('pure sparse',pure,(986,10,4)),('prior160 sparse',sparse,(981,2,17))]:
        for col,value in zip(('true_retained','wrong','both_rejected'),expected):check('5.2/'+label+'/'+col,int(r[col]),value,'metrics.csv')
    check('5.2/prior160 uniform retention',int(uniform['true_retained']),987,'metrics.csv')
    for col,value in [('true_retained_wilson_lo',97.052),('true_retained_wilson_hi',98.780)]:check('5.2/retention CI/'+col,100*float(sparse[col]),value,'metrics.csv',.0005)
    for realization,ret,clas in [('1',257,266),('2',244,265)]:
        r=next(r for r in real if r['mode']=='full66' and r['rule']=='envelope' and r['realization']==realization)
        for col,value in [('cells',270),('retention_confidence_pass',ret),('classification_confidence_pass',clas)]:check('5.2/realization'+realization+'/'+col,int(r[col]),value,'realization_summary.csv')

    def contrast(other):return next(r for r in info if r['first_rule']=='envelope' and r['second_mode']==other)
    for other,expected in [('full34',(424,9,107,-.4,7.1)),('no_diagnostics66',(122,103,315,-1.5,1)),('covariates_only66',(540,0,0,42.8,100)),('conditional66',(318,11,211,-.3,9.6))]:
        r=contrast(other)
        for col,value in zip(('correct_increased_cells','correct_decreased_cells','correct_unchanged_cells'),expected[:3]):check('5.3/5.4 contrast '+other+'/'+col,int(r[col]),value,'information_summary.csv')
        for col,value in zip(('correct_minimum_difference','correct_maximum_difference'),expected[3:]):check('5.3/5.4 contrast '+other+'/'+col,100*float(r[col]),value,'information_summary.csv',1e-10)
    paired34=[r for r in pairs if r['comparison_type']=='information' and r['first_rule']=='envelope' and r['second_mode']=='full34']
    check('5.3/full66-full34 mean pp',100*sum(float(r['correct_difference']) for r in paired34)/len(paired34),1.0128,'paired.csv',.00005)
    means34=[sum(float(r['correct_difference']) for r in paired34 if r['realization']==j)/270 for j in ('1','2')]
    check('5.3/full66-full34 means positive in both realizations',all(x>0 for x in means34),True,'paired.csv')
    nodpairs=[r for r in pairs if r['comparison_type']=='information' and r['first_rule']=='envelope' and r['second_mode']=='no_diagnostics66']
    for col,value in [('correct_first_only',451),('correct_second_only',440),('wrong_first_only',148),('wrong_second_only',2)]:check('5.3/diagnostics/'+col,sum(int(r[col]) for r in nodpairs),value,'paired.csv')
    check('5.3/diagnostic net correct',sum(int(r['correct_first_only'])-int(r['correct_second_only']) for r in nodpairs),11,'paired.csv')
    nodmeans=[sum(float(r['correct_difference']) for r in nodpairs if r['realization']==j)/270 for j in ('1','2')]
    check('5.3/diagnostic mean opposite signs',nodmeans[0]*nodmeans[1]<0,True,'paired.csv')
    check('5.3/no-diagnostic Wilson retention',int(nod['retention_confidence_pass']),493,'method_summary.csv')
    check('5.4/covariate count pass neither rule',all(int(lookup['covariates_only66',r]['classification_count_pass'])==0 for r in ('baseline','envelope')),True,'method_summary.csv')
    cond=lookup['conditional66','envelope']
    for col,val in [('classification_count_pass',528),('classification_confidence_pass',508)]:check('5.4/conditional '+col,int(cond[col]),val,'method_summary.csv')
    cmin=cell('conditional66','QUMOND',1.6,.3,'40','clustered_correlated',1);cwrong=cell('conditional66','Newton',1,0,'','sparse_floor',2)
    check('5.4/conditional correct worst cell',float(cmin['correct_rate']),.901,'metrics.csv',1e-12)
    check('5.4/conditional wrong worst cell',float(cwrong['wrong_rate']),.021,'metrics.csv',1e-12)
    for col,val in [('wrong_wilson_lo',1.378),('wrong_wilson_hi',3.189)]:check('5.4/conditional wrong CI '+col,100*float(cwrong[col]),val,'metrics.csv',.0005)

    ef=next(r for r in env if r['mode']=='full66')
    for col,val in [('correct_increased_cells',61),('correct_decreased_cells',0),('correct_paired_gained',126),('correct_paired_lost',0),('wrong_paired_gained',15),('wrong_paired_lost',0),('wrong_increased_cells',15),('double_rejected_to_wrong_total',15)]:check('5.5/full66 envelope '+col,int(ef[col]),val,'envelope_summary.csv')
    poswrong=[r for r in effects if r['mode']=='full66' and int(r['wrong_gained'])]
    check('5.5/one extra wrong in each15 cells',len(poswrong)==15 and all(int(r['wrong_gained'])==1 and int(r['wrong_lost'])==0 and int(r['baseline_double_rejected_to_wrong'])==1 for r in poswrong),True,'envelope_effects.csv')
    for mode,correct_loss,wrong_loss in [('full34',850,100),('conditional66',41,57)]:
        r=next(r for r in env if r['mode']==mode)
        check('5.5/'+mode+' net correct lost',int(r['correct_paired_lost'])-int(r['correct_paired_gained']),correct_loss,'envelope_summary.csv')
        check('5.5/'+mode+' net wrong removed',int(r['wrong_paired_lost'])-int(r['wrong_paired_gained']),wrong_loss,'envelope_summary.csv')
    check('5.5/all individual true retention losses zero',all(int(r['true_retained_lost'])==0 for r in effects),True,'envelope_effects.csv')
    check('5.5/5400 metric rows',len(metrics),5400,'metrics.csv');check('5.5/7020 paired rows',len(pairs),7020,'paired.csv')
    claim('meets both stronger classification and retention targets in501',[
        ('Conclusion/501 combined gates',int(full['overall_confidence_pass']),501,0),
        ('Conclusion/no rule passes all combined gates',max(int(r['overall_confidence_pass']) for r in groups)<540,True,0)],'method_summary.csv')

    old=proof('../dr4_chain/review/population_verification.json')
    fixture_rows=sum(p['rows'] for p in old['populations'] if p['label'].startswith('trajectory_initial/'))
    check('Historical441024 includes1024 validation fixtures',fixture_rows,1024,'old population_verification.json')
    check('Historical scientific parents440000',old['total_accepted_rows']-fixture_rows,440000,'old population_verification.json')
    check('Current880000 =640000new+240000oldsubset',physical['total_accepted_rows'],physical['fresh_rows']+12*20000,'new physical_banks_verification.json')
    oldpath=ROOT.parent/'dr4_chain/decision/results/metrics.csv';inputs[str(oldpath.relative_to(ROOT.parent))]=sha(oldpath);oldmetrics=rows(oldpath)
    check('4.1/old216 cells',len(oldmetrics),216,'old metrics.csv')
    check('4.1/old208 Wilson classification',sum(r['confidence_target_passed']=='True' for r in oldmetrics),208,'old metrics.csv')
    fail=[r for r in oldmetrics if r['confidence_target_passed']=='False']
    check('4.1/old8 exceptions Newton beta1.3',len(fail)==8 and all(r['gravity']=='Newton' and float(r['beta'])==1.3 for r in fail),True,'old metrics.csv')
    oldfamily=[r for r in oldmetrics if r['family']=='endpoint' and r['beta_family']=='endpoint']
    check('4.1/old fitted-family minimum retention',min(float(r['true_retained_rate']) for r in oldfamily),.983,'old metrics.csv',1e-12)
    check('4.1/old22 scientific libraries',sum(not p['label'].startswith('trajectory_initial/') for p in old['populations']),22,'old population_verification.json')
    check('Current statistical count arithmetic',cal['catalogue_views']==702000 and test['catalogue_views']==540000 and len(exe['calibration_cases'])==234 and len(exe['test_cases'])==180,True,'stage manifests')
    check('Manuscript unchanged during check',sha(MAN),manhash,'SHA256 before/after')
    doc=dict(passed=all(r['passed'] for r in checks),checks_count=len(checks),manuscript_sha256=manhash,checks=checks,input_sha256=inputs,
        source_sha256=sha(__file__),scope='Editorial agreement of the abstract, all numerical results in section5, conclusion, and relevant design/old-new-count accounting with saved CSV/JSON. No production import or scientific reconstruction.',
        limitations='This checks recorded numerical claims and rounding, not grammar, publication suitability, bibliography, physical adequacy or the independent recalculation still owned by the separate verifier.')
    (OUT/'MANUSCRIPT_NUMBERS_CHECK.json').write_text(json.dumps(doc,indent=2)+'\n',encoding='utf-8')
    md=['# Manuscript numerical editorial check','',f"Status: {'PASS' if doc['passed'] else 'FAIL'}; {len(checks)} checks, {sum(not r['passed'] for r in checks)} discrepancies.",'',
        f'Manuscript SHA-256: `{manhash}`.',f'Checker SHA-256: `{doc["source_sha256"]}`.','',
        doc['scope'],'','No manuscript text, frozen source, protocol, score or decision was edited. The companion JSON records every comparison and all input hashes.','',
        'The two result tables were parsed and all80 numerical cells compared. Counts and transitions require exact integer agreement; displayed percentages/ranges are checked at their stated decimal precision. The four-decimal mean uses0.00005 percentage-point tolerance; three-decimal interval endpoints use0.0005 points. The effective size and maximum weight use half a last displayed decimal. Text anchors normalize whitespace only, to tolerate the concurrent editorial spacing correction. The first anchor-only failure is preserved in MANUSCRIPT_NUMBERS_CHECK_pre_spacing_normalization.json/.md; it contained no numerical discrepancy.','',
        'The declared501 combined full66-envelope passes equal257+244 across realizations. Its classification Wilson counts are266+265=531. Both retention and classification labels match their CSV columns. The diagnostic contrast has451 correct gains minus440 losses=11 net, but148 wrong gains minus2 losses=146 additional wrong views. Full34 envelope loses973−123=850 correct views net and removes106−6=100 wrong views net; conditional loses111−70=41 correct and removes57 wrong. All are paired catalogue views, not independent astrophysical systems.','',
        'Historical count reconciliation: the earlier independent population proof covers441024 rows =440000 scientific parents in22 libraries +1024 trajectory fixtures. The final experiment uses880000 scientific parents =640000 fresh +240000 from12 historical inputs; it does not add all22 earlier libraries or include the old fixtures.','',
        'The abstract and conclusion retain the distinction between observed95%/1% classification, marginal Wilson gates and99% true-model retention. No numerical discrepancy was found in the quoted extrema, identities, ranges, signed means or transitions. A physical or inference PASS is not inferred from this editorial check.','',
        '## Discrepancies','']
    errors=[r for r in checks if not r['passed']]
    md += [json.dumps(r,ensure_ascii=False) for r in errors] if errors else ['None.']
    md += ['','## Hashed inputs','']+[f'- `{name}`: `{h}`' for name,h in inputs.items()]
    (OUT/'MANUSCRIPT_NUMBERS_CHECK.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=doc['passed'],checks_count=len(checks),manuscript_sha256=manhash,errors=errors),indent=2))
    if not doc['passed']:raise SystemExit(1)


if __name__=='__main__':main()
