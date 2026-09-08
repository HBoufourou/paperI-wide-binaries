"""Software-only invariance, algebra, leakage and calibration-order tests."""
import copy,unittest
import numpy as np
import engine as e
import features as f


def fixture(n=8):
    p=dict(distance_pc=np.full(n,100.),Mphot1=np.full(n,.7),Mphot2=np.full(n,.3),
        rproj_kau=np.tile([3.,4.],(n,1)),energy=np.full(n,-1.),true_gravity=np.full(n,999.))
    obs={}
    for b in (34,66):
        obs[f'pm{b}_masyr']=np.tile([.2,.3],(n,1))
        obs[f'cov{b}_masyr2']=np.tile(np.eye(2)*(.03 if b==34 else .02),(n,1,1))
        obs[f'pm_source{b}_masyr']=np.tile([[-.1,-.15],[.1,.15]],(n,1,1))
        obs[f'residual_score_source{b}']=np.zeros((n,2))
        obs[f'acceleration_score_source{b}']=np.ones((n,2))
    obs['cov_delta_source_masyr2']=np.tile(np.eye(2)*.01,(n,2,1,1))
    obs['cov_delta_relative_masyr2']=np.tile(np.eye(2)*.02,(n,1,1))
    return p,obs


def artificial_template():
    x=np.arange(2048)[None,None,None,:]
    h=np.arange(2)[None,:,None,None];t=np.arange(39)[None,None,:,None]
    p=1.+((x+7*t+3*h)%17)/20;p=p/p.sum(axis=-1,keepdims=True);p=np.tile(p,(3,1,1,1))
    return dict(full66=p,full34=p[...,::-1],no_diagnostics66=f.project_histogram(p,'no_diagnostics66'),
        covariates_only66=f.project_histogram(p,'covariates_only66'))


class StatsTests(unittest.TestCase):
    def test_units_and_source_relabel(self):
        p,o=fixture();x=f.measured_features(p,o,'full66')
        vn=np.sqrt(1.3271244e20/(5000*149597870700))/1000
        vel=np.array([.2,.3])*4.740470463533349*.1/vn
        self.assertAlmostEqual(x['parallel'][0],vel@np.array([.6,.8]))
        self.assertAlmostEqual(x['perpendicular'][0],vel[1]*.6-vel[0]*.8)
        q=copy.deepcopy(p);a=copy.deepcopy(o);q['rproj_kau']*=-1
        q['Mphot1'],q['Mphot2']=q['Mphot2'],q['Mphot1']
        for k in a:
            if '_source' in k:a[k]=a[k][:,::-1]
            elif k in ('pm34_masyr','pm66_masyr'):a[k]*=-1
        np.testing.assert_array_equal(f.measured_codes(q,a,'full66'),f.measured_codes(p,o,'full66'))

    def test_full34_never_accesses_66_or_difference(self):
        p,o=fixture();expected=f.measured_codes(p,o,'full34')
        class Guard(dict):
            def __getitem__(self,k):
                if '66' in k or 'delta' in k or 'cross' in k:raise AssertionError('Forbidden future access: '+k)
                return super().__getitem__(k)
        a=Guard({k:v for k,v in o.items() if k in f.OBS_FIELDS['full34']})
        p['energy'][:]=np.nan;p['true_gravity'][:]=np.nan
        np.testing.assert_array_equal(expected,f.measured_codes(p,a,'full34'))
        self.assertTrue(all('34' in k for k in f.OBS_FIELDS['full34']))

    def test_minimal_ablation_interfaces(self):
        p,o=fixture()
        full=f.measured_codes(p,o,'full66')
        for mode in ('no_diagnostics66','covariates_only66'):
            minimal={k:o[k] for k in f.OBS_FIELDS[mode]}
            np.testing.assert_array_equal(f.measured_codes(p,minimal,mode),f.project_codes(full,mode))

    def test_invalid_shapes_covariances(self):
        p,o=fixture()
        for mode in f.MODES:
            b='34' if mode=='full34' else '66';a=copy.deepcopy(o);a['cov'+b+'_masyr2'][0,0,1]=.01
            with self.assertRaises(ValueError):f.measured_codes(p,a,mode)
            a=copy.deepcopy(o);a['cov'+b+'_masyr2']=a['cov'+b+'_masyr2'][:,:1]
            with self.assertRaises(ValueError):f.measured_codes(p,a,mode)
        a=copy.deepcopy(o);a['residual_score_source34']=np.zeros((8,))
        with self.assertRaises(ValueError):f.measured_codes(p,a,'full34')

    def test_marginal_projection_and_smoothing_commute(self):
        rng=np.random.default_rng(310000001);codes=rng.integers(0,2048,701);weights=rng.uniform(.2,4,701)
        p,_,_=e.component_probability(codes,weights,2048)
        for mode,bins in [('no_diagnostics66',512),('covariates_only66',8)]:
            q,_,_=e.component_probability(f.project_codes(codes,mode),weights,bins)
            np.testing.assert_allclose(f.project_histogram(p,mode),q,rtol=1e-14,atol=2e-16)
            counts=np.bincount(codes,minlength=2048)
            np.testing.assert_array_equal(f.project_histogram(counts,mode),np.bincount(f.project_codes(codes,mode),minlength=bins))

    def test_conditional_surface_matches_direct_formula(self):
        rng=np.random.default_rng(310000002);template=artificial_template()
        c=np.array([[rng.multinomial(1000,np.full(2048,1/2048)) for _ in range(3)] for _ in range(3)],np.uint16)
        counts=dict(full66=c,full34=c[...,::-1],no_diagnostics66=f.project_histogram(c,'no_diagnostics66'),covariates_only66=f.project_histogram(c,'covariates_only66'))
        stats=e.profile(counts,template);z=f.project_codes(np.arange(2048),'covariates_only66')
        for si in range(3):
            for ri in range(3):
                cc=c[si,ri].astype(float);nz=np.bincount(z,weights=cc,minlength=8);positive=cc>0
                for hi,ti in [(0,0),(1,1),(1,38)]:
                    pp=template['full66'][si,hi,ti];pz=np.bincount(z,weights=pp,minlength=8)
                    expected=2*np.sum(cc[positive]*np.log(cc[positive]/(nz[z[positive]]*pp[positive]/pz[z[positive]])))
                    self.assertAlmostEqual(stats['surface'][4,si,ri,hi,ti],expected,places=9)
        np.testing.assert_array_equal(stats['best'][4],stats['surface'][4].argmin(axis=-1))

    def test_conditioning_after_mixture_is_required(self):
        # Different z masses make mix(conditionals) differ from condition(mix).
        p0=np.array([[.81,.09],[.05,.05]]);p1=np.array([[.02,.08],[.18,.72]])
        mix=.7*p0+.3*p1;correct=mix/mix.sum(axis=1,keepdims=True)
        wrong=.7*p0/p0.sum(axis=1,keepdims=True)+.3*p1/p1.sum(axis=1,keepdims=True)
        self.assertGreater(np.max(abs(correct-wrong)),.1)

    def test_envelope_after_combination_not_before(self):
        stats=dict(D=np.full((5,3,1,2),5.),R=np.full((5,3,1,2),5.))
        reference=np.ones((5,3,3,2,39,4,2))
        reference[:,0,...,1]=10;reference[:,1,...,0]=10;reference[:,2,...,1]=10
        a=e.calibrated(stats,reference)
        np.testing.assert_allclose(a['p_library'],.4)
        np.testing.assert_allclose(a['p_joint'],.4)
        incorrect=np.minimum(1,2*a['p_marginal'].max(axis=1).min(axis=-1))
        np.testing.assert_allclose(incorrect,1.)
        self.assertTrue(np.all(a['p_joint'][1]>=a['p_joint'][0]))

    def test_ties_and_direct_rank(self):
        stats=dict(D=np.full((5,3,2,2),2.),R=np.full((5,3,2,2),3.))
        refs=np.broadcast_to(np.array([1.,2.,2.,3.])[None,None,None,None,None,:,None],(5,3,3,2,39,4,2)).copy()
        a=e.calibrated(stats,refs)
        np.testing.assert_allclose(a['p_marginal'][...,0],4/5)
        np.testing.assert_allclose(a['p_marginal'][...,1],2/5)
        np.testing.assert_allclose(a['p_joint'],.8)

    def test_declared_cases_and_namespace(self):
        cal=list(e.calibration_cases());test=list(e.test_cases());s=[x['seed'] for x in cal+test]
        self.assertEqual(len(cal),234);self.assertEqual(len(test),180);self.assertEqual(len(set(s)),414)
        self.assertTrue(min(s)>=320000000 and max(s)<327000000)
        self.assertEqual(len(e.required_records()),44)
        for r in (1,2):
            for h in e.GRAVITIES:
                cases=[x for x in test if x['realization']==r and x['gravity']==h]
                self.assertEqual(len(cases),45)
                self.assertEqual(sum(x['f']==0 for x in cases),5)
                self.assertTrue(all(x['prior'] is None for x in cases if x['f']==0))

    def test_shared_indices_and_raw_hashes(self):
        rng=np.random.default_rng(310000003);banks={}
        for c in ('pure','comp40'):
            code=rng.integers(0,2048,(3,21),dtype=np.uint16)
            banks['Newton','test',1,c]=dict(codes={'full66':code,'full34':code[:,::-1]},weights={1.:np.arange(1,22)},N=21)
        case=dict(gravity='Newton',role='test',realization=1,beta=1.,f=.3,eta=None,prior=40,seed=310000004)
        counts,index,lineage=e.sample(banks,case,4)
        joined=np.concatenate([banks['Newton','test',1,c]['codes']['full66'] for c in ('pure','comp40')],axis=1)
        for si in range(3):
            for ri in range(4):np.testing.assert_array_equal(counts['full66'][si,ri],np.bincount(joined[si,index[ri]],minlength=2048))
        np.testing.assert_array_equal(counts['covariates_only66'],f.project_histogram(counts['full66'],'covariates_only66'))
        self.assertEqual(e.array_hash(counts['full66']),e.array_hash(counts['full66'].copy()))
        self.assertNotEqual(e.array_hash(counts['full66']),e.array_hash(counts['full66'].astype(np.uint32)))

    def test_classification_and_retention_distinct(self):
        lo,hi=e.wilson(990,1000);self.assertLess(lo,.99)
        self.assertLess(e.wilson(950,1000)[0],.95)
        self.assertGreater(e.wilson(10,1000)[1],.01)
        a=e.decisions(np.array([[[.005,.005],[.02,.005],[.005,.02],[.02,.02]]]))
        np.testing.assert_array_equal(a['decision'],[[-1,0,1,-1]])
        np.testing.assert_array_equal(a['reason'],[[1,2,2,0]])

    def test_public_api_unsupported_and_34_whitelist(self):
        p,o=fixture(1000);template=artificial_template()
        for mode,bins in zip(f.MODES,f.BINS):template['raw_'+mode]=np.ones((3,bins),bool)
        refs=np.zeros((5,3,3,2,39,2,2))
        minimal={k:o[k] for k in f.OBS_FIELDS['full34']}
        a=e.predict_catalogue(p,minimal,'uniform_white',template,refs,mode='full34',declared_scope_confirmed=True)
        self.assertEqual(a['status'],'conditional');self.assertEqual(a['decision'],'indeterminate')
        for pp,oo,confirmed in [(p,minimal,False),(p,{},True),({'distance_pc':np.ones(3)},minimal,True)]:
            a=e.predict_catalogue(pp,oo,'uniform_white',template,refs,mode='full34',declared_scope_confirmed=confirmed)
            self.assertEqual(a['status'],'unsupported');self.assertEqual(a['decision'],'indeterminate')


if __name__=='__main__':unittest.main()
