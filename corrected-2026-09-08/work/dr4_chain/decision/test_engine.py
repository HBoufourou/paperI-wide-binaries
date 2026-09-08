"""Software-only checks. No physical orbital library or old test is reused."""
import copy,unittest
import numpy as np
import engine as e

def fixture(n=10):
    p=dict(rproj_kau=np.tile([3.,4.],(n,1)),distance_pc=np.full(n,100.),Mphot1=np.full(n,1.),Mphot2=np.full(n,.3))
    o={k:np.ones((n,2)) for k in e.SOURCE_SCORES}
    o.update(pm34_masyr=np.tile([.1,.2],(n,1)),pm66_masyr=np.tile([.2,.3],(n,1)),
        pm_source34_masyr=np.tile([[-.05,-.1],[.05,.1]],(n,1,1)),pm_source66_masyr=np.tile([[-.1,-.15],[.1,.15]],(n,1,1)),
        cov66_masyr2=np.tile(np.eye(2)*.01,(n,1,1)),cov_delta_relative_masyr2=np.tile(np.eye(2)*.03,(n,1,1)),
        cov_delta_source_masyr2=np.tile(np.eye(2)*.015,(n,2,1,1)))
    return p,o

class EngineTests(unittest.TestCase):
    def test_unit_conversion_and_vector(self):
        p,o=fixture();x=e.measured_features(p,o)
        vc=np.sqrt(e.GM_SUN*1.3/(5000*e.AU))/1000
        v=np.array([.2,.3])*100*e.K/1000
        self.assertAlmostEqual(x[0,1],v@np.array([.6,.8])/vc,places=13)
        self.assertAlmostEqual(x[0,2],(v[1]*.6-v[0]*.8)/vc,places=13)
        self.assertAlmostEqual(x[0,0],5000*e.AU/np.sqrt(e.GM_SUN*1.3/e.A0),places=13)
    def test_source_relabel_and_hidden_fields(self):
        p,o=fixture();expected=e.measured_features(p,o);p['rproj_kau']*=-1;p['Mphot1'],p['Mphot2']=p['Mphot2'],p['Mphot1']
        for k in o:
            if '_source' in k:o[k]=o[k][:,::-1]
            elif k.startswith('pm'):o[k]*=-1
        p['gravity']=np.full(10,np.nan);o['true_velocity']=np.full((10,2),np.inf);p['has_companion']=np.ones(10)
        np.testing.assert_allclose(e.measured_features(p,o),expected,rtol=1e-13,atol=1e-13)
    def test_wrong_shapes_covariance_rejected(self):
        p,o=fixture()
        for k in ('cov66_masyr2','cov_delta_relative_masyr2','cov_delta_source_masyr2'):
            bad=copy.deepcopy(o);bad[k][...,0,1]+=1
            with self.assertRaises(ValueError):e.measured_features(p,bad)
        bad=copy.deepcopy(o);bad['pm_source34_masyr']=bad['pm_source34_masyr'][:,0]
        with self.assertRaises(ValueError):e.measured_features(p,bad)
    def test_joint_encoding_and_probabilities(self):
        x=np.array([[.1,-5,-5,0,.01],[4,5,5,1000,.1]])
        code=e.encode_features(x);np.testing.assert_array_equal(code,[0,2047])
        p,a=e.probabilities(code,np.array([1.,3.]),True)
        self.assertAlmostEqual(a['effective_rows'],1.6);self.assertAlmostEqual(p.sum(),1.);self.assertTrue((p>0).all())
        self.assertEqual(a['nonempty_bins'],2)
    def test_profile_against_scalar_deviance(self):
        rng=np.random.default_rng(279000001)
        p=rng.dirichlet(np.ones(e.NBIN),size=(3,2,len(e.THETA)))
        counts=np.array([[rng.multinomial(1000,np.full(e.NBIN,1/e.NBIN)) for _ in range(3)] for _ in range(3)])
        got=e.profile(counts,dict(probability=p))
        for si in range(3):
            for ri in range(3):
                c=counts[si,ri];active=c>0
                expected=np.array([[2*sum(float(n)*np.log(float(n)/(1000*float(pp))) for n,pp in zip(c[active],p[si,h,t,active])) for t in range(len(e.THETA))] for h in range(2)])
                np.testing.assert_allclose(expected,got['surface'][si,ri],rtol=1e-12,atol=1e-9)
        np.testing.assert_allclose(got['R'][:,:,0],-got['R'][:,:,1])
    def test_composite_ties_against_direct_rank(self):
        rng=np.random.default_rng(279000002)
        cal=np.sort(rng.integers(0,100,(3,2,len(e.THETA),1000,2)).astype(float),axis=3)
        stats={k:rng.integers(0,120,(3,5,2)).astype(float) for k in ('D','R')};result=e.composite_pvalues(stats,cal)
        for si in range(3):
            for hi in range(2):
                for ri in range(5):
                    vals=[]
                    for ai,k in enumerate(('D','R')):vals.append(max((1+np.count_nonzero(cal[si,hi,t,:,ai]>=stats[k][si,ri,hi]))/1001 for t in range(len(e.THETA))))
                    self.assertAlmostEqual(result['p_joint'][si,ri,hi],min(1.,2*min(vals)))
    def test_seed_namespaces_and_count_vs_confidence(self):
        records=list(e.calibration_cases())+list(e.test_cases());seeds=[r['seed'] for r in records]
        self.assertEqual(len(seeds),124);self.assertEqual(len(set(seeds)),124)
        self.assertTrue(min(seeds)>269999999 and max(seeds)<279000001)
        self.assertLess(e.wilson(950,1000)[0],.95);self.assertGreater(e.wilson(10,1000)[1],.01)
        self.assertLess(e.wilson(0,1000)[1],.01)
    def test_weighted_sampling_pairs_scenarios(self):
        banks={};codes=np.vstack((np.arange(8),np.arange(8)+10,np.arange(8)+20)).astype(np.uint16)
        for c in ('pure','comp40','comp80'):banks[('Newton','calibration',c)]=dict(N=8,codes=codes,weights={1.:np.arange(1,9,dtype=float)})
        counts,index,lineage=e.sample_counts(banks,'Newton','calibration',1.,.3,.5,279000003,repeats=4)
        self.assertEqual(counts.shape,(3,4,2048));self.assertTrue(np.all(counts.sum(axis=2)==1000))
        np.testing.assert_array_equal(counts[0,:,:8],counts[1,:,10:18]);np.testing.assert_array_equal(counts[0,:,:8],counts[2,:,20:28])
        again=e.sample_counts(banks,'Newton','calibration',1.,.3,.5,279000003,repeats=4)
        np.testing.assert_array_equal(index,again[1]);self.assertTrue(lineage['scenarios_paired'])

if __name__=='__main__':unittest.main()
