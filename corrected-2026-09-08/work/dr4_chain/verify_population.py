"""Synthetic Newton reference fixtures for population algebra and admissibility."""
from pathlib import Path
import hashlib,json
import numpy as np
from population import generate_library,inner_orbits

HERE=Path(__file__).resolve().parent

def main():
    potential=lambda x,q:-1/np.linalg.norm(x,axis=1)
    records=[]
    for j,component in enumerate(('pure','comp40','comp160')):
        p,meta=generate_library(256,268000011+j,potential,2.,component)
        radius=np.linalg.norm(p['position'],axis=1)
        def check(name,value):
            records.append(dict(component=component,name=name,passed=bool(value)))
            if not value:raise AssertionError(records[-1])
        check('finite arrays',all(np.isfinite(v).all() for v in p.values()))
        check('Newton baseline importance weight',np.allclose(p['weight_beta1'],1,rtol=1e-14,atol=1e-14))
        check('energy identity',np.max(abs(.5*np.sum(p['velocity']**2,axis=1)-1/radius-p['energy']))<1e-12)
        check('Lz identity',np.max(abs(np.cross(p['position'],p['velocity'])[:,2]-p['lz']))<1e-12)
        check('binding band',np.all((p['energy']>=-5)&(p['energy']<=-.1)))
        check('inner barrier',np.all(p['barrier_rlo']>.035))
        check('outer barrier',np.all(p['energy']<-2/30))
        check('projected selection',np.all((np.linalg.norm(p['rproj_kau'],axis=1)>=2)&(np.linalg.norm(p['rproj_kau'],axis=1)<=30)))
        check('beta reweight',np.allclose(p['weight_beta1p6']/p['weight_beta1'],(-p['energy'])**.6))
        check('true mass ratio',np.allclose(p['Mtrue_source'][:,1]/p['Mtrue_source'][:,0],p['q']))
        check('photometric anchor',np.allclose(p['Mphot1']+p['Mphot2'],p['Mphot_total_anchor'],rtol=1e-13,atol=1e-13))
        check('true mass flux conversion',np.allclose(p['Mtrue_total']*p['true_mass_flux_factor'],p['Mphot_total_anchor'],rtol=1e-13,atol=1e-13))
        check('photocentre separation identity',np.allclose(p['rproj_kau'],p['rproj_bary_kau']+p['initial_photocentre_offset_kau']))
        if component!='pure':
            check('hierarchy',np.all(p['inner_a_au']*(1+p['inner_e'])<p['barrier_rlo']*p['rM_au']/20))
            check('flux mass',np.all(p['Mphot1']+p['Mphot2']<=p['Mtrue_total']+1e-12))
        check('roles not science',meta['seed']<269000000)
    empty=inner_orbits(0,np.random.default_rng(268000020),np.empty((0,2)),np.empty(0),np.empty(0),40)
    assert empty[1].shape==(0,2) and len(empty[2])==0
    records.append(dict(name='empty companion proposal',passed=True))
    result=dict(passed=True,checks=len(records),records=records,software_fixture_only=True,
        scientific_libraries_generated=False,physical_potential_validated=False,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'population.py',HERE/'POPULATION_PROTOCOL.md',Path(__file__)]})
    (HERE/'POPULATION_IMPLEMENTATION_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=True,checks=len(records))))

if __name__=='__main__':main()
