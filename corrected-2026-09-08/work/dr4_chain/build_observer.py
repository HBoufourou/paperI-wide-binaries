"""Derive a joint trajectory observer from the frozen prior epoch operator."""
from pathlib import Path
import hashlib,json

HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'dr4_binary/companions/epoch_operator_extended.py'
EXPECTED='57d6b2f1ba761fae481650cb3cb483508f73cd644804ff527703fa4fdaa78ed2'

def main():
    assert hashlib.sha256(BASE.read_bytes()).hexdigest()==EXPECTED
    code=BASE.read_text(encoding='utf-8')
    changes=[
      ("scenario='uniform_white',seed=1,ra_deg=None,dec_deg=None,batch_size=384):",
       "scenario='uniform_white',seed=1,ra_deg=None,dec_deg=None,batch_size=384,\n                  outer_relative_xy_au=None,mass_fraction=None):"),
      ("names=['pm','bias_noise','inner_bias','noise']", "names=['pm','bias_noise','inner_bias','noise','outer']"),
      ("L=np.linalg.cholesky(C);rng=np.random.default_rng(seed)",
       """L=np.linalg.cholesky(C);rng=np.random.default_rng(seed)
    trajectory=np.zeros((n,len(t),2)) if outer_relative_xy_au is None else np.asarray(outer_relative_xy_au,float)
    fractions=np.full((n,2),.5) if mass_fraction is None else np.asarray(mass_fraction,float)
    if trajectory.shape!=(n,len(t),2) or not np.isfinite(trajectory).all():
        raise ValueError('outer trajectory must be finite (N,epochs,2) AU')
    if fractions.shape!=(n,2) or not np.isfinite(fractions).all() or np.any(fractions<=0) or not np.allclose(fractions.sum(axis=1),1,rtol=0,atol=1e-12):
        raise ValueError('positive true source mass fractions summing to one required')
    if outer_relative_xy_au is not None and wide_v_kms is not None:
        raise ValueError('Supply a trajectory or a straight velocity, never both')"""),
      ("xy=inner_positions(subinner,t,d[sl]);al=xy[:,:,:,0]*np.cos(angle)+xy[:,:,:,1]*np.sin(angle)",
       """xy=inner_positions(subinner,t,d[sl]);al=xy[:,:,:,0]*np.cos(angle)+xy[:,:,:,1]*np.sin(angle)
        relative=trajectory[sl]*1000/d[sl,None,None]
        coefficients=np.column_stack((-fractions[sl,1],fractions[sl,0]))
        wide_al=coefficients[:,:,None]*(relative[:,:,0]*np.cos(angle)+relative[:,:,1]*np.sin(angle))[:,None,:]"""),
      ("y=al[:,:,:m]+noise[:,:,:m]", "y=al[:,:,:m]+wide_al[:,:,:m]+noise[:,:,:m]"),
      ("fit5=np.einsum('npt,nst->nsp',g[5][1],yw)",
       """fit5=np.einsum('npt,nst->nsp',g[5][1],yw)
            ww=np.einsum('ij,nsj->nsi',g['invL'],wide_al[:,:,:m])
            fitw=np.einsum('npt,nst->nsp',g[5][1],ww)
            wb=fitw[:,1,2:4]-fitw[:,0,2:4]"""),
      ("out[f'bias_noise{b}_masyr'][sl]=offset", "out[f'bias_noise{b}_masyr'][sl]=offset-wb"),
      ("out[f'noise{b}_masyr'][sl]=offset-ib", "out[f'noise{b}_masyr'][sl]=offset-ib-wb\n            out[f'outer{b}_masyr'][sl]=truepm[sl]+wb"),
    ]
    for old,new in changes:
        assert code.count(old)==1,(old,code.count(old));code=code.replace(old,new)
    prefix='''# Derived from the frozen dr4_binary epoch operator; see OBSERVER_PROTOCOL.md.
# Both outer and inner deterministic signals enter the same GLS residuals.
'''
    target=HERE/'joint_observer.py';target.write_text(prefix+code,encoding='utf-8')
    result=dict(base_sha256=EXPECTED,derived_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),replacements=len(changes))
    (HERE/'OBSERVER_BUILD.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))

if __name__=='__main__':main()
