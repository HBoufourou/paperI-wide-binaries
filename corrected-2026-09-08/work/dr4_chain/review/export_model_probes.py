"""Export production API outputs for a separate independent comparison.

This exporter intentionally calls the production API. The independent checker
does not import it. require_validation=False is used only for numerical probes;
this program never generates scientific populations or authorizes a model.
"""
import argparse,hashlib,json
from pathlib import Path
import sys
import numpy as np


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--model',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    physics=args.model.resolve().parent.parent
    sys.path.insert(0,str(physics))
    from interaction_model import InteractionModel
    model=InteractionModel.load(args.model,require_validation=False)
    rng=np.random.default_rng(268900101);n=1200
    x=rng.normal(size=(n,3));x/=np.linalg.norm(x,axis=1)[:,None]
    x*=np.exp(rng.uniform(np.log(.03)+1e-6,np.log(30)-1e-6,n))[:,None]
    q=rng.choice([.1,.3,1.],n)
    vectors=np.array([[0.,0.,.03],[0.,0.,-.03],[0.,0.,30.],[0.,0.,-30.],
                      [.03,0.,0.],[30.,0.,0.]])
    cases=[(x,q), (x[:120].reshape(12,10,3),q[:12,None]),
           (x[0],.3),(x[0],np.array([.1,.3,1.])),
           (vectors,np.array([.1,.3,1.,1.,.3,.1]))]
    meta={'purpose':'numerical_api_probes_only','require_validation':False,
          'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'production_source_sha256':hashlib.sha256((physics/'interaction_model.py').read_bytes()).hexdigest(),
          'model_sha256':hashlib.sha256(args.model.read_bytes()).hexdigest(),
          'cases':[],'invalid_inputs':[]}
    arrays={}
    for gravity in ('newton','qumond'):
        for position,qq in cases:
            k=f'case{len(meta["cases"]):02d}'
            p=model.potential(position,qq,gravity);a=model.acceleration(position,qq,gravity)
            for field,value in [('position',position),('q',qq),('potential',p),('acceleration',a)]:arrays[k+'_'+field]=np.asarray(value)
            meta['cases'].append({'key':k,'gravity':gravity,'position_shape':list(np.shape(position)),
                                  'q_shape':list(np.shape(qq)),'potential_shape':list(np.shape(p)),
                                  'acceleration_shape':list(np.shape(a))})
    for label,pos,qq,g in [('zero_radius',[0,0,0],.3,'qumond'),('too_small',[.029,0,0],.3,'qumond'),
                          ('too_large',[30.1,0,0],.3,'qumond'),('unknown_q',[1,0,0],.4,'qumond'),
                          ('nan_position',[np.nan,0,0],.3,'qumond'),('unknown_gravity',[1,0,0],.3,'MOND')]:
        try:model.potential(pos,qq,g);rejected=False
        except ValueError:rejected=True
        meta['invalid_inputs'].append({'label':label,'rejected':rejected})
    arrays['metadata']=np.array(json.dumps(meta,sort_keys=True))
    args.out.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(args.out,**arrays)
    print(json.dumps(meta),flush=True)

if __name__=='__main__':main()
