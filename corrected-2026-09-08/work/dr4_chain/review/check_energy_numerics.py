"""Declared auxiliary symmetry/tail checks, separate from the frozen gradient run."""
import hashlib,json,time
from pathlib import Path
import numpy as np
from independent_energy import potential

HERE=Path(__file__).resolve().parent

def main():
    specs=[dict(q=1.,theta_deg=45.),dict(q=1.,theta_deg=135.),
           dict(q=.3,theta_deg=135.),dict(q=1/.3,theta_deg=45.),
           dict(q=.3,theta_deg=45.,threshold=.0005),
           dict(q=.3,theta_deg=45.,radial_max=100.)]
    out={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'energy_source_sha256':hashlib.sha256((HERE/'independent_energy.py').read_bytes()).hexdigest(),
         'predeclared_tolerance_relative':1e-6,'cases':[],'complete':False}
    path=HERE/'energy_numerics.json'
    for spec in specs:
        t=time.time();p=potential(order=8,**spec);p['seconds']=time.time()-t
        out['cases'].append(p);path.write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(p),flush=True)
    values=[p['U'] for p in out['cases']]
    base=json.loads((HERE/'energy_pilot.json').read_text())['results'][-1]['U']
    out['checks']={
        'equal_mass_reflection_relative':abs(values[0]-values[1])/abs(values[0]),
        'unequal_mass_label_exchange_relative':abs(values[2]-values[3])/abs(values[2]),
        'mixed_hessian_switch_relative':abs(values[4]-base)/abs(base),
        'infinite_tail_breakpoint_relative':abs(values[5]-base)/abs(base)}
    out['passed']=max(out['checks'].values())<1e-6;out['complete']=True
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':out['passed'],'checks':out['checks']}),flush=True)

if __name__=='__main__': main()
