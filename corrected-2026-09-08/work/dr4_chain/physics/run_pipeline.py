"""Resume the declared development/refinement/final validation sequence."""
import os
for _name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[_name]='1'
import json,hashlib
from pathlib import Path
from build_grid import run as build
from validate_model import run as validate
HERE=Path(__file__).resolve().parent

if __name__=='__main__':
    build(False)
    initial=HERE/'results'/'model_initial.npz'
    validate('development',initial)
    development=json.loads((HERE/'results'/'development_model_initial.json').read_text())
    selected=initial
    if development['refine_required']:
        build(True);selected=HERE/'results'/'model_refined.npz';validate('development',selected)
    choice=dict(selected_model=selected.name,model_sha256=hashlib.sha256(selected.read_bytes()).hexdigest(),
                initial_development_max_error=development['max_force_error'],refinement_trigger=.005,
                fixed_before_final=True,protocol_sha256=hashlib.sha256((HERE/'PROTOCOL.md').read_bytes()).hexdigest())
    (HERE/'results'/'grid_choice_before_final.json').write_text(json.dumps(choice,indent=2))
    print(json.dumps(choice),flush=True)
    validate('final',selected)
