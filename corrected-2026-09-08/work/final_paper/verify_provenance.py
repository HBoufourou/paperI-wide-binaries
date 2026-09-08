"""Repeat the original-mock identity proof without writing historical outputs."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(out):
    sources=ROOT/'work/sources';repo=ROOT/'work/github/paperI-wide-binaries-af8ddf07/data'
    official=sources/'10652994_Newton_dr3_MSMS_d200pc_5.csv';mp=sources/'zenodo_10652994.json'
    meta=json.loads(mp.read_text(encoding='utf-8'));file=next(x for x in meta['files'] if x['key']=='Newton_dr3_MSMS_d200pc_5.csv')
    md5=hashlib.md5(official.read_bytes()).hexdigest();assert file['checksum']=='md5:'+md5
    ids=['source_id1','source_id2'];dtype={k:str for k in ids}
    a=pd.read_csv(official,dtype=dtype)
    bp=repo/'chae_colonnes_utiles.csv';cp=repo/'chae_complement.csv'
    b=pd.read_csv(bp,dtype=dtype).merge(pd.read_csv(cp,dtype=dtype),on=ids,validate='one_to_one')
    assert not a.duplicated(ids).any() and not b.duplicated(ids).any()
    aligned=a.set_index(ids).loc[pd.MultiIndex.from_frame(b[ids])].reset_index()
    result=[]
    for column in b:
        if column in ids:continue
        ok=np.isclose(b[column],aligned[column],rtol=1e-9,atol=0,equal_nan=True)
        result.append(dict(column=column,n_not_close=int(np.sum(~ok)),max_absolute_difference=float(np.nanmax(abs(b[column]-aligned[column])))))
    assert len(a)==len(b)==81088 and len(result)==36 and all(x['n_not_close']==0 for x in result)
    observed=sources/'10986733_gaia_dr3_MSMS_d200pc_ruwe.csv';om=sources/'zenodo_10986733.json'
    ometa=json.loads(om.read_text(encoding='utf-8'));of=next(x for x in ometa['files'] if x['key']=='gaia_dr3_MSMS_d200pc_ruwe.csv')
    observed_md5=hashlib.md5(observed.read_bytes()).hexdigest();assert of['checksum']=='md5:'+observed_md5
    observed_rows=len(pd.read_csv(observed,usecols=ids,dtype=dtype));assert observed_rows==81880
    record=dict(passed=True,official_rows=len(a),repository_rows=len(b),checked_non_id_columns=36,
        rtol=1e-9,atol=0,equal_nan=True,numerical_catalogue_identity=True,csv_byte_identity_claimed=False,
        mock_md5=md5,mock_sha256=sha(official),observed_rows=observed_rows,observed_md5=observed_md5,observed_sha256=sha(observed),
        source_sha256=sha(Path(__file__)),inputs={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in (official,mp,bp,cp,observed,om)},columns=result)
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k not in ('inputs','columns')}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=Path(__file__).resolve().parent/'review/provenance_verification.json');run(p.parse_args().out)
