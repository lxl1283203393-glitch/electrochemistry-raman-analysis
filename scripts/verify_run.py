"""Verify source hashes, raw-array/text fidelity and repeat-run numerical artifacts."""
from __future__ import annotations
import argparse
import hashlib
import json
import csv
import io
import re
from pathlib import Path
import sys
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from optosky_raman.storage import write_json


def verify(results_dir: Path,repeat: Path|None=None) -> dict:
    root=results_dir.resolve(); saved=json.loads((root/'results.json').read_text(encoding='utf-8'))
    checks={'sources':[],'figures':[],'repeat_arrays':[],'repeat_tables':[]}
    for record in saved['records']:
        source=Path(record['path']); content=source.read_bytes()
        text=content.decode(record['encoding'])
        snapshot=(root/record['snapshot']).read_bytes().decode('utf-8')
        raw_arrays_equal=None
        if record['arrays']:
            rows=list(csv.reader(io.StringIO(text),delimiter=record['metadata']['delimiter']))
            header_row=record['metadata']['data_header_row']-1
            header=[re.sub(r'[^a-z]','',v.lower()) for v in rows[header_row]]
            array_keys={'ramanshift':'x','raw':'channel_raw','dark':'channel_dark',
                'darksubtracted':'channel_dark_subtracted','baselinesubtracted':'channel_baseline_subtracted'}
            data=[row for row in rows[header_row+1:] if row and any(v.strip() for v in row)]
            with np.load(root/record['arrays']) as stored:
                raw_arrays_equal=all(np.array_equal(np.array([float(row[header.index(name)]) for row in data]),stored[key])
                    for name,key in array_keys.items())
        checks['sources'].append({'record_id':record['record_id'],
            'hash_match':hashlib.sha256(content).hexdigest()==record['sha256'],'snapshot_text_equal':snapshot==text,
            'raw_arrays_exact':raw_arrays_equal})
    for item in saved['figures']:
        if item['status']!='generated': continue
        for kind,path in {**item['paths'],'data':item['data_path']}.items():
            file=root/path
            checks['figures'].append({'category':item['category'],'format':kind,'exists_nonempty':file.exists() and file.stat().st_size>0})
    if repeat:
        again=json.loads((repeat/'results.json').read_text(encoding='utf-8'))
        other={r['record_id']:r for r in again['records']}
        for record in saved['records']:
            counterpart=other.get(record['record_id'])
            if not counterpart or not record['arrays'] or not counterpart['arrays']: continue
            with np.load(root/record['arrays']) as a,np.load(repeat/counterpart['arrays']) as b:
                equal=set(a.files)==set(b.files) and all(np.array_equal(a[k],b[k]) for k in a.files)
            checks['repeat_arrays'].append({'record_id':record['record_id'],'all_arrays_exact':equal})
        # Analyze and compare contain identical per-spectrum calculations; QC similarity
        # and per-batch summary tables intentionally depend on the batch collection.
        for name in ('peaks','preprocessing_sensitivity','detection_sensitivity','fit_components','fit_models'):
            a=pd.read_csv(root/'tables'/f'{name}.csv'); b=pd.read_csv(repeat/'tables'/f'{name}.csv')
            if 'record_id' in a and 'record_id' in b:
                b=b[b.record_id.isin(a.record_id)].reset_index(drop=True)
            equal=a.reset_index(drop=True).equals(b)
            checks['repeat_tables'].append({'table':name,'dataframe_exact':equal})
    all_bools=[v for rows in checks.values() for row in rows for v in row.values() if isinstance(v,bool)]
    checks['passed']=bool(all_bools) and all(all_bools)
    write_json(root/'verification.json',checks)
    return checks


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',required=True,type=Path)
    parser.add_argument('--compare',type=Path)
    args=parser.parse_args()
    result=verify(args.results,args.compare)
    print(json.dumps({'passed':result['passed'],**{k:len(v) for k,v in result.items() if isinstance(v,list)}}))
    raise SystemExit(0 if result['passed'] else 2)
