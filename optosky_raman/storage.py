"""JSON and table writers with explicit numerical/null encoding."""
from __future__ import annotations
from datetime import datetime
from pathlib import Path
from typing import Any
import json
import numpy as np
import pandas as pd


def serializable(value: Any) -> Any:
    if isinstance(value,dict): return {str(k):serializable(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [serializable(v) for v in value]
    if isinstance(value,np.ndarray): return serializable(value.tolist())
    if isinstance(value,np.generic): return serializable(value.item())
    if isinstance(value,float) and not np.isfinite(value): return None
    if isinstance(value,(datetime,Path)): return str(value)
    if value is pd.NA: return None
    return value


def write_json(path: Path,value: Any) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(serializable(value),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')


def write_table(path: Path,frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    if len(frame.columns)==0: frame=pd.DataFrame(columns=['status'])
    frame.to_csv(path,index=False,encoding='utf-8-sig',float_format='%.12g')
