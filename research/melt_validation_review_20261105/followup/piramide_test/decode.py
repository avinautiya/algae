"""Narrow, validated reader for this archive's MATLAB5 MCOS table layout.

No MATLAB code is executed. Unsupported layouts are refused.
"""
import io
import numpy as np
import scipy.io
from scipy.io.matlab._mio5 import varmats_from_mat

def objects(raw):
    root=scipy.io.loadmat(io.BytesIO(raw),struct_as_record=False,squeeze_me=True)
    ws=root['__function_workspace__'].tobytes()
    assert ws[:4]==b'\x00\x01IM'
    parts=varmats_from_mat(io.BytesIO(raw[:124]+ws[:4]+ws[8:]))
    assert len(parts)==2
    sub=scipy.io.loadmat(parts[0][1],struct_as_record=False,squeeze_me=True)
    arr=sub['__function_workspace__'][0,0].MCOS['arr'][0]
    return root,arr

def table(raw):
    root,arr=objects(raw);props=arr[3];time=np.asarray(arr[2],dtype=float)
    assert hasattr(props,'varNames') and props.numRows==len(time)
    assert np.isfinite(time).all() and (np.diff(time)>0).all()
    cols={str(k):np.asarray(v,dtype=float) for k,v in zip(props.varNames,props.data)}
    assert all(v.shape==time.shape for v in cols.values())
    return time,cols

def raw_weather(raw):
    root,arr=objects(raw)
    # Explicitly verify this object is the first/raw Pir table, not the filled tables.
    reference=root['Hour'].Pir['arr'][0]
    assert int(reference[4])==1
    props=arr[3];time=np.asarray(arr[2],dtype=float)
    assert props.numRows==len(time) and (np.diff(time)>0).all()
    cols={str(k):np.asarray(v,dtype=float) for k,v in zip(props.varNames,props.data)}
    assert {'Ta','SW_in','SW_out'}.issubset(cols)
    assert all(v.shape==time.shape for v in cols.values())
    return time,cols
