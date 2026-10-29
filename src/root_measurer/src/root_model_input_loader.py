"""Shared mesh-only file input; complete native binary XYZ is readonly memmap."""
import hashlib,os
from pathlib import Path
import numpy as np
from root_triangle_measurement_input import TriangleInput,read_triangle_ply_bytes

def file_sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1048576),b''):h.update(chunk)
    return h.hexdigest()

def load_model_input(path,*,execution_domain='local',expected_sha=None):
    path=Path(path);native=execution_domain=='registered_slurm_complete_source'
    if execution_domain!='local' and (not native or not os.environ.get('SLURM_JOB_ID','').isdigit()):raise ValueError('REGISTERED_SLURM_BEFORE_COMPLETE_INPUT')
    size=path.stat().st_size
    if not 0<size<=268435456:raise ValueError('UNCHANGED_SOURCE_BYTE_GUARD')
    before=file_sha(path)
    if expected_sha is not None and before!=expected_sha:raise ValueError('SOURCE_SHA_MISMATCH')
    if not native:
        result=read_triangle_ply_bytes(path.read_bytes(),before,max_vertices=3000,max_faces=1000)
    else:
        lines=[]
        with path.open('rb') as stream:
            for _ in range(256):
                line=stream.readline();lines.append(line)
                if sum(map(len,lines))>65536 or not line:raise ValueError('PLY_HEADER_GUARD')
                if line.strip()==b'end_header':break
            else:raise ValueError('PLY_HEADER_GUARD')
        header=b''.join(lines);tokens=[line.decode('ascii').split() for line in lines if line.strip()]
        tokens=[t for t in tokens if t[0] not in ('comment','obj_info')]
        if len(tokens)!=9 or tokens[0]!=['ply'] or tokens[1]!=['format','binary_little_endian','1.0'] or tokens[-1]!=['end_header']:raise ValueError('NATIVE_READONLY_BINARY_XYZ_TRIANGLES_LAYOUT_REQUIRED')
        if tokens[2][:2]!=['element','vertex'] or tokens[6][:2]!=['element','face']:raise ValueError('NATIVE_READONLY_BINARY_XYZ_TRIANGLES_LAYOUT_REQUIRED')
        nv,nf=int(tokens[2][2]),int(tokens[6][2])
        if not 0<nv<=1200000 or not 0<nf<=2000000:raise ValueError('UNCHANGED_COMPLETE_SOURCE_COUNTS')
        if tokens[3:6]!=[['property','double',a] for a in ('x','y','z')]:raise ValueError('NATIVE_READONLY_BINARY_XYZ_TRIANGLES_LAYOUT_REQUIRED')
        prop=tokens[7]
        if len(prop)!=5 or prop[:3]!=['property','list','uchar'] or prop[3] not in ('int','uint') or prop[4] not in ('vertex_indices','vertex_index'):raise ValueError('NATIVE_READONLY_BINARY_XYZ_TRIANGLES_LAYOUT_REQUIRED')
        if size!=len(header)+24*nv+13*nf:raise ValueError('EXACT_BINARY_TRIANGLE_FILE_BYTES_REQUIRED')
        v=np.memmap(path,dtype='<f8',mode='r',offset=len(header),shape=(nv,3))
        dt=np.dtype([('count','u1'),('vertices','<i4' if prop[3]=='int' else '<u4',(3,))],align=False)
        table=np.memmap(path,dtype=dt,mode='r',offset=len(header)+24*nv,shape=(nf,));f=table['vertices']
        if not np.all(table['count']==3):raise ValueError('NONTRIANGULAR_FACE_UNSUPPORTED_NO_TRIANGULATION')
        if not np.isfinite(v).all() or f.min()<0 or f.max()>=nv:raise ValueError('FINITE_XYZ_AND_ORIGINAL_FACE_INDICES_REQUIRED')
        assert not v.flags.writeable and not f.flags.writeable
        result=TriangleInput(v,f,before,dict(format='binary_little_endian',source_vertices=nv,source_faces=nf,readonly_memmap=True,source_arrays_preserved=True,geometry_repaired=False))
    if path.stat().st_size!=size or file_sha(path)!=before:raise ValueError('ORIGINAL_SOURCE_FILE_CHANGED')
    return result
