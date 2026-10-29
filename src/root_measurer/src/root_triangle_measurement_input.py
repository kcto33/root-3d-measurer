"""Read immutable triangle PLY geometry without repair, welding or voxelization."""
from dataclasses import dataclass
import hashlib
import io
import math
import os
from pathlib import Path
import struct
import numpy as np

TYPES={'char':'b','int8':'b','uchar':'B','uint8':'B','short':'h','int16':'h',
       'ushort':'H','uint16':'H','int':'i','int32':'i','uint':'I','uint32':'I',
       'float':'f','float32':'f','double':'d','float64':'d'}
INTEGERS={'b','B','h','H','i','I'}


@dataclass(frozen=True)
class TriangleInput:
    vertices: np.ndarray
    faces: np.ndarray
    source_sha256: str
    metadata: dict


def read_triangle_ply_bytes(data,expected_sha=None,max_vertices=10000000,max_faces=20000000):
    sha=hashlib.sha256(data).hexdigest()
    if expected_sha is not None and sha!=expected_sha:raise ValueError('SOURCE_SHA_MISMATCH')
    stream=io.BytesIO(data)
    if stream.readline().strip()!=b'ply':raise ValueError('PLY_REQUIRED')
    elements=[];current=None;encoding=None
    for _ in range(256):
        raw=stream.readline()
        if not raw or stream.tell()>65536:raise ValueError('PLY_HEADER_GUARD')
        p=raw.decode('ascii').split()
        if not p:continue
        if p[0]=='end_header':break
        if p[0] in ('comment','obj_info'):continue
        if p[0]=='format':
            if len(p)!=3 or p[2]!='1.0':raise ValueError('PLY_VERSION')
            encoding=p[1]
        elif p[0]=='element':
            if len(p)!=3 or p[1] not in ('vertex','face'):raise ValueError('UNSUPPORTED_PLY_ELEMENT')
            current=dict(name=p[1],count=int(p[2]),properties=[]);elements.append(current)
        elif p[0]=='property':
            if current is None:raise ValueError('PROPERTY_WITHOUT_ELEMENT')
            current['properties'].append(p[1:])
        else:raise ValueError('UNSUPPORTED_PLY_HEADER')
    else:raise ValueError('PLY_HEADER_GUARD')
    if [e['name'] for e in elements]!=['vertex','face']:raise ValueError('TRIANGLE_SURFACE_REQUIRED_NOT_POINT_CLOUD')
    if encoding not in ('ascii','binary_little_endian','binary_big_endian'):raise ValueError('PLY_ENCODING')
    nv,nf=[e['count'] for e in elements]
    if not 0<nv<=max_vertices or not 0<nf<=max_faces:raise ValueError('PLY_COUNT_GUARD')
    for element in elements:
        names=[]
        for p in element['properties']:
            if p[0]=='list':
                if len(p)!=4 or element['name']!='face' or p[3] not in ('vertex_indices','vertex_index'):
                    raise ValueError('UNSUPPORTED_LIST_PROPERTY')
                if p[1] not in TYPES or p[2] not in TYPES or TYPES[p[1]] not in INTEGERS or TYPES[p[2]] not in INTEGERS:
                    raise ValueError('INTEGER_FACE_LIST_REQUIRED')
            elif len(p)!=2 or p[0] not in TYPES:raise ValueError('UNSUPPORTED_SCALAR_PROPERTY')
            names.append(p[-1])
        if len(set(names))!=len(names):raise ValueError('DUPLICATE_PROPERTY')
    vp=elements[0]['properties'];fp=elements[1]['properties']
    vnames=[p[-1] for p in vp]
    if not all(x in vnames for x in ['x','y','z']):raise ValueError('XYZ_REQUIRED')
    if sum(p[0]=='list' for p in fp)!=1:raise ValueError('ONE_FACE_INDEX_LIST_REQUIRED')
    endian='<' if encoding=='binary_little_endian' else '>'
    def exact(n):
        value=stream.read(n)
        if len(value)!=n:raise ValueError('TRUNCATED_PLY_BODY')
        return value
    def scalar(kind):
        fmt=endian+TYPES[kind];return struct.unpack(fmt,exact(struct.calcsize(fmt)))[0]
    vertices=np.empty((nv,3),float);faces=np.empty((nf,3),np.int64)
    if encoding=='ascii':
        for i in range(nv):
            tokens=stream.readline().decode('ascii').split()
            if len(tokens)!=len(vp):raise ValueError('ASCII_VERTEX_PROPERTY_COUNT')
            vertices[i]=[float(tokens[vnames.index(a)]) for a in ('x','y','z')]
    else:
        fmt=endian+''.join(TYPES[p[0]] for p in vp);size=struct.calcsize(fmt)
        for i in range(nv):
            row=struct.unpack(fmt,exact(size));vertices[i]=[row[vnames.index(a)] for a in ('x','y','z')]
    for i in range(nf):
        tokens=stream.readline().decode('ascii').split() if encoding=='ascii' else None
        cursor=0
        for prop in fp:
            if prop[0]=='list':
                count=int(tokens[cursor]) if tokens is not None else scalar(prop[1]);cursor+=1
                if count!=3:raise ValueError('NONTRIANGULAR_FACE_UNSUPPORTED_NO_TRIANGULATION')
                if tokens is not None:
                    if cursor+3>len(tokens):raise ValueError('TRUNCATED_FACE')
                    face=[int(v) for v in tokens[cursor:cursor+3]];cursor+=3
                else:face=[scalar(prop[2]) for _ in range(3)]
                faces[i]=face
            else:
                if tokens is not None:
                    if cursor>=len(tokens):raise ValueError('TRUNCATED_FACE')
                    cursor+=1
                else:scalar(prop[0])
        if tokens is not None and cursor!=len(tokens):raise ValueError('ASCII_FACE_PROPERTY_COUNT')
    if stream.read().strip():raise ValueError('UNPARSED_PLY_BODY')
    if not np.isfinite(vertices).all():raise ValueError('NONFINITE_COORDINATES')
    if faces.min()<0 or faces.max()>=nv:raise ValueError('FACE_INDEX_OUT_OF_RANGE')
    vertices.flags.writeable=False;faces.flags.writeable=False
    return TriangleInput(vertices,faces,sha,dict(format=encoding,vertex_properties=vnames,
        source_vertices=nv,source_faces=nf,coordinate_units='raw_source_units_not_inferred',
        geometry_repaired=False,voxelized=False,source_arrays_preserved=True))


def read_triangle_ply(path,expected_sha):
    path=Path(path)
    if os.name=='nt' and not str(path).startswith('\\\\?\\'):path=Path('\\\\?\\'+str(path.absolute()))
    data=path.read_bytes()
    result=read_triangle_ply_bytes(data,expected_sha)
    assert hashlib.sha256(path.read_bytes()).hexdigest()==expected_sha
    return result


def surface_measurements(mesh):
    v,f=mesh.vertices,mesh.faces
    repeated=int(np.count_nonzero((f[:,0]==f[:,1])|(f[:,1]==f[:,2])|(f[:,2]==f[:,0])))
    duplicate=len(f)-len(np.unique(np.sort(f,axis=1),axis=0))
    cross=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]])
    twice_area=np.linalg.norm(cross,axis=1);zero=int(np.count_nonzero(twice_area==0))
    edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
    edge,labels,counts=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    orientations=np.bincount(labels,weights=np.where(edges[:,0]<edges[:,1],1.,-1.),minlength=len(edge))
    boundary=int(np.count_nonzero(counts==1));nonmanifold=int(np.count_nonzero(counts>2))
    orientation_bad=int(np.count_nonzero((counts==2)&(orientations!=0)))
    closed_edge_consistent=not (boundary or nonmanifold or orientation_bad or duplicate or repeated or zero)
    used=np.unique(f);center=v[used].mean(axis=0);shift=v-center
    signed=math.fsum(np.einsum('ij,ij->i',shift[f[:,0]],np.cross(shift[f[:,1]],shift[f[:,2]])).tolist())/6
    return dict(status='SOURCE_TRIANGLE_GEOMETRY_DIAGNOSTICS',vertices=len(v),faces=len(f),
        referenced_vertices=len(used),unreferenced_vertices=len(v)-len(used),
        source_arrays_unchanged=True,surface_area=math.fsum(twice_area.tolist())/2,
        boundary_edges=boundary,nonmanifold_edges=nonmanifold,orientation_inconsistent_edges=orientation_bad,
        duplicate_faces=duplicate,repeated_index_faces=repeated,zero_area_faces=zero,
        closed_edge_consistent=closed_edge_consistent,declared_complex_euler=len(used)-len(edge)+len(f),
        oriented_enclosed_volume_integral=signed if closed_edge_consistent else None,
        absolute_net_volume_integral=abs(signed) if closed_edge_consistent else None,
        material_union_volume=None,vertex_manifold_qualified=False,self_intersections_checked=False,
        physical_volume_qualified=False,measurement_units='raw_source_units',
        limitation='closed oriented edge incidence is necessary, not certification of nonoverlap/material union or physical volume')
