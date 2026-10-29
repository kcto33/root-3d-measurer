"""Exact center-sampling of closed oriented triangle chains, without mesh repair.

The solid is the UNION of each chain's nonzero winding region, including its
boundary. Attribute seam IDs remain distinct. This is a declared geometric
interpretation, not a selfintersection or physical-material certification.
"""
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import math
import numpy as np


def _cross2(a, b):
    return a[0]*b[1]-a[1]*b[0]


def _orient(a, b, p):
    return _cross2((b[0]-a[0], b[1]-a[1]), (p[0]-a[0], p[1]-a[1]))


def _ceil(a, b):
    return -((-a)//b)


def encode_dyadic_mesh(vertices, pitch, origin=None):
    v=np.asarray(vertices)
    if v.ndim!=2 or v.shape[1]!=3 or not len(v) or not np.isfinite(v).all():
        raise ValueError('FINITE_XYZ_REQUIRED')
    if isinstance(pitch, (bool,np.bool_)) or not math.isfinite(pitch) or pitch<=0:
        raise ValueError('POSITIVE_FINITE_PITCH_REQUIRED')
    pairs=[float(x).as_integer_ratio() for x in v.flat]
    origin=[0.,0.,0.] if origin is None else origin
    if len(origin)!=3 or not all(math.isfinite(float(x)) for x in origin):raise ValueError('FINITE_GRID_ORIGIN_REQUIRED')
    shifts=[float(x).as_integer_ratio() for x in origin]
    pn,pd=float(pitch).as_integer_ratio()
    power=max([pd.bit_length()-1]+[d.bit_length()-1 for n,d in pairs+shifts])+1
    ticks=[n << (power-(d.bit_length()-1)) for n,d in pairs]
    offset=[n << (power-(d.bit_length()-1)) for n,d in shifts]
    points=[tuple(x-y for x,y in zip(ticks[i:i+3],offset)) for i in range(0,len(ticks),3)]
    half=(pn << (power-(pd.bit_length()-1)))//2
    assert half>0
    return points,half,power


def _validate(v, faces, groups):
    f=np.asarray(faces)
    if f.ndim!=2 or f.shape[1]!=3 or not len(f) or f.dtype.kind not in 'iu':
        raise ValueError('INTEGER_TRIANGLES_REQUIRED')
    if f.min()<0 or f.max()>=len(v):raise ValueError('FACE_INDEX_OUT_OF_RANGE')
    if groups is None:groups=[(0,len(f))]
    start=0;proof=[]
    for a,b in groups:
        if a!=start or not a<b<=len(f):raise ValueError('EXACT_FACE_PARTITION_REQUIRED')
        counts=Counter();sign=Counter()
        for row in f[a:b]:
            p,q,r=[v[int(i)] for i in row]
            u=tuple(q[i]-p[i] for i in range(3));w=tuple(r[i]-p[i] for i in range(3))
            if not any((u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0])):
                raise ValueError('DEGENERATE_FACE')
            for x,y in ((p,q),(q,r),(r,p)):
                key=(x,y) if x<y else (y,x);counts[key]+=1;sign[key]+=1 if x<y else -1
        if any(n!=2 for n in counts.values()) or any(sign.values()):
            raise ValueError('GEOMETRIC_ORIENTED_EDGE_PAIR_CLOSURE_REQUIRED')
        proof.append(dict(face_range=[a,b],geometric_edge_pairs=len(counts),
                          geometric_oriented_boundary=0,source_IDs_merged=False))
        start=b
    if start!=len(f):raise ValueError('EXACT_FACE_PARTITION_REQUIRED')
    return f,list(groups),proof


def _vertical_span(tri, x, y):
    p,q,r=tri
    # The projected triangle is a line. A 3D nondegenerate triangle has a
    # nonconstant projected x or y; intersect its (t,z) triangle exactly.
    axis=0 if len({p[0],q[0],r[0]})>1 else 1
    t=(x,y)[axis]
    a=min(tri,key=lambda v:v[axis]);b=max(tri,key=lambda v:v[axis])
    if _orient(a,b,(x,y))!=0 or not a[axis]<=t<=b[axis]:return None
    zs=[]
    for u,v in ((p,q),(q,r),(r,p)):
        if u[axis]==v[axis]:
            if t==u[axis]:zs.extend((Fraction(u[2]),Fraction(v[2])))
        elif min(u[axis],v[axis])<=t<=max(u[axis],v[axis]):
            zs.append(Fraction(u[2])+Fraction(t-u[axis],v[axis]-u[axis])*(v[2]-u[2]))
    return (min(zs),max(zs)) if zs else None


def voxelize_closed_chains(vertices, faces, pitch, face_ranges=None, *,
                          max_xy_tests=20000000,max_output_cells=5000000,
                          max_abs_grid_index=2**30,max_abs_winding=None,grid_origin=None):
    """Return canonical lower-corner indices; source arrays are never modified."""
    source_v=np.asarray(vertices);source_f=np.asarray(faces)
    before=hashlib.sha256(source_v.tobytes()+source_f.tobytes()).hexdigest()
    v,half,power=encode_dyadic_mesh(source_v,pitch,grid_origin)
    f,groups,proof=_validate(v,source_f,face_ranges)
    stride=2*half;tests=0;cells=set();chain_rows=[]
    for chain,(a,b) in enumerate(groups):
        events=defaultdict(lambda:defaultdict(int));surfaces=defaultdict(list)
        for row in f[a:b]:
            tri=[v[int(i)] for i in row];area=_orient(*tri)
            xlo=_ceil(min(p[0] for p in tri)-half,stride)
            xhi=(max(p[0] for p in tri)-half)//stride
            ylo=_ceil(min(p[1] for p in tri)-half,stride)
            yhi=(max(p[1] for p in tri)-half)//stride
            if max(map(abs,(xlo,xhi,ylo,yhi)))>max_abs_grid_index:raise ValueError('GRID_INDEX_GUARD')
            tests+=max(0,xhi-xlo+1)*max(0,yhi-ylo+1)
            if tests>max_xy_tests:raise ValueError('TRIANGLE_XY_TEST_GUARD')
            sign=1 if area>0 else -1
            if area<0:tri=[tri[0],tri[2],tri[1]];area=-area
            for i in range(xlo,xhi+1):
                x=i*stride+half
                for j in range(ylo,yhi+1):
                    y=j*stride+half;key=(i,j)
                    if not area:
                        span=_vertical_span(tri,x,y)
                        if span is not None:surfaces[key].append(span)
                        continue
                    ws=[_orient(tri[1],tri[2],(x,y)),_orient(tri[2],tri[0],(x,y)),_orient(tri[0],tri[1],(x,y))]
                    if min(ws)<0:continue
                    z=Fraction(sum(w*p[2] for w,p in zip(ws,tri)),area)
                    surfaces[key].append((z,z))
                    owned=True
                    for w,u,t in zip(ws,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])):
                        if w==0:
                            dx=t[0]-u[0];dy=t[1]-u[1]
                            if not (dy>0 or (dy==0 and dx<0)):owned=False;break
                    if owned:events[key][z]+=-sign
        chain_count=0
        for key in surfaces.keys()|events.keys():
            spans=list(surfaces[key]);winding=0;previous=None
            for z,delta in sorted(events[key].items()):
                if winding and previous is not None:spans.append((previous,z))
                winding+=delta;previous=z
                if max_abs_winding is not None and abs(winding)>max_abs_winding:
                    raise ValueError('SAMPLED_CHAIN_NONUNIT_WINDING')
            if winding:raise ValueError('UNBALANCED_ORIENTED_RAY')
            ranges=[]
            for z0,z1 in spans:
                low=_ceil(z0.numerator-z0.denominator*half,z0.denominator*stride)
                high=(z1.numerator-z1.denominator*half)//(z1.denominator*stride)
                if low<=high:ranges.append((low,high))
            end=None
            for low,high in sorted(ranges):
                if max(abs(low),abs(high))>max_abs_grid_index:raise ValueError('GRID_INDEX_GUARD')
                low=max(low,end+1) if end is not None else low
                if low>high:continue
                if high-low+1>max_output_cells:raise ValueError('OUTPUT_CELL_GUARD')
                for k in range(low,high+1):cells.add((key[0],key[1],k));chain_count+=1
                if len(cells)>max_output_cells:raise ValueError('OUTPUT_CELL_GUARD')
                end=high
        chain_rows.append(dict(chain_id=chain,closed_sampled_cells=chain_count,**proof[chain]))
    result=np.asarray(sorted(cells),dtype=np.int64).reshape(-1,3);result.flags.writeable=False
    assert hashlib.sha256(source_v.tobytes()+source_f.tobytes()).hexdigest()==before
    return result,dict(status='COMPLETE_EXACT_CLOSED_CHAIN_CENTER_SAMPLING',
        pitch=float(pitch),dyadic_coordinate_power=power,xy_triangle_tests=tests,
        occupied_cells=len(result),chains=chain_rows,source_arrays_sha256=before,
        source_ID_arrays_unchanged=True,closed_boundary_included=True,
        solid_definition='union of declared chains nonzero winding regions including boundary',
        coordinate_rule='exact dyadic (index+1/2)*pitch; returned index is lower-corner cell index',
        welding_repair_unknown_fill_or_physical_admission=False)
