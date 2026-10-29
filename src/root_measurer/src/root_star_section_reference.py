"""Independent sufficient reference for a star-shaped unit-winding plane disk.

This module is a scorer, never an input or fallback for the radius estimator.
It integrates oriented source-triangle cuts without constructing contour paths.
"""
from collections import Counter,defaultdict
from fractions import Fraction as F
import hashlib
import math
import numpy as np
from root_closed_triangle_voxelization import encode_dyadic_mesh,_validate


def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def orient(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def hull(points):
    points=sorted(set(points))
    if len(points)<3:return points
    def half(seq):
        out=[]
        for p in seq:
            while len(out)>=2 and orient(out[-2],out[-1],p)<=0:out.pop()
            out.append(p)
        return out
    return half(points)[:-1]+half(points[::-1])[:-1]


def strictly_separated(a,b):
    for poly in [a,b]:
        for u,v in zip(poly,poly[1:]+poly[:1]):
            n=(-(v[1]-u[1]),v[0]-u[0]);pa=[dot(p,n) for p in a];pb=[dot(p,n) for p in b]
            if max(pa)<min(pb) or max(pb)<min(pa):return True
    return False


class StarSectionReference:
    def __init__(self,mesh,face_ranges=None):
        self.mesh=mesh;self.before=hashlib.sha256(mesh.vertices.tobytes()+mesh.faces.tobytes()).hexdigest()
        self.vertices,_,self.power=encode_dyadic_mesh(mesh.vertices,1.);self.scale=1<<self.power
        faces,self.groups,self.closure=_validate(self.vertices,mesh.faces,face_ranges)
        self.triangles=[tuple(self.vertices[int(i)] for i in face) for face in faces]
        self.normals=[cross(sub(b,a),sub(c,a)) for a,b,c in self.triangles]

    def section(self,point,normal):
        point=np.asarray(point,float);normal=np.asarray(normal,float)
        if point.shape!=(3,) or normal.shape!=(3,) or not np.isfinite(point).all() or not np.isfinite(normal).all() or not np.any(normal):
            raise ValueError('FINITE_REFERENCE_QUERY_PLANE_REQUIRED')
        p=tuple(F(float(x))*self.scale for x in point);n=tuple(F(float(x)) for x in normal)
        drop=max(range(3),key=lambda i:abs(n[i]));keep=[i for i in range(3) if i!=drop]
        project=lambda q:tuple(q[k] for k in keep);q=project(p)
        base=dict(area=None,radius=None,reference_qualified=False,source_sha256=self.mesh.source_sha256,
            definition='oriented-cut integration with strict star-wedge, unit-winding and other-chain hull separation certificates',
            contour_path_reconstruction_or_measured_radius_used=False,formal_float_interval=False)
        groups=[]
        for group,(begin,end) in enumerate(self.groups):
            cuts=defaultdict(list)
            for face_id in range(begin,end):
                tri=self.triangles[face_id];ds=[dot(n,sub(v,p)) for v in tri]
                if all(x==0 for x in ds):return dict(base,status='REFERENCE_MISSING_COPLANAR_FACE',face_id=face_id)
                if min(ds)>0 or max(ds)<0:continue
                hits={tri[i] for i in range(3) if ds[i]==0}
                for i,j in [(0,1),(1,2),(2,0)]:
                    if ds[i]*ds[j]<0:
                        t=ds[i]/(ds[i]-ds[j]);hits.add(tuple(tri[i][k]+t*(tri[j][k]-tri[i][k]) for k in range(3)))
                if len(hits)==1:return dict(base,status='REFERENCE_MISSING_POINT_CONTACT',face_id=face_id)
                if not hits:continue
                if len(hits)!=2:return dict(base,status='REFERENCE_MISSING_UNRESOLVED_CUT')
                a,b=sorted(hits);tangent=cross(n,self.normals[face_id])
                direction=dot(sub(b,a),tangent)
                if not direction:return dict(base,status='REFERENCE_MISSING_CUT_ORIENTATION')
                oriented=(a,b) if direction>0 else (b,a)
                nonzero=[x for x in ds if x]
                cuts[(a,b)].append(dict(oriented=oriented,edge_in_plane=sum(x==0 for x in ds)==2,
                    third_sign=(1 if nonzero[0]>0 else -1) if len(nonzero)==1 else None))
            segments=[]
            for owners in cuts.values():
                if len(owners)>1:
                    if (len(owners)!=2 or not all(o['edge_in_plane'] for o in owners)
                        or owners[0]['third_sign']==owners[1]['third_sign']
                        or owners[0]['oriented']!=owners[1]['oriented']):
                        return dict(base,status='REFERENCE_MISSING_NONTRANSVERSE_EDGE_CONTACT')
                elif owners[0]['edge_in_plane']:return dict(base,status='REFERENCE_MISSING_OPEN_PLANE_EDGE')
                segments.append(owners[0]['oriented'])
            if not segments:continue
            balance=Counter();winding=0;wedges=[];vector=[F(0),F(0),F(0)];on_boundary=False
            for a,b in segments:
                balance[a]+=1;balance[b]-=1;aa,bb=project(a),project(b);w=orient(q,aa,bb);wedges.append(w)
                if w==0 and all(min(aa[k],bb[k])<=q[k]<=max(aa[k],bb[k]) for k in range(2)):on_boundary=True
                if aa[1]<=q[1]<bb[1] and w>0:winding+=1
                elif bb[1]<=q[1]<aa[1] and w<0:winding-=1
                area_vector=cross(sub(a,p),sub(b,p))
                for k in range(3):vector[k]+=area_vector[k]/(2*self.scale**2)
            if any(balance.values()):return dict(base,status='REFERENCE_MISSING_UNBALANCED_ORIENTED_CUTS')
            if on_boundary:return dict(base,status='REFERENCE_MISSING_QUERY_ON_BOUNDARY')
            star=all(w>0 for w in wedges) or all(w<0 for w in wedges)
            groups.append(dict(chain_id=group,cut_segments=len(segments),winding=winding,strict_star_wedges=star,
                hull=hull([project(v) for segment in segments for v in segment]),vector=vector))
        selected=[g for g in groups if abs(g['winding'])==1 and g['strict_star_wedges']]
        evidence=[dict(chain_id=g['chain_id'],cut_segments=g['cut_segments'],winding=g['winding'],strict_star_wedges=g['strict_star_wedges']) for g in groups]
        if len(selected)!=1:return dict(base,status='REFERENCE_MISSING_UNIQUE_STAR_UNIT_WINDING_CHAIN',chains=evidence)
        target=selected[0]
        for other in groups:
            if other is target:continue
            if len(target['hull'])<3 or len(other['hull'])<3 or not strictly_separated(target['hull'],other['hull']):
                return dict(base,status='REFERENCE_MISSING_OTHER_CHAIN_HULL_SEPARATION',chains=evidence)
        area=math.hypot(*(float(x) for x in target['vector']))
        if not math.isfinite(area) or area<=0:return dict(base,status='REFERENCE_MISSING_POSITIVE_FINITE_AREA')
        assert self.before==hashlib.sha256(self.mesh.vertices.tobytes()+self.mesh.faces.tobytes()).hexdigest()
        return dict(base,status='REFERENCE_COMPLETE_STRICT_STAR_UNIT_WINDING_DISK',reference_qualified=True,
            area=area,radius=math.sqrt(area/math.pi),chains=evidence,selected_chain_id=target['chain_id'],
            exact_oriented_area_vector=[dict(numerator=str(x.numerator),denominator=str(x.denominator)) for x in target['vector']],
            other_chain_sections_strictly_hull_separated=True,source_geometry_unchanged=True,
            final_norm_and_sqrt_rounded=True)
