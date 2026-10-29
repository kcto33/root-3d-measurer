"""Exact binary-plane intersection contours of immutable triangle geometry.

Only one simple closed positive-area contour admits an area-equivalent radius.
No surface repair, global material union, coplanar filling or source-ID weld.
"""
from collections import defaultdict
from fractions import Fraction
import math
import numpy as np
from root_triangle_measurement_input import TriangleInput


def _orient(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def _intersects(a,b,c,d):
    if any(max(a[k],b[k])<min(c[k],d[k]) or max(c[k],d[k])<min(a[k],b[k]) for k in range(2)):return False
    aa,bb,cc,dd=_orient(a,b,c),_orient(a,b,d),_orient(c,d,a),_orient(c,d,b)
    return aa*bb<=0 and cc*dd<=0


class TrianglePlaneSectionMeasurer:
    def __init__(self,mesh,max_segments=10000):
        if not isinstance(mesh,TriangleInput):raise ValueError('QUALIFIED_TRIANGLE_INPUT_REQUIRED')
        self.mesh=mesh;self.vertices=[tuple(Fraction(float(x)) for x in row) for row in mesh.vertices]
        self.faces=mesh.faces;self.max_segments=max_segments

    def section(self,point,normal):
        p=np.asarray(point,float);n=np.asarray(normal,float)
        if p.shape!=(3,) or n.shape!=(3,) or not np.isfinite(p).all() or not np.isfinite(n).all() or not np.any(n):
            raise ValueError('FINITE_NONZERO_SECTION_PLANE_REQUIRED')
        p=tuple(Fraction(float(x)) for x in p);n=tuple(Fraction(float(x)) for x in n)
        distance=[sum(q*(x-y) for q,x,y in zip(n,v,p)) for v in self.vertices]
        segments=defaultdict(list);contacts=set()
        base=dict(area=None,radius=None,source_sha256=self.mesh.source_sha256,source_ID_arrays_modified=False,
            global_solid_or_union_admitted=False,physical_radius_qualified=False,
            geometry='exact binary64 query plane/triangle intersection; source coordinates unchanged')
        for face_id,face in enumerate(self.faces):
            ids=list(map(int,face));v=[self.vertices[i] for i in ids];ds=[distance[i] for i in ids]
            if all(d==0 for d in ds):return dict(base,status='MISSING_COPLANAR_TRIANGLE_PLANE',face_id=face_id)
            if all(d>0 for d in ds) or all(d<0 for d in ds):continue
            hits={v[i] for i in range(3) if ds[i]==0}
            for i,j in ((0,1),(1,2),(2,0)):
                if ds[i]*ds[j]<0:
                    t=ds[i]/(ds[i]-ds[j]);hits.add(tuple(v[i][k]+t*(v[j][k]-v[i][k]) for k in range(3)))
            if len(hits)==1:contacts.update(hits);continue
            if not hits:continue
            if len(hits)!=2:return dict(base,status='MISSING_UNRESOLVED_TRIANGLE_INTERSECTION',face_id=face_id)
            a,b=sorted(hits)
            third=[d for d in ds if d!=0]
            in_plane_edge=sum(d==0 for d in ds)==2
            segments[(a,b)].append(dict(face_id=face_id,in_plane_edge=in_plane_edge,third_sign=(1 if third[0]>0 else -1) if in_plane_edge else None))
            if len(segments)>self.max_segments:return dict(base,status='MISSING_SECTION_SEGMENT_GUARD')
        if not segments:return dict(base,status='NO_POSITIVE_TRIANGLE_SECTION',isolated_contacts=len(contacts))
        adj=defaultdict(set);edge_ownership=[]
        for (a,b),owners in segments.items():
            if len(owners)>1:
                if len(owners)!=2 or not all(o['in_plane_edge'] for o in owners):
                    return dict(base,status='MISSING_COINCIDENT_INTERSECTION_SEGMENTS')
                if owners[0]['third_sign']==owners[1]['third_sign']:
                    return dict(base,status='MISSING_TANGENT_EDGE_CONTACT')
                edge_ownership.append(dict(source_faces=[o['face_id'] for o in owners],opposite_plane_sides=True))
            elif owners[0]['in_plane_edge']:return dict(base,status='MISSING_OPEN_EDGE_COINCIDENCE')
            adj[a].add(b);adj[b].add(a)
        if any(len(v)!=2 for v in adj.values()):return dict(base,status='MISSING_NONMANIFOLD_SECTION_CONTACT',degree_counts=[len(v) for v in adj.values()])
        if contacts-set(adj):return dict(base,status='MISSING_ISOLATED_SECTION_CONTACT')
        unseen=set(adj);loops=[]
        while unseen:
            start=min(unseen);path=[start];previous=None;current=start
            while True:
                following=sorted(adj[current]-({previous} if previous is not None else set()))
                nxt=following[0]
                if nxt==start:break
                if nxt in path:return dict(base,status='MISSING_REPEATED_SECTION_CONTOUR')
                path.append(nxt);previous,current=current,nxt
            if len(path)<3:return dict(base,status='NO_POSITIVE_TRIANGLE_SECTION')
            unseen.difference_update(path);loops.append(path)
        if len(loops)!=1:return dict(base,status='MISSING_MULTIPLE_SECTION_CONTOURS',contours=len(loops),holes_or_union_not_filled=True)
        path=loops[0];axis=max(range(3),key=lambda i:abs(n[i]));other=[i for i in range(3) if i!=axis]
        xy=[tuple(v[i] for i in other) for v in path];count=len(path)
        for i in range(count):
            for j in range(i+1,count):
                if j==i+1 or (i==0 and j==count-1):continue
                if _intersects(xy[i],xy[(i+1)%count],xy[j],xy[(j+1)%count]):
                    return dict(base,status='MISSING_SELF_INTERSECTING_SECTION_CONTOUR')
        vector=[Fraction(0),Fraction(0),Fraction(0)]
        for a,b in zip(path,path[1:]+path[:1]):
            for k in range(3):vector[k]+=(a[(k+1)%3]*b[(k+2)%3]-a[(k+2)%3]*b[(k+1)%3])/2
        area=math.hypot(*(float(x) for x in vector))
        if not math.isfinite(area) or area<=0:return dict(base,status='MISSING_POSITIVE_FINITE_SECTION_AREA')
        return dict(base,status='COMPLETE_SIMPLE_CLOSED_TRIANGLE_PLANE_SECTION',area=area,radius=math.sqrt(area/math.pi),
            contour_vertices=[[float(x) for x in v] for v in path],contour_edges=count,
            exact_area_vector=[dict(numerator=str(x.numerator),denominator=str(x.denominator)) for x in vector],
            geometric_edge_plane_ownership=edge_ownership,
            radius_definition='sqrt(area/pi) of exactly intersected single simple plane contour',
            final_area_and_radius_sqrt_rounded=True)
