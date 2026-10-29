"""Isolated local closed disks while preserving all open source components."""
import ast
import copy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path
from root_triangle_plane_section import _orient,_intersects
from root_local_triangle_section import _contains


def _hull(points):
    points=sorted(set(points))
    if len(points)<3:return points
    def half(seq):
        out=[]
        for p in seq:
            while len(out)>=2 and _orient(out[-2],out[-1],p)<=0:out.pop()
            out.append(p)
        return out
    return half(points)[:-1]+half(points[::-1])[:-1]


def _separated(a,b):
    for poly in [a,b]:
        for u,v in zip(poly,poly[1:]+poly[:1]):
            n=(-(v[1]-u[1]),v[0]-u[0])
            if n==(0,0):continue
            aa=[x*n[0]+y*n[1] for x,y in a];bb=[x*n[0]+y*n[1] for x,y in b]
            if max(aa)<min(bb) or max(bb)<min(aa):return True
    return False


def _boundary(poly,q):
    return any(_orient(a,b,q)==0 and all(min(a[k],b[k])<=q[k]<=max(a[k],b[k]) for k in [0,1])
        for a,b in zip(poly,poly[1:]+poly[:1]))


@lru_cache(maxsize=1)
def _intersection_class():
    path=Path(__file__).with_name('root_triangle_plane_section.py');raw=path.read_bytes();tree=ast.parse(raw)
    original=ast.dump(tree,include_attributes=False)
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='TrianglePlaneSectionMeasurer')
    method=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='section')
    at=next(i for i,n in enumerate(method.body) if isinstance(n,ast.If) and ast.unparse(n.test).startswith('any((len(v) != 2'))
    old=copy.deepcopy(method.body[at]);method.body[at]=ast.parse("return dict(base,status='COMPLETE_EXACT_PLANE_GRAPH',_adj=adj,_contacts=contacts,_segments=segments,_p=p,_n=n)").body[0]
    restored=copy.deepcopy(tree);rc=next(n for n in restored.body if isinstance(n,ast.ClassDef) and n.name==cls.name)
    rm=next(n for n in rc.body if isinstance(n,ast.FunctionDef) and n.name=='section');rm.body[at]=old
    assert ast.dump(restored,include_attributes=False)==original
    namespace={};exec(compile(ast.fix_missing_locations(tree),str(path),'exec'),namespace)
    kind=namespace['TrianglePlaneSectionMeasurer'];kind.original_source_sha256=hashlib.sha256(raw).hexdigest()
    return kind


class FragmentTriangleSectionMeasurer:
    def __init__(self,mesh):
        self.mesh=mesh;self.kernel=_intersection_class()(mesh);self.prepared=None
        self.before=hashlib.sha256(mesh.vertices.tobytes()+mesh.faces.tobytes()).hexdigest()

    def prepare(self,origin,normal):
        """Extract each original source face intersection once for this plane."""
        raw=self.kernel.section(origin,normal);self.prepared=raw
        if '_adj' not in raw:return dict(status=raw['status'],graph_available=False)
        p,n=raw['_p'],raw['_n'];drop=max(range(3),key=lambda k:abs(n[k]));self.keep=[k for k in range(3) if k!=drop]
        self.drop=drop;adj=raw['_adj'];unseen=set(adj);components=[]
        while unseen:
            start=min(unseen);nodes=set();stack=[start]
            while stack:
                node=stack.pop()
                if node in nodes:continue
                nodes.add(node);stack.extend(adj[node]-nodes)
            unseen-=nodes;path=None;simple=False
            if len(nodes)>=3 and all(len(adj[v])==2 for v in nodes):
                path=[min(nodes)];previous=None;current=path[0]
                while True:
                    nxt=min(adj[current]-({previous} if previous is not None else set()))
                    if nxt==path[0]:break
                    if nxt in path:path=None;break
                    path.append(nxt);previous,current=current,nxt
                if path is not None:
                    xy=[self.project(v) for v in path];simple=True
                    for i in range(len(xy)):
                        for j in range(i+1,len(xy)):
                            if j==i+1 or (i==0 and j==len(xy)-1):continue
                            if _intersects(xy[i],xy[(i+1)%len(xy)],xy[j],xy[(j+1)%len(xy)]):simple=False
            components.append(dict(component_id=len(components),nodes=sorted(nodes),path=path,simple_closed=simple,
                edges=sorted({tuple(sorted((v,w))) for v in nodes for w in adj[v]})))
        self.components=components
        return dict(status='COMPLETE_ALL_SOURCE_PLANE_COMPONENTS_RETAINED',graph_available=True,
            components=len(components),simple_closed_components=sum(c['simple_closed'] for c in components),
            open_or_unqualified_components=sum(not c['simple_closed'] for c in components),isolated_contacts=len(raw['_contacts']))

    def project(self,q):return tuple(q[k] for k in self.keep)

    def _world(self,q):
        raw=self.prepared;p,n=raw['_p'],raw['_n'];out=[F(0)]*3
        for k,v in zip(self.keep,q):out[k]=v
        out[self.drop]=p[self.drop]-sum(n[k]*(out[k]-p[k]) for k in self.keep)/n[self.drop]
        return tuple(out)

    def candidate_queries(self):
        if self.prepared is None or '_adj' not in self.prepared:return []
        out=[]
        for c in self.components:
            if not c['simple_closed']:continue
            poly=[self.project(v) for v in c['path']]
            cross=[a[0]*b[1]-a[1]*b[0] for a,b in zip(poly,poly[1:]+poly[:1])];area2=sum(cross)
            if not area2:continue
            q=tuple(sum((a[k]+b[k])*w for a,b,w in zip(poly,poly[1:]+poly[:1],cross))/(3*area2) for k in [0,1]);method='source_contour_area_centroid'
            if not _contains(poly,q) or _boundary(poly,q):
                q=None;method='first_exact_source_polygon_ear_centroid'
                for i in range(len(poly)):
                    a,b,d=poly[i-1],poly[i],poly[(i+1)%len(poly)]
                    if _orient(a,b,d)*area2<=0:continue
                    def inside(v):return all(w*area2>=0 for w in [_orient(a,b,v),_orient(b,d,v),_orient(d,a,v)])
                    if any(inside(v) for j,v in enumerate(poly) if j not in {(i-1)%len(poly),i,(i+1)%len(poly)}):continue
                    q=tuple((a[k]+b[k]+d[k])/3 for k in [0,1]);break
                if q is None:continue
            out.append(dict(component_id=c['component_id'],query=self._world(q),query_method=method))
        return out

    def section_in_prepared_plane(self,query):
        raw=self.prepared;base=dict(area=None,radius=None,source_sha256=self.mesh.source_sha256,
            full_source_surface_closed_or_material_volume_admitted=False,physical_root_radius=None,
            source_geometry_or_outside_components_deleted=False,original_kernel_prefix_inverse_AST_exact=True)
        if raw is None or '_adj' not in raw:return dict(base,status='MISSING_PREPARED_SOURCE_PLANE_GRAPH')
        point=tuple(F(v) for v in query)
        if sum(n*(x-p) for n,x,p in zip(raw['_n'],point,raw['_p']))!=0:raise ValueError('QUERY_MUST_LIE_ON_THE_EXACT_PREPARED_PLANE')
        q=self.project(point);selected=[c for c in self.components if c['simple_closed'] and _contains([self.project(v) for v in c['path']],q)]
        evidence=dict(total_source_components=len(self.components),outside_source_components=len(self.components)-len(selected),
            exact_query=[[v.numerator,v.denominator] for v in point],query_display=list(map(float,point)))
        if len(selected)!=1:return dict(base,**evidence,status='MISSING_UNIQUE_QUERY_CONTAINING_CLOSED_COMPONENT')
        c=selected[0];poly=[self.project(v) for v in c['path']];target=_hull(poly)
        if _boundary(poly,q):return dict(base,**evidence,status='MISSING_QUERY_STRICTLY_INSIDE_LOCAL_DISK')
        for other in self.components:
            if other is c:continue
            other_hull=_hull([self.project(v) for v in other['nodes']])
            if not _separated(target,other_hull):
                return dict(base,**evidence,status='MISSING_STRICT_LOCAL_DISK_OUTSIDE_COMPONENT_SEPARATION',other_component_id=other['component_id'])
        for contact in raw['_contacts']-set(raw['_adj']):
            if not _separated(target,[self.project(contact)]):return dict(base,**evidence,status='MISSING_ISOLATED_CONTACT_SEPARATION')
        vector=[F(0),F(0),F(0)];path=c['path']
        for a,b in zip(path,path[1:]+path[:1]):
            for k in range(3):vector[k]+=(a[(k+1)%3]*b[(k+2)%3]-a[(k+2)%3]*b[(k+1)%3])/2
        area=math.hypot(*(float(x) for x in vector))
        if not math.isfinite(area) or area<=0:return dict(base,**evidence,status='MISSING_POSITIVE_FINITE_LOCAL_AREA')
        assert hashlib.sha256(self.mesh.vertices.tobytes()+self.mesh.faces.tobytes()).hexdigest()==self.before
        return dict(base,**evidence,status='COMPLETE_STRICTLY_ISOLATED_LOCAL_SOURCE_DISK',area=area,radius=math.sqrt(area/math.pi),
            selected_component_id=c['component_id'],selected_contour_vertices=[list(map(float,v)) for v in path],
            exact_area_vector=[[v.numerator,v.denominator] for v in vector],outside_components_strict_hull_separation=True,
            radius_definition='sqrt(area/pi) of one strictly isolated local source-plane disk',final_norm_sqrt_rounded=True)

    def graph_receipt(self):
        raw=self.prepared
        if raw is None or '_adj' not in raw:return dict(status=None if raw is None else raw['status'])
        vertices=sorted(raw['_adj']);ids={v:i for i,v in enumerate(vertices)}
        return dict(exact_vertices=[[[x.numerator,x.denominator] for x in v] for v in vertices],
            segments=[dict(nodes=[ids[a],ids[b]],source_face_ids=[o['face_id'] for o in owners]) for (a,b),owners in raw['_segments'].items()],
            components=[dict(component_id=c['component_id'],node_ids=[ids[v] for v in c['nodes']],simple_closed=c['simple_closed']) for c in self.components],
            isolated_contacts=[[[x.numerator,x.denominator] for x in v] for v in sorted(raw['_contacts']-set(raw['_adj']))],
            all_source_components_retained=True,source_vertices_or_IDs_welded=False)
