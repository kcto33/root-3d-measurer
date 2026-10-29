"""Exact plane sections of reconstructed cubical solids, with missingness.

No surface fitting, triangle clipping heuristic, hole closing or observation
weights are introduced. Shared global grid-edge/vertex identities identify
polygon boundaries. Degenerate and non-disk sections are explicitly rejected
as radius estimates; their area may remain a geometric descriptor.
"""
from collections import defaultdict
import itertools
import math
import numpy as np

_CORNERS = np.array(list(itertools.product((0,1),repeat=3)),np.int64)
_EDGES = [(i,j,int(np.flatnonzero(_CORNERS[j]-_CORNERS[i])[0]))
          for i in range(8) for j in range(i+1,8)
          if np.sum(np.abs(_CORNERS[j]-_CORNERS[i]))==1]


class CubicalSectionMeasurer:
    def __init__(self, solid):
        self.solid=solid
        self.cells=solid.cells
        self.centers=self.cells.astype(float)+.5
        self.occupied=set(map(tuple,self.cells))

    def section(self, point, normal):
        point=np.asarray(point,float)/self.solid.pitch
        normal=np.asarray(normal,float)
        if point.shape!=(3,) or normal.shape!=(3,) or not np.isfinite(point).all() or not np.isfinite(normal).all():
            raise ValueError('FINITE_SECTION_PLANE_REQUIRED')
        norm=np.linalg.norm(normal)
        if not np.isfinite(norm) or norm<=0:
            raise ValueError('NONZERO_SECTION_NORMAL_REQUIRED')
        normal=normal/norm
        scale=max(1.,float(np.max(np.abs(point))),float(np.max(np.abs(self.cells))) if len(self.cells) else 1.)
        tol=256*np.finfo(float).eps*scale
        base=dict(radius=None,area=None,physical_radius_qualified=False,
                  observation_confidence_available=False,unknown_cells_filled=0,
                  definition='area-equivalent radius of one disk-shaped cubical-solid plane component')
        if not len(self.cells):
            return dict(base,status='EMPTY_SOLID')
        active_axis=np.flatnonzero(np.abs(normal)>256*np.finfo(float).eps)
        if len(active_axis)==1 and abs(point[active_axis[0]]-round(point[active_axis[0]]))<=tol:
            return dict(base,status='DEGENERATE_PLANE_COINCIDENT_WITH_GRID_FACE')
        signed=np.einsum('ij,j->i',self.centers-point,normal)
        candidates=np.flatnonzero(np.abs(signed)<=.5*np.sum(np.abs(normal))+tol)
        axis=np.eye(3)[np.argmin(np.abs(normal))]
        u=np.cross(normal,axis);u/=np.linalg.norm(u);v=np.cross(normal,u)
        basis=np.array([u,v]).T
        polygons=[];areas=[];cells=[];contains=[];on_boundary=[]
        edges=defaultdict(list)
        for source in candidates:
            cell=self.cells[source]
            corners=cell+_CORNERS
            distances=np.einsum('ij,j->i',corners-point,normal)
            vertices={}
            for i,j,dimension in _EDGES:
                for at in (i,j):
                    if abs(distances[at])<=tol:
                        key=(0,)+tuple(map(int,corners[at]))+(0,)
                        vertices[key]=corners[at].astype(float)
                if (distances[i]<-tol and distances[j]>tol) or (distances[j]<-tol and distances[i]>tol):
                    t=distances[i]/(distances[i]-distances[j])
                    start=corners[i]
                    key=(1,)+tuple(map(int,start))+(dimension,)
                    vertices[key]=corners[i]+t*(corners[j]-corners[i])
            if len(vertices)<3:
                continue
            keys=list(vertices)
            xy=np.array([(vertices[k]-point)@basis for k in keys])
            center=xy.mean(axis=0);order=np.argsort(np.arctan2(xy[:,1]-center[1],xy[:,0]-center[0]))
            xy=xy[order];keys=[keys[i] for i in order]
            area=.5*float(np.sum(xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1)))
            if area<=tol*tol:
                continue
            delta=np.roll(xy,-1,axis=0)-xy
            cross=delta[:,0]*(-xy[:,1])-delta[:,1]*(-xy[:,0])
            cross_tol=tol*np.maximum(1.,np.linalg.norm(delta,axis=1))
            inside=bool(np.all(cross>=-cross_tol))
            boundary=bool(inside and np.any(np.abs(cross)<=cross_tol))
            index=len(polygons);polygons.append(keys);areas.append(area);cells.append(cell)
            contains.append(inside);on_boundary.append(boundary)
            for a,b in zip(keys,keys[1:]+keys[:1]):
                if a!=b:edges[tuple(sorted((a,b)))].append(index)
        if not polygons:
            return dict(base,status='NO_POSITIVE_AREA_SECTION',candidate_cells=len(candidates))
        hits=np.flatnonzero(contains)
        if not len(hits):
            return dict(base,status='QUERY_OUTSIDE_SECTION',candidate_cells=len(candidates),polygons=len(polygons))
        if len(hits)!=1 or on_boundary[int(hits[0])]:
            return dict(base,status='QUERY_ON_CELL_SECTION_BOUNDARY',candidate_cells=len(candidates),polygons=len(polygons))
        parent=list(range(len(polygons)));contact_edges=[];nonmanifold=[]
        def find(i):
            while parent[i]!=i:
                parent[i]=parent[parent[i]];i=parent[i]
            return i
        def join(a,b):
            a,b=find(a),find(b)
            if a!=b:parent[a]=b
        for edge,owners in edges.items():
            if len(owners)==2:
                a,b=owners;delta=cells[b]-cells[a];manhattan=int(np.sum(np.abs(delta)))
                face_path=manhattan==1
                if manhattan==2 and np.max(np.abs(delta))==1:
                    # An oblique plane can meet two diagonal cells while only
                    # touching intervening occupied cells in zero area. Their
                    # occupied face path resolves the artificial cut-cell split.
                    face_path=any(tuple(cells[a]+np.eye(3,dtype=int)[axis]*delta[axis]) in self.occupied
                                  for axis in np.flatnonzero(delta))
                if face_path:join(a,b)
                else:contact_edges.append(owners)
            elif len(owners)>2:
                nonmanifold.append(owners)
        selected=find(int(hits[0]));ids={i for i in range(len(polygons)) if find(i)==selected}
        area=math.fsum(areas[i] for i in ids)*self.solid.pitch**2
        boundary_graph=defaultdict(set)
        bad=any(any(i in ids for i in group) for group in nonmanifold+contact_edges)
        for (a,b),owners in edges.items():
            n=sum(i in ids for i in owners)
            if n==1:
                boundary_graph[a].add(b);boundary_graph[b].add(a)
            elif n>2:bad=True
        components=len({find(i) for i in range(len(polygons))})
        result=dict(base,area=area,candidate_cells=len(candidates),polygons=len(polygons),
                    selected_polygons=len(ids),plane_components=components,
                    contact_ambiguity_present=bad,exact_cut_polygons=True,
                    numerical_plane_tolerance_lattice_units=tol)
        if bad or not boundary_graph or any(len(ns)!=2 for ns in boundary_graph.values()):
            return dict(result,status='NONMANIFOLD_OR_CONTACT_SECTION')
        visited=set();loops=0
        for start in boundary_graph:
            if start in visited:continue
            loops+=1;todo=[start];visited.add(start)
            while todo:
                at=todo.pop()
                for nxt in boundary_graph[at]-visited:
                    visited.add(nxt);todo.append(nxt)
        result['boundary_loops']=loops
        if loops!=1:
            return dict(result,status='NONDISK_SECTION_HOLES_NOT_FILLED')
        return dict(result,status='CLOSED_DISK_CUBICAL_SECTION',radius=math.sqrt(area/math.pi))
