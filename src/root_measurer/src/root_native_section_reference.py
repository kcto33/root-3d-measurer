"""Independent cube/plane polygon union reference, using Shapely/GEOS.

No predictor clipping, disk, ownership or moment routines are imported.
All occupied cubes intersecting the infinite query plane are considered;
unrelated nonmanifold regions do not invalidate the queried disk.
See https://shapely.readthedocs.io/en/2.1.2/reference/shapely.union_all.html
"""
import itertools
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

from shapely import Polygon,Point,union_all,__version__ as SHAPELY_VERSION

CORNERS=np.array(list(itertools.product([0,1],repeat=3)),np.int64)
EDGES=[(i,j) for i,a in enumerate(CORNERS) for j,b in enumerate(CORNERS) if i<j and np.sum(a!=b)==1]
NEIGHBORS=[d for d in itertools.product([-1,0,1],repeat=3) if d!=(0,0,0)]


class NativePlaneReference:
    def __init__(self,solid,*,whole_plane=False):
        self.solid=solid
        self.whole_plane=whole_plane
        self.centers=(solid.cells+.5)*solid.pitch

    def section(self,point,normal,*,include_witnesses=False,connected_area=False):
        point=np.asarray(point,float);normal=np.asarray(normal,float)
        if point.shape!=(3,) or normal.shape!=(3,) or not np.isfinite(np.r_[point,normal]).all():
            raise ValueError('FINITE_REFERENCE_PLANE_REQUIRED')
        size=math.hypot(*normal)
        if size==0:raise ValueError('NONZERO_REFERENCE_NORMAL_REQUIRED')
        normal=normal/size;h=self.solid.pitch
        aligned=np.flatnonzero(normal)
        exact_axis=int(aligned[0]) if len(aligned)==1 else None
        axis=np.eye(3)[int(np.argmin(np.abs(normal)))];u=np.cross(normal,axis);u/=math.hypot(*u);v=np.cross(normal,u)
        distances=sum((self.centers[:,a]-point[a])*normal[a] for a in range(3))
        selected=self.solid.cells[np.abs(distances)<=h*np.abs(normal).sum()/2+32*np.finfo(float).eps*max(1.,math.hypot(*point))]
        global_candidate_count=len(selected)
        candidates=set(map(tuple,selected));seed=np.floor(point/h).astype(int)
        component={tuple(seed+d) for d in itertools.product([-1,0,1],repeat=3) if tuple(seed+d) in candidates}
        queue=list(component)
        while queue:
            cell=queue.pop()
            for delta in NEIGHBORS:
                neighbor=tuple(cell[a]+delta[a] for a in range(3))
                if neighbor in candidates and neighbor not in component:component.add(neighbor);queue.append(neighbor)
        # Every continuously connected cut path traverses touching cubes;
        # those cube indices are26-neighbors. This is a complete superset,
        # not a radius-limited ROI or a predictor section implementation.
        if not self.whole_plane:selected=np.array(sorted(component),np.int64).reshape(-1,3)
        polygons=[];edge_cache={};vertex_cache={}
        for cell in selected:
            grid=cell+CORNERS;vertices=grid*h;values=(vertices-point)@normal;cuts={}
            if exact_axis is not None:
                # The lattice is integer*pitch, not its rounded world vertex.
                # A binary64 query infinitesimally off a face must not combine
                # both adjacent layers as a falsely coplanar section.
                values=np.array([float(Fraction(int(index))*Fraction(h)-Fraction(float(point[exact_axis])))
                                 *normal[exact_axis] for index in grid[:,exact_axis]])
            for i,j in EDGES:
                for k in [i,j]:
                    if values[k]==0:
                        key=tuple(map(int,grid[k]))
                        if key not in vertex_cache:vertex_cache[key]=tuple((vertices[k]-point)@basis for basis in [u,v])
                        cuts[('v',key)]=vertex_cache[key]
                if values[i]*values[j]<0:
                    # Same global directed lattice edge always has the same
                    # interpolation and projected coordinates in both cubes.
                    lower=np.minimum(grid[i],grid[j]);direction=int(np.flatnonzero(grid[i]!=grid[j])[0]);key=(*map(int,lower),direction)
                    if key not in edge_cache:
                        a=lower*h;b=a.copy();b[direction]+=h
                        da=float((a-point)@normal);db=float((b-point)@normal)
                        q=a+da/(da-db)*(b-a)
                        edge_cache[key]=tuple(float((q-point)@basis) for basis in [u,v])
                    cuts[('e',key)]=edge_cache[key]
            xy=np.unique(np.array(list(cuts.values()),float).reshape(-1,2),axis=0)
            if len(xy)<3:continue
            center=xy.mean(axis=0);order=np.argsort(np.arctan2(xy[:,1]-center[1],xy[:,0]-center[0]))
            polygon=Polygon(xy[order])
            if polygon.area>0:
                if not polygon.is_valid:return dict(radius=None,status='REFERENCE_CELL_CLIP_INVALID')
                polygons.append(polygon)
        if not polygons:return dict(radius=None,status='REFERENCE_NO_POSITIVE_AREA')
        try:region=union_all(polygons)
        except Exception as exc:return dict(radius=None,status='REFERENCE_UNION_UNQUALIFIED',reason=str(exc))
        if not region.is_valid:return dict(radius=None,status='REFERENCE_UNION_INVALID')
        if region.geom_type=='Polygon':regions=[region]
        elif region.geom_type=='MultiPolygon':regions=list(region.geoms)
        else:return dict(radius=None,status='REFERENCE_NONAREA_COMPLEX')
        containing=[g for g in regions if g.contains(Point(0,0))]
        if len(containing)!=1:return dict(radius=None,status='REFERENCE_QUERY_NOT_STRICTLY_INSIDE_ONE_COMPONENT')
        chosen=containing[0]
        if len(chosen.interiors) and not connected_area:return dict(radius=None,status='REFERENCE_NONDISK_HOLE')
        if any(chosen.intersects(g) for g in regions if g is not chosen):
            return dict(radius=None,status='REFERENCE_POINT_CONTACT')
        result=dict(radius=math.sqrt(chosen.area/math.pi),area=float(chosen.area),
                    status='REFERENCE_COMPLETE_CUBE_PLANE_GEOS_DISK',
                    independent_of_predictor_disk_kernel=True,shapely_version=SHAPELY_VERSION,
                    plane_moved=False,source_cells_considered=len(selected),
                    whole_plane_hits=global_candidate_count,whole_plane_query_component_complete=True,
                    finite_roi_clipped=False,connected_area=connected_area,
                    holes=len(chosen.interiors),holes_filled=False)
        if include_witnesses:
            # Each convex cube cut has an interior centroid. Select only cuts
            # in the queried component; disconnected branches are irrelevant.
            result['interior_cube_cut_witnesses']=[
                (point+g.centroid.x*u+g.centroid.y*v).tolist()
                for g in polygons if chosen.contains(g.centroid)]
        return result

    def owned_section(self,point,normal,ownership,reference_path):
        result=self.section(point,normal,include_witnesses=True)
        witnesses=result.pop('interior_cube_cut_witnesses',[])
        if result.get('radius') is None:return result
        foreign=sum(ownership.nearest_owners(q)!={reference_path} for q in witnesses)
        result.update(reference_path=reference_path,reference_branch_witnesses=len(witnesses),
                      foreign_reference_branch_witnesses=foreign,
                      reference_ownership='nearest reference path at every interior cube-cut centroid; finite witnesses, not whole-region proof')
        if foreign or not witnesses:
            result.update(radius=None,status='REFERENCE_MIXED_BRANCH_SECTION')
        return result
