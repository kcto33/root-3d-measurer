"""Whole cube-union section area/moments; never a disk or radius admission."""
import math
import numpy as np
from root_cubical_sections import _CORNERS,_EDGES
from root_reconstruction_measurement import neighbor_graph


def cube_union_plane_integrals(solid,point,normal):
    p=np.asarray(point,float)/solid.pitch;n=np.asarray(normal,float)
    if p.shape!=(3,) or n.shape!=(3,) or not np.isfinite(p).all() or not np.isfinite(n).all():
        raise ValueError('FINITE_SECTION_PLANE')
    norm=math.hypot(*n)
    if not norm>0:raise ValueError('NONZERO_SECTION_NORMAL')
    n=n/norm;cells=solid.cells;h=solid.pitch
    base=dict(area=None,centroid=None,disk_or_radius_qualified=False,source_material_modified=False)
    if not len(cells):return dict(base,status='EMPTY_SOURCE_SOLID')
    tol=256*np.finfo(float).eps*max(1.,float(np.abs(p).max()),float(np.abs(cells).max()))
    active=np.flatnonzero(np.abs(n)>256*np.finfo(float).eps)
    if len(active)==1:
        axis=int(active[0]);near=abs(p[axis]-round(p[axis]))<=tol
        if near and np.any(n[np.arange(3)!=axis]!=0):
            return dict(base,status='MISSING_NUMERICALLY_NEAR_AXIS_FACE_SECTION_AREA')
        if near:
            level=round(p[axis]);selected=cells[(cells[:,axis]==level)|(cells[:,axis]==level-1)]
            other=[a for a in range(3) if a!=axis]
            xy=np.unique(selected[:,other],axis=0)
            if not len(xy):return dict(base,status='NO_POSITIVE_SECTION_AREA')
            witnesses=np.zeros((len(xy),3));witnesses[:,axis]=p[axis];witnesses[:,other]=xy+.5
            _,labels=neighbor_graph(selected)
            center=np.mean(witnesses,axis=0)*h
            return dict(base,status='COMPLETE_CLOSED_GRID_FACE_UNION_AREA',area=len(xy)*h*h,centroid=center.tolist(),
                polygon_witnesses=(witnesses*h).tolist(),positive_cut_cells=selected.tolist(),
                plane_cell_components=len(set(map(int,labels))),cut_polygons=len(xy),
                axis_face_squares_counted_once=True,numerical_tolerance_lattice=tol)
    basis_axis=np.eye(3)[int(np.argmin(abs(n)))];u=np.cross(n,basis_axis);u/=math.hypot(*u)
    basis=np.array([u,np.cross(n,u)]).T
    signed=(cells+.5-p)@n;candidate=np.flatnonzero(abs(signed)<=.5*np.sum(abs(n))+tol)
    areas=[];moments=[];selected=[];witnesses=[]
    for index in candidate:
        cell=cells[index];corners=cell+_CORNERS;distance=(corners-p)@n;vertices={}
        for i,j,dimension in _EDGES:
            for at in [i,j]:
                if abs(distance[at])<=tol:vertices[(0,)+tuple(map(int,corners[at]))]=corners[at].astype(float)
            if (distance[i]<-tol and distance[j]>tol) or (distance[j]<-tol and distance[i]>tol):
                t=distance[i]/(distance[i]-distance[j])
                vertices[(1,)+tuple(map(int,corners[i]))+(dimension,)]=corners[i]+t*(corners[j]-corners[i])
        if len(vertices)<3:continue
        xy=np.array([(v-p)@basis for v in vertices.values()]);center=xy.mean(0)
        xy=xy[np.argsort(np.arctan2(xy[:,1]-center[1],xy[:,0]-center[0]))]
        nxt=np.roll(xy,-1,axis=0);cross=xy[:,0]*nxt[:,1]-nxt[:,0]*xy[:,1]
        area=math.fsum(map(float,cross))/2
        if area<=tol*tol:continue
        moment=np.array([math.fsum(map(float,(xy[:,a]+nxt[:,a])*cross))/6 for a in range(2)])
        areas.append(area);moments.append(moment);selected.append(cell)
        witnesses.append((p+basis@xy.mean(0))*h)
    if not areas:return dict(base,status='NO_POSITIVE_SECTION_AREA')
    area=math.fsum(areas);center=np.array([math.fsum(float(m[a]) for m in moments)/area for a in range(2)])
    _,labels=neighbor_graph(np.asarray(selected,np.int64))
    return dict(base,status='COMPLETE_WHOLE_CUBE_UNION_SECTION_AREA',area=area*h*h,
        centroid=((p+basis@center)*h).tolist(),polygon_witnesses=np.asarray(witnesses).tolist(),
        positive_cut_cells=np.asarray(selected).tolist(),plane_cell_components=len(set(map(int,labels))),
        cut_polygons=len(areas),axis_face_squares_counted_once=False,numerical_tolerance_lattice=tol,
        components_and_holes_not_filled_or_removed=True)
