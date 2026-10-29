"""Conditional flat-cap plane location from terminal section-volume conservation."""
import math
import numpy as np
from root_graph_measurement_ledger import graph_measurement_ledger
from root_cube_plane_integrals import cube_union_plane_integrals
from root_branch_section_measurement import BranchOwnership
from root_cubical_ray_measurement import cube_union_exit
from root_section_centroid import _normal,_covered
from root_cube_slab_volume import cube_union_slab_volume


def fit_flat_cap_profile(interior_positions,interior_areas,begin,end,section_volume,pitch):
    """Fit sqrt(A)=a+b*s inside; conserve the supplied terminal slab volume."""
    s=np.asarray(interior_positions,float);a=np.asarray(interior_areas,float)
    if s.ndim!=1 or a.shape!=s.shape or len(s)<3 or not np.isfinite(s).all() or not np.isfinite(a).all() or np.any(a<=0):
        raise ValueError('POSITIVE_INTERIOR_AREA_PROFILE_REQUIRED')
    if not math.isfinite(pitch) or pitch<=0 or not np.all(np.diff(s)>0):raise ValueError('RESOLVED_PROFILE_SCALE_REQUIRED')
    if not all(math.isfinite(x) for x in (begin,end,section_volume)) or not begin<end or section_volume<=0:
        raise ValueError('POSITIVE_FINITE_TERMINAL_SLAB_REQUIRED')
    x=(s-begin)/pitch;y=np.sqrt(a)/pitch
    xm=float(np.mean(x));ym=float(np.mean(y));den=float(np.dot(x-xm,x-xm))
    if den<=0:raise ValueError('RESOLVED_INTERIOR_SPAN_REQUIRED')
    slope=float(np.dot(x-xm,y-ym)/den);intercept=ym-slope*xm
    limit=(end-begin)/pitch
    if min(intercept,intercept+slope*limit)<=0:raise ValueError('POSITIVE_EXTRAPOLATED_AREA_MODEL_REQUIRED')
    def integral(t):return intercept**2*t+intercept*slope*t*t+slope*slope*t**3/3
    volume=section_volume/pitch**3
    if volume>integral(limit):raise ValueError('TERMINAL_VOLUME_EXCEEDS_REGISTERED_MODEL_WINDOW')
    lo=0.;hi=limit
    for _ in range(80):
        mid=(lo+hi)/2
        if integral(mid)<volume:lo=mid
        else:hi=mid
    distance=begin+(lo+hi)*pitch/2
    prediction=(intercept+slope*x)**2*pitch*pitch
    return dict(status='COMPLETE_CONDITIONAL_FLAT_CAP_VOLUME_FIT',cap_distance=float(distance),
        sqrt_area_intercept_at_begin=intercept*pitch,sqrt_area_slope=slope,
        model='A(s)=square(linear sqrt(area)) before flat cap; zero beyond',
        interior_area_residuals=(a-prediction).tolist(),terminal_section_volume=section_volume,
        flat_cap_and_linear_radial_taper_assumed_not_inferred=True,
        GT_endpoint_or_radius_used=False,statistical_confidence_interval=None,physical_tip_qualified=False)


def flat_cap_graph_measurement(solid,points,adjacency,*,flat_caps_assumed,body_points=None):
    p=np.asarray(points,float);q=p.copy() if body_points is None else np.asarray(body_points,float).copy()
    if q.shape!=p.shape or not np.isfinite(q).all():raise ValueError('BODY_GRAPH_EXACT_ID_SHAPE_REQUIRED')
    source=graph_measurement_ledger(p,adjacency,coordinate_units='caller_solid_units')
    ownership=BranchOwnership(p,adjacency);h=solid.pitch;paths={};rows=[]
    for branch in source['branches']:
        ids=branch['source_node_ids']
        if len(adjacency[ids[0]])==1:paths[ids[0]]=ids
        if len(adjacency[ids[-1]])==1:paths[ids[-1]]=list(reversed(ids))
    for leaf in source['endpoint_node_ids']:
        row=dict(node_id=leaf,original_position=p[leaf].tolist(),estimated_geometric_cap_center=None,
            physical_root_tip=None,flat_cap_model_assumed=bool(flat_caps_assumed))
        if not flat_caps_assumed:rows.append(dict(row,status='MISSING_DECLARED_FLAT_CAP_MODEL'));continue
        inward=_normal(p[paths[leaf]],0,h,False)
        if inward is None:rows.append(dict(row,status='MISSING_TERMINAL_GEOMETRIC_NORMAL'));continue
        direction=-inward;owner=ownership.node_branches[leaf]
        if len(owner)!=1:rows.append(dict(row,status='MISSING_UNIQUE_TERMINAL_BRANCH'));continue
        def section(s,allow_empty=False):
            result=cube_union_plane_integrals(solid,p[leaf]+s*direction,direction)
            if result['area'] is None or result['area']<=0:
                if allow_empty and result['status'] in ('NO_POSITIVE_SECTION_AREA','EMPTY_SOURCE_SOLID'):
                    return dict(area=0.,centroid=None,status='EMPTY_CAP_PROFILE_SECTION')
                raise ValueError('MISSING_POSITIVE_CAP_PROFILE:'+result['status'])
            if result['plane_cell_components']!=1:raise ValueError('MISSING_SINGLE_CAP_PLANE_CELL_GROUP')
            if any(ownership.nearest_owners(w)!=owner for w in result['polygon_witnesses']):
                raise ValueError('MISSING_SAME_BRANCH_CAP_PROFILE')
            return result
        try:
            seed=section(0.);radius=math.sqrt(seed['area']/math.pi)
            begin=-radius;exit=cube_union_exit(solid.cells,h,p[leaf],direction);end=exit['boundary_distance']+2*h
            interior_s=np.linspace(begin,0.,21);interior=[section(s) for s in interior_s]
            positions=np.linspace(begin,end,129);profiles=[section(s,True) for s in positions]
            values=[r['area'] for r in profiles]
            slab=cube_union_slab_volume(solid,direction,p[leaf]+begin*direction,p[leaf]+end*direction)
            if any(ownership.nearest_owners((solid.cells[i]+.5)*h)!=owner for i in slab['intersecting_source_cell_indices']):
                raise ValueError('MISSING_SAME_BRANCH_WHOLE_TERMINAL_SLAB')
            volume=slab['volume']/math.hypot(*direction)
            fit=fit_flat_cap_profile(interior_s,[r['area'] for r in interior],begin,end,volume,h)
            offset=np.mean([np.asarray(r['centroid'])-(p[leaf]+s*direction) for s,r in zip(interior_s,interior)],axis=0)
            offset-=np.dot(offset,direction)*direction
            center=p[leaf]+fit['cap_distance']*direction+offset
            row.update(fit=fit,geometric_exit_only_sampling_bound=exit,interior_positions=interior_s.tolist(),
                interior_areas=[r['area'] for r in interior],integration_positions=positions.tolist(),
                integration_areas=values,exact_terminal_slab=slab,conditional_fitted_center=center.tolist())
            neighbor=next(iter(adjacency[leaf]))
            if not _covered(solid.cells,h,center,q[neighbor]):raise ValueError('MISSING_FITTED_FLAT_CAP_EDGE_SUPPORT')
            q[leaf]=center
            rows.append(dict(row,status='COMPLETE_CONDITIONAL_GEOMETRIC_FLAT_CAP_CENTER',estimated_geometric_cap_center=center.tolist()))
        except ValueError as error:rows.append(dict(row,status=str(error)))
    candidate=graph_measurement_ledger(q,adjacency,coordinate_units='caller_solid_units')
    assert candidate['source_edge_ids']==source['source_edge_ids']
    return dict(status='COMPLETE_CONDITIONAL_FLAT_CAP_GRAPH_MEASUREMENT',source_graph_ledger=source,
        conditional_graph_ledger=candidate,points=q.tolist(),adjacency=[sorted(map(int,a)) for a in adjacency],
        endpoints=rows,available_cap_centers=sum(r['estimated_geometric_cap_center'] is not None for r in rows),
        source_material_modified=False,original_endpoint_positions_retained=True,all_source_IDs_and_edges_preserved=True,
        GT_reference_accepted_by_inference=False,ray_exit_adopted_as_endpoint=False,flat_cap_model_assumed_not_verified=bool(flat_caps_assumed),
        fixed_129_samples_for_diagnostics_and_gates_only=True,exact_source_cube_slab_volume=True,production_metric_adopted=False,
        physical_root_length=None,physical_root_tip=None,physical_tip_extension=0.)
