"""Conditional geometric round-cap centers from entity section-area profiles."""
import math
import numpy as np
from root_graph_measurement_ledger import graph_measurement_ledger
from root_cubical_ray_measurement import cube_union_exit
from root_section_centroid import centroid_section_class,_normal,_covered


def fit_round_cap_profile(positions,areas,pitch):
    s=np.asarray(positions,float);a=np.asarray(areas,float)
    if s.ndim!=1 or a.shape!=s.shape or len(s)<3 or not np.isfinite(s).all() or not np.isfinite(a).all() or np.any(a<=0):
        raise ValueError('FINITE_POSITIVE_ROUND_CAP_PROFILE_REQUIRED')
    if not math.isfinite(pitch) or pitch<=0:raise ValueError('POSITIVE_PROFILE_PITCH_REQUIRED')
    if np.any(np.diff(s)<=0):raise ValueError('INCREASING_PROFILE_POSITIONS_REQUIRED')
    origin=float(s[0]);x=(s-origin)/pitch;y=a/(math.pi*pitch*pitch)+x*x
    xm=math.fsum(map(float,x))/len(x);ym=math.fsum(map(float,y))/len(y)
    denominator=math.fsum(float(v-xm)**2 for v in x)
    if not denominator>0:raise ValueError('UNRESOLVED_CAP_PROFILE_SPAN')
    slope=math.fsum(float(u-xm)*float(v-ym) for u,v in zip(x,y))/denominator
    c=slope/2;b=ym-slope*xm;radius_squared=b+c*c
    if not math.isfinite(radius_squared) or radius_squared<=0:raise ValueError('NO_POSITIVE_ROUND_CAP_RADIUS')
    center=origin+c*pitch;radius=math.sqrt(radius_squared)*pitch
    predicted=math.pi*(radius*radius-(s-center)**2)
    residual=a-predicted
    # A free quadratic is diagnostic only. It never selects or tunes the model.
    scaled=(s-float(s.mean()))/(s[-1]-s[0])
    coefficients=np.linalg.lstsq(np.c_[scaled*scaled,scaled,np.ones(len(s))],a,rcond=None)[0]
    curvature=float(coefficients[0]/(s[-1]-s[0])**2)
    return dict(status='COMPLETE_CONDITIONAL_ROUND_CAP_PROFILE_FIT',center_distance=float(center),radius=float(radius),
        area_residuals=residual.tolist(),maximum_absolute_area_residual=float(np.max(np.abs(residual))),
        free_area_quadratic_coefficient=curvature,assumed_area_quadratic_coefficient=-math.pi,
        fit_samples=len(s),model='A(s)/pi=R^2-(s-c)^2',
        round_cap_geometry_assumed_not_inferred=True,statistical_confidence_interval=None,
        physical_tip_qualified=False,GT_center_or_radius_used=False)


def round_cap_graph_measurement(solid,points,adjacency,*,round_caps_assumed,body_points=None):
    p=np.asarray(points,float);source=graph_measurement_ledger(p,adjacency,coordinate_units='caller_solid_units')
    q=p.copy() if body_points is None else np.asarray(body_points,float).copy()
    if q.shape!=p.shape or not np.isfinite(q).all():raise ValueError('BODY_GRAPH_EXACT_ID_SHAPE_REQUIRED')
    engine=centroid_section_class()(solid);h=solid.pitch;rows=[];occupied=set(map(tuple,solid.cells))
    leaf_paths={}
    for branch in source['branches']:
        ids=branch['source_node_ids']
        if len(adjacency[ids[0]])==1:leaf_paths[ids[0]]=ids
        if len(adjacency[ids[-1]])==1:leaf_paths[ids[-1]]=list(reversed(ids))
    for leaf in source['endpoint_node_ids']:
        row=dict(node_id=leaf,original_position=p[leaf].tolist(),estimated_geometric_cap_center=None,
            physical_root_tip=None,physical_tip_extension=0.,round_cap_model_assumed=bool(round_caps_assumed))
        if not round_caps_assumed:
            rows.append(dict(row,status='MISSING_DECLARED_ROUND_CAP_MODEL'));continue
        path=leaf_paths[leaf];inward=_normal(p[path],0,h,False)
        if inward is None:
            rows.append(dict(row,status='MISSING_TERMINAL_GEOMETRIC_NORMAL'));continue
        direction=-inward
        seed=engine.section(p[leaf],direction)
        if seed['radius'] is None:
            rows.append(dict(row,status='MISSING_SEED_DISK_SECTION',seed_section_status=seed['status']));continue
        try:exit=cube_union_exit(solid.cells,h,p[leaf],direction)
        except ValueError as error:
            rows.append(dict(row,status='MISSING_GEOMETRIC_PROFILE_BOUNDARY',reason=str(error)));continue
        end=exit['boundary_distance']-.25*h;begin=exit['boundary_distance']-seed['radius']
        if not end>begin:
            rows.append(dict(row,status='MISSING_RESOLVED_CAP_PROFILE_WINDOW'));continue
        distances=np.linspace(begin,end,21);sections=[engine.section(p[leaf]+s*direction,direction) for s in distances]
        row.update(profile_positions=distances.tolist(),section_statuses=[r['status'] for r in sections],
            geometric_exit_only_sampling_bound=exit,seed_radius=seed['radius'])
        if any(r['radius'] is None for r in sections):
            rows.append(dict(row,status='MISSING_COMPLETE_DISK_CAP_PROFILE'));continue
        areas=[r['area'] for r in sections]
        try:fit=fit_round_cap_profile(distances,areas,h)
        except ValueError as error:
            rows.append(dict(row,status='MISSING_POSITIVE_ROUND_CAP_FIT',reason=str(error),profile_areas=areas));continue
        transverse=np.mean([np.array(r['section_centroid_evidence']['centroid'])-(p[leaf]+s*direction)
                            for r,s in zip(sections,distances)],axis=0)
        transverse-=np.dot(transverse,direction)*direction
        center=p[leaf]+fit['center_distance']*direction+transverse
        row.update(fit=fit,profile_areas=areas,conditional_fitted_center=center.tolist())
        if tuple(np.floor(center/h).astype(np.int64)) not in occupied:
            rows.append(dict(row,status='MISSING_FITTED_CENTER_INSIDE_SOURCE'));continue
        neighbor=next(iter(adjacency[leaf]))
        if not _covered(solid.cells,h,center,q[neighbor]):
            rows.append(dict(row,status='MISSING_FITTED_TERMINAL_EDGE_SUPPORT'));continue
        q[leaf]=center
        rows.append(dict(row,status='COMPLETE_CONDITIONAL_GEOMETRIC_CAP_CENTER',estimated_geometric_cap_center=center.tolist(),
            geometric_displacement=math.dist(center,p[leaf]),geometric_cap_radius=fit['radius']))
    result=graph_measurement_ledger(q,adjacency,coordinate_units='caller_solid_units')
    assert result['source_edge_ids']==source['source_edge_ids']
    return dict(status='COMPLETE_CONDITIONAL_ROUND_CAP_GRAPH_MEASUREMENT',source_graph_ledger=source,
        conditional_graph_ledger=result,points=q.tolist(),adjacency=[sorted(map(int,ns)) for ns in adjacency],endpoints=rows,
        available_cap_centers=sum(r['estimated_geometric_cap_center'] is not None for r in rows),
        original_endpoint_positions_retained_in_receipts=True,all_source_edges_and_IDs_preserved=True,
        source_material_modified=False,GT_reference_accepted_by_inference=False,
        ray_exit_adopted_as_endpoint=False,round_cap_model_assumed_not_verified=bool(round_caps_assumed),
        no_residual_based_model_selection=True,production_metric_adopted=False,
        physical_root_length=None,physical_root_tip=None,physical_tip_extension=0.)
