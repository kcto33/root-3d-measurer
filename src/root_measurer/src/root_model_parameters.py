"""One parameter contract for exact input geometry and measured axis samples.

Missing intervals are never filled. Radius statistics describe valid measured
output arc only; a partial centerline is never relabelled complete model length.
"""
import csv,io,json,math
from pathlib import Path
import numpy as np


def radius_statistics(samples):
    total=math.fsum(s['arc_weight'] for s in samples)
    valid=[s for s in samples if s.get('available') and s.get('radius') is not None]
    if any(not math.isfinite(s['radius']) or s['radius']<0 or s['arc_weight']<0 for s in valid):
        raise ValueError('FINITE_NONNEGATIVE_RADIUS_AND_ARC_REQUIRED')
    weight=math.fsum(s['arc_weight'] for s in valid)
    values=sorted((s['radius'],s['arc_weight']) for s in valid if s['arc_weight']>0)
    def quantile(q):
        cumulative=0.
        for value,w in values:
            cumulative+=w
            if cumulative>=q*weight:return value
        return values[-1][0]
    return dict(scope='valid measured output arc only',full_output_arc=total,valid_output_arc=weight,
        valid_arc_fraction=weight/total if total else None,missing_output_arc=max(0.,total-weight),
        mean=math.fsum(v*w for v,w in values)/weight if weight else None,
        minimum=values[0][0] if values else None,maximum=values[-1][0] if values else None,
        median=quantile(.5) if weight else None,p95=quantile(.95) if weight else None,
        weighting='edge sample arc weights; quantiles use left-continuous weighted empirical inverse CDF',
        missing_interpolated=False)


def build_parameters(report,*,solid=None,boundary_face_count=None,mesh=None,morphology_options=None):
    from root_model_delivery_exports import measurement_tables
    from root_graph_measurement_ledger import graph_measurement_ledger
    from root_model_delivery_validation import complete_length_state
    paths,samples,coverage=measurement_tables(report)
    graph=report.get('candidate_centerline') or dict(points=[],adjacency=[])
    ledger=graph_measurement_ledger(graph['points'],graph['adjacency'],coordinate_units=report['coordinate_units']) if graph['points'] else None
    state=complete_length_state(report)
    blocked=state['complete_length_status'].startswith('BLOCKED') or bool(state['isolated_component_node_ids']) or 'BLOCKED' in report.get('status','')
    scale=report.get('scale_mm_per_model_unit')
    if scale is not None and (not math.isfinite(scale) or scale<=0):raise ValueError('POSITIVE_FINITE_SCALE_REQUIRED')
    unit='mm' if scale is not None else report['coordinate_units'];factor=1. if scale is None else scale
    radii=radius_statistics(samples)
    connected=radius_statistics([dict(s,radius=s.get('connected_section_radius'),
        available=s.get('connected_section_radius') is not None) for s in samples])
    radius_display={k:(None if radii[k] is None else radii[k]*factor) for k in ['mean','minimum','maximum','median','p95']}
    geometry=dict(kind='NOT_PROVIDED',volume=None,surface_area=None)
    if solid is not None:
        if boundary_face_count is None:
            from root_reconstruction_measurement import boundary_faces
            boundary_face_count=len(boundary_faces(solid)[0])
        geometry=dict(kind='EXACT_NATIVE_CUBE_UNION',volume=len(solid.cells)*solid.pitch**3,
            surface_area=int(boundary_face_count)*solid.pitch**2,occupied_cells=len(solid.cells),
            boundary_faces=int(boundary_face_count),pitch=solid.pitch,
            definition='all original occupied cubes and exposed faces, including cavity boundaries; no tube approximation')
    elif mesh is not None:
        triangles=mesh.vertices[mesh.faces]
        area=float(np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1).sum()/2)
        geometry=dict(kind='ORIGINAL_TRIANGLE_AREA_SUM',volume=None,surface_area=None,
            triangle_area_sum=area,triangles=len(mesh.faces),
            definition='original triangle area sum; contact/overlap means union surface/solid volume not automatically qualified')
    bounds=report.get('model_bounds');dimensions=None if bounds is None else [b-a for a,b in zip(bounds['minimum'],bounds['maximum'])]
    rows=[]
    def metric(key,label,value,dimension,status,definition,coverage_fraction=None):
        rows.append(dict(parameter=key,label=label,value=value,unit='count' if dimension==0 else unit if dimension==1 else unit+'^'+str(dimension),
            status=status,covered_arc_fraction=coverage_fraction,definition=definition))
    length=coverage['output_total_length']
    metric('model_length','完整几何轴长度',None if blocked else length*factor,1,
        'BLOCKED' if blocked else 'GEOMETRIC_ESTIMATE_ACCURACY_SEPARATE','unique final graph edges; all geometric endpoints resolved and no isolated components required')
    metric('extracted_axis_length','已提取唯一轴长',length*factor,1,'PARTIAL' if blocked else 'AVAILABLE','unique final graph edges; not a biological root-to-tip cumulative length')
    for key,label in [('mean','弧长加权平均半径'),('minimum','最小实测半径'),('maximum','最大实测半径'),('median','弧长加权半径中位数'),('p95','弧长加权半径95分位')]:
        metric('radius_'+key,label,radius_display[key],1,'CONDITIONAL_VALID_ARC' if radii['valid_arc_fraction']!=1 else 'AVAILABLE',
            'area-equivalent radius of unchanged original-entity sections; valid measured output arc only',radii['valid_arc_fraction'])
    metric('diameter_mean','单分支平均直径',None if radius_display['mean'] is None else 2*radius_display['mean'],1,'CONDITIONAL_VALID_ARC','twice area-equivalent mean radius',radii['valid_arc_fraction'])
    metric('connected_section_radius_mean','整连通截面等效半径均值',None if connected['mean'] is None else connected['mean']*factor,1,
        'ORIGINAL_CONNECTED_SECTION_GEOMETRY','sqrt(entire connected section area/pi), holes deducted; includes labelled mixed/junction sections',connected['valid_arc_fraction'])
    metric('connected_section_diameter_mean','整连通截面等效直径均值',None if connected['mean'] is None else 2*connected['mean']*factor,1,
        'ORIGINAL_CONNECTED_SECTION_GEOMETRY','twice entire connected-section equivalent radius',connected['valid_arc_fraction'])
    for key,label,power in [('volume','实体体积',3),('surface_area','实体边界面积',2)]:
        value=geometry.get(key)
        metric(key,label,None if value is None else value*factor**power,power,'EXACT_INPUT_GEOMETRY' if value is not None else 'UNQUALIFIED',geometry.get('definition','input geometry unavailable'))
    for i,label in enumerate(['包围盒X尺寸','包围盒Y尺寸','包围盒Z尺寸']):
        metric('bbox_'+'xyz'[i],label,None if dimensions is None else dimensions[i]*factor,1,'AVAILABLE','axis-aligned original-model bounding dimensions')
    counts=dict(paths=len(paths),components=0 if ledger is None else ledger['source_components'],
        endpoints=state['geometric_endpoint_count'],resolved_endpoints=state['estimated_geometric_endpoint_count'],
        junctions=sum(len(a)>=3 for a in graph['adjacency']),cycles=0 if ledger is None else ledger['source_cycle_rank'])
    for key,label in [('paths','几何路径数'),('components','中心线分量数'),('endpoints','图端点数'),('resolved_endpoints','已估计几何端点数'),('junctions','分叉节点数'),('cycles','实际保留环数')]:
        metric(key,label,counts[key],0,'GRAPH_DESCRIPTOR','finite measurement graph; no biological primary/lateral identity')
    per_path=[]
    for path in paths:
        stats=radius_statistics(path['samples'])
        per_path.append(dict(path_id=path['path_id'],component_id=path['component_id'],closed=path['closed'],
            length=path['length']*factor,start_status=path['start_status'],end_status=path['end_status'],
            radius_mean=None if stats['mean'] is None else stats['mean']*factor,
            radius_min=None if stats['minimum'] is None else stats['minimum']*factor,
            radius_max=None if stats['maximum'] is None else stats['maximum']*factor,
            diameter_mean=None if stats['mean'] is None else 2*stats['mean']*factor,
            valid_radius_arc_fraction=stats['valid_arc_fraction'],radius_scope=stats['scope']))
    from root_model_morphology import build_morphology
    morphology=build_morphology(report,solid=solid,mesh=mesh,options=morphology_options,geometry=geometry)
    rows.extend(morphology['metrics'])
    return dict(schema_version='root-model-parameters-v001',model=report.get('case_label'),unit=unit,morphology=morphology,
        physical_scale_supplied=scale is not None,scale_mm_per_model_unit=scale,
        geometry=geometry,radius_model_units=radii,connected_section_radius_model_units=connected,coverage=coverage,length_state=state,
        complete_model_length_available=not blocked,metrics=rows,paths=per_path,
        reference_precision_separate=True,source_geometry_unchanged=True,
        measured_radius_scope='valid output arc; not complete-reference radius coverage')


def parameter_csv(parameters):
    output=io.StringIO(newline='');fields=['parameter','label','value','unit','status','covered_arc_fraction','definition']
    writer=csv.DictWriter(output,fieldnames=fields,lineterminator='\r\n');writer.writeheader();writer.writerows(parameters['metrics'])
    return output.getvalue().encode('utf-8-sig')


def write_parameters(report,prefix,**geometry):
    from root_model_delivery_exports import file_reference,write_json
    value=build_parameters(report,**geometry);prefix=Path(prefix)
    path=prefix.with_name(prefix.name+'_parameters.json');write_json(path,value)
    csv_path=prefix.with_name(prefix.name+'_parameters.csv')
    with csv_path.open('xb') as f:f.write(parameter_csv(value))
    return dict(json=file_reference(path,prefix.parent),csv=file_reference(csv_path,prefix.parent))
