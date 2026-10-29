"""Mesh-only automatic geometric prototype; proxy and Source evidence separate."""
import hashlib,math,os
from fractions import Fraction as F
import numpy as np
from root_triangle_measurement_input import TriangleInput
from root_closed_triangle_voxelization import voxelize_closed_chains
from root_reconstruction_measurement import ReconstructedSolid
from root_tube_measurement import measure_tube_solid
from root_body_measurement_api import measure_centerline_body
from root_axis_radius_profile import axis_radius_profile
from root_graph_measurement_ledger import graph_measurement_ledger
from root_exact_query_sections import ExactQuerySectionMeasurer
from root_graph_source_length import pack

def measure_model(mesh,*,coordinate_units='model_units_unspecified_not_mm',grid_divisions=16,execution_domain='local'):
    if not isinstance(mesh,TriangleInput):raise ValueError('IMMUTABLE_TRIANGLE_MODEL_REQUIRED')
    if not isinstance(coordinate_units,str) or not coordinate_units.strip():raise ValueError('MODEL_UNIT_LABEL_REQUIRED')
    if type(grid_divisions)is not int or not 8<=grid_divisions<=32:raise ValueError('FIXED_GLOBAL_RESOLUTION_GUARD_8_TO_32')
    nv,nf=len(mesh.vertices),len(mesh.faces)
    if not 0<nv<=1200000 or not 0<nf<=2000000:raise ValueError('COMPLETE_SOURCE_GUARDS')
    if execution_domain=='local':
        if nv>3000 or nf>1000:raise ValueError('UNCHANGED_LOCAL_SOURCE_GUARDS')
    elif execution_domain!='registered_slurm_complete_source' or not os.environ.get('SLURM_JOB_ID','').isdigit():raise ValueError('REGISTERED_SLURM_FOR_COMPLETE_MODEL_REQUIRED')
    source_before=hashlib.sha256(mesh.vertices.tobytes()+mesh.faces.tobytes()).hexdigest()
    v=np.asarray(mesh.vertices);f=np.asarray(mesh.faces)
    if v.ndim!=2 or v.shape[1]!=3 or not np.isfinite(v).all() or f.ndim!=2 or f.shape[1]!=3 or f.dtype.kind not in 'iu' or f.min()<0 or f.max()>=nv:raise ValueError('FINITE_ORIGINAL_TRIANGLE_MODEL_REQUIRED')
    bounds=dict(minimum=v.min(axis=0).tolist(),maximum=v.max(axis=0).tolist());extent=float(np.ptp(v,axis=0).max())
    if not math.isfinite(extent) or extent<=0:raise ValueError('POSITIVE_FINITE_MODEL_EXTENT_REQUIRED')
    h=extent/grid_divisions;failures=[];stages={};graph=None;length=None;radius=None;source_checks=[];stage='EXACT_MODEL_CENTRE_SAMPLING'
    def fail(code,*,position=None,location_kind='whole_model_bounds',node_id=None,edge_id=None,details=None):
        failures.append(dict(code=code,position=position,location_kind=location_kind,model_bounds=bounds if position is None else None,generated_node_id=node_id,generated_edge_id=edge_id,details=details))
    try:
        cells,receipt=voxelize_closed_chains(v,f,h,max_xy_tests=10000000,max_output_cells=100000)
        stages['voxelization']=receipt
        if not len(cells):raise ValueError('EMPTY_CENTRE_SAMPLED_MODEL_AT_FIXED_RESOLUTION')
        solid=ReconstructedSolid(cells,h,dict(original_mesh_sha256=mesh.source_sha256,derived_proxy=True,grid_origin=[0,0,0]))
        stage='AUTOMATIC_GEOMETRIC_SKELETON'
        kernel=measure_tube_solid(solid,max_dense_cells=1000000)
        stages['skeleton']=dict(status=kernel['status'],thinning=kernel['thinning'],representation=kernel['representation'],generated_nodes=len(kernel['points']))
        if execution_domain=='local' and len(kernel['points'])>256:raise ValueError('UNCHANGED_LOCAL_256_GENERATED_NODE_GUARD')
        if not kernel['points']:raise ValueError('EMPTY_AUTOMATIC_SKELETON')
        stage='AUTOMATIC_PROXY_BODY_CENTRING'
        body=measure_centerline_body(solid,coordinate_units=coordinate_units,kernel_graph=kernel,include_radius_profile=False)
        stages['proxy_body']=dict(supported=body['body_geometry'] is not None,source_representation=body['source_solid_measurement']['source_representation_receipt'])
        if body['body_geometry'] is not None:
            points=np.asarray(body['body_geometry']['points'],float);adj=body['body_geometry']['adjacency'];mode='AUTOMATIC_VOXEL_SECTION_CENTROID_CANDIDATE'
        else:
            points=np.asarray(kernel['points'],float);adj=kernel['adjacency'];mode='UNSUPPORTED_AUTOMATIC_SKELETON_DIAGNOSTIC'
            fail('VOXEL_COMPLEX_NOT_ADMITTED_AS_1D_CENTRELINE')
        ledger=graph_measurement_ledger(points,adj,coordinate_units=coordinate_units)
        graph=dict(points=points.tolist(),adjacency=adj,automatic_node_ids=True,caller_node_ids_used=False,definition=mode,components=ledger['source_components'],cycle_rank=ledger['source_cycle_rank'],branches=ledger['branches'])
        length=dict(value=ledger['measurement_graph_length'],definition='sum of unique automatically derived candidate graph segments in original model units',coordinate_units=coordinate_units,
            full_model_terminal_coverage=None,physical_length=None,resolution_error_bound=None)
        for i,ns in enumerate(adj):
            if len(ns)!=2:fail('ENDPOINT_OR_JUNCTION_NO_UNIQUE_TWO_SIDED_NORMAL',position=points[i].tolist(),location_kind='generated_node',node_id=i)
        stage='PROXY_AREA_RADIUS_CURVE'
        proxy=axis_radius_profile(solid,points,adj)
        radius=dict(provider='DERIVED_VOXEL_PROXY_AREA_EQUIVALENT_SECTION',definition=proxy['radius_definition'],pitch=h,rows=proxy['rows'],radii=proxy['radii'],coverage=proxy['new_coverage'],
            original_triangle_radius_qualified=False,proxy_to_original_geometry_error_bound=None,interpolated_values=None)
        for row in proxy['rows']:
            if row['new_radius'] is None:fail(row['new_status'],position=row['position'],location_kind='generated_node',node_id=row['node_id'])
        # Freeze all five station choices before observing any Source outcome.
        # This is an explicitly sparse original-mesh check, never a complete field.
        eligible=[row for row in proxy['rows'] if row['normal'] is not None]
        selected=[]
        if eligible:
            selected=sorted({min(len(eligible)-1,(2*i+1)*len(eligible)//10) for i in range(min(5,len(eligible)))})
        stations=[eligible[i] for i in selected]
        stages['fixed_automatic_source_station_ids']=[r['node_id'] for r in stations]
        stage='FRESH_ORIGINAL_MESH_SPARSE_SECTIONS'
        engine=ExactQuerySectionMeasurer(mesh,execution_domain=execution_domain,max_queries=5)
        for row in stations:
            point=[pack(F(float(x))) for x in row['position']];normal=[pack(F(float(x))) for x in row['normal']]
            section=engine.section(point,normal)
            source_checks.append(dict(generated_node_id=row['node_id'],position=row['position'],exact_query_point=point,exact_query_normal=normal,section=section))
            if section.get('radius') is None:fail(section['status'],position=row['position'],location_kind='original_mesh_section',node_id=row['node_id'])
        selected_nodes={r['node_id'] for r in stations}
        for row in proxy['rows']:
            if row['node_id'] not in selected_nodes:fail('ORIGINAL_TRIANGLE_SECTION_UNMEASURED_AT_THIS_NODE',position=row['position'],location_kind='source_radius_coverage_gap',node_id=row['node_id'])
        fail('MODEL_TERMINAL_COVERAGE_NOT_CERTIFIED')
        fail('CONTINUOUS_ORIGINAL_MESH_CENTRELINE_INTERIORITY_UNPROVEN')
    except (ValueError,RuntimeError,AssertionError,OverflowError) as exc:
        fail('STAGE_REFUSED_OR_GUARDED',details=dict(stage=stage,error_type=type(exc).__name__,reason=str(exc)))
    if hashlib.sha256(mesh.vertices.tobytes()+mesh.faces.tobytes()).hexdigest()!=source_before:raise ValueError('ORIGINAL_MODEL_ARRAYS_CHANGED')
    return dict(schema_version='automatic-new-model-geometric-measurer-prototype-v001',status='PARTIAL_AUTOMATIC_MODEL_MEASUREMENT' if graph is not None else 'FAILED_AUTOMATIC_MODEL_MEASUREMENT',
        source_sha256=mesh.source_sha256,coordinate_units=coordinate_units,model_bounds=bounds,parameters=dict(grid_divisions=grid_divisions,pitch=h,grid_origin=[0,0,0],radius_validation_limit=5),
        candidate_centerline=graph,length=length,proxy_radius_curve=radius,original_triangle_sparse_radius_checks=source_checks,failures=failures,stages=stages,
        no_model_specific_node_ids_required=True,original_model_modified=False,original_model_array_sha_before_after=source_before,
        endpoint_extension_or_unknown_fill=False,reference_model_or_GT_used=False,original_mesh_radius_curve_complete=False,
        new_model_execution_contract_satisfied=graph is not None and length is not None and radius is not None,
        full_acceptance_satisfied=False,physical_length=None,physical_radius=None,actual_GT_accuracy=None)
