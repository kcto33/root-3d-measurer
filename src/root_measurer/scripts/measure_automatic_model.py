"""New triangle PLY -> automatic graph, length, sampled Source radii and failures."""
import argparse,hashlib,json,sys,time
from pathlib import Path
P=Path(__file__).resolve().parents[1];sys.path[:0]=[str(P/'src'),str(P/'scripts')];sys.dont_write_bytecode=True
from root_triangle_measurement_input import read_triangle_ply_bytes
from root_model_input_loader import load_model_input,file_sha
from root_automatic_model_measurement import measure_model
from root_automatic_source_profile import original_mesh_edge_profile
from measure_centerline_body import _output

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--mesh',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--coordinate-units',default='model_units_unspecified_not_mm');ap.add_argument('--grid-divisions',type=int,default=16)
    ap.add_argument('--execution-domain',choices=['local','registered_slurm_complete_source'],default='local');args=ap.parse_args(argv)
    output=args.output.resolve();_output(output);source=args.mesh.resolve();started=time.monotonic()
    if source.stat().st_size>268435456:raise ValueError('SOURCE_BYTE_GUARD')
    size=source.stat().st_size;digest=file_sha(source)
    mesh=load_model_input(source,execution_domain=args.execution_domain,expected_sha=digest)
    prototype=measure_model(mesh,coordinate_units=args.coordinate_units,grid_divisions=args.grid_divisions,execution_domain=args.execution_domain)
    failures=list(prototype['failures']);profile=None;graph=prototype['candidate_centerline']
    if graph is not None:
        try:
            profile=original_mesh_edge_profile(mesh,graph['points'],graph['adjacency'],execution_domain=args.execution_domain)
            for row in profile['rows']:
                if not row['available']:failures.append(dict(code=row['status'],location_kind='generated_edge_midpoint',generated_edge_id=row['generated_edge_id'],position=row['position'],details=row.get('original_mesh_section')))
        except (ValueError,RuntimeError,AssertionError,OverflowError) as exc:
            failures.append(dict(code='ORIGINAL_MESH_RADIUS_PROFILE_GUARDED',location_kind='whole_model_bounds',model_bounds=prototype['model_bounds'],position=None,details=dict(error_type=type(exc).__name__,reason=str(exc))))
    assert file_sha(source)==digest and source.stat().st_size==size
    complete_execution=bool(prototype['new_model_execution_contract_satisfied'] and profile is not None)
    result=dict(schema_version='new-model-automatic-geometric-measurer-v001',status='AUTOMATIC_OUTPUT_WITH_EXPLICIT_UNCERTAINTIES' if complete_execution else 'PARTIAL_OR_FAILED_AUTOMATIC_MODEL',
        input_file_binding=dict(sha256=digest,bytes=size,format='triangle_ply',original_file_unchanged=True),coordinate_units=args.coordinate_units,parameters=prototype['parameters'],
        candidate_centerline=graph,length=prototype['length'],radius_curve=profile,voxel_proxy_radius_curve=prototype['proxy_radius_curve'],failures=failures,
        automatic_generation_evidence=prototype,new_model_without_manual_node_ids=complete_execution,all_radius_station_outcomes_reported=profile is not None,
        all_original_triangle_station_radii_available=bool(profile is not None and profile['station_denominator']>0 and profile['available_radii']==profile['station_denominator']),
        centreline_accuracy_or_uniqueness_qualified=False,model_terminal_coverage_qualified=False,continuous_original_mesh_interiority_qualified=False,
        physical_length=None,physical_radius=None,actual_GT_accuracy=None,overall_scientific_goal_complete=False,worker_seconds=time.monotonic()-started)
    with output.open('x',encoding='utf8') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps(dict(status=result['status'],output=output.name,automatic_nodes=len(graph['points']) if graph else 0,source_radius_stations=profile['station_denominator'] if profile else 0,available_Source_radii=profile['available_radii'] if profile else 0)),flush=True)
    return result

if __name__=='__main__':main()
