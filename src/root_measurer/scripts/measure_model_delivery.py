"""Integrated first-version CLI; no input graph or model-specific node IDs."""
import argparse,hashlib,json,sys
import ctypes,subprocess,threading,time
from pathlib import Path
P=Path(__file__).resolve().parents[1];W=P.parents[1];sys.path[:0]=[str(P/'src'),str(P/'scripts')];sys.dont_write_bytecode=True
from root_model_delivery import write_delivery
from measure_centerline_body import _output
from root_model_delivery_runtime import measure_integrated
from root_model_input_loader import file_sha

def _worker_memory(process):
    import psutil
    try:
        proc=psutil.Process(process.pid);info=proc.memory_info()
        private=getattr(info,'private',None)
        if private is None:
            try:private=proc.memory_full_info().uss
            except (psutil.AccessDenied,AttributeError):private=info.rss
        return dict(working_set=info.rss,private_bytes=private,
                    peak_working_set=getattr(info,'peak_wset',info.rss))
    except (psutil.NoSuchProcess,psutil.ZombieProcess):
        return dict(working_set=0,private_bytes=0,peak_working_set=0)


def _supervise_native(args,argv):
    from root_model_delivery_exports import file_reference,write_json,local_date
    prefix=args.output_prefix.resolve()
    root=next((p for p in prefix.parents if (p/'EXPERIMENT_LOG.md').is_file()),None)
    if root is None or not prefix.is_relative_to(root/'results'):
        raise ValueError('REGISTERED_EXPERIMENT_RESULTS_REQUIRED')
    if args.max_memory_mb<=0:
        raise ValueError('POSITIVE_WORKER_MEMORY_BUDGET_REQUIRED')
    prefix.parent.mkdir(parents=True,exist_ok=True)
    if any(prefix.with_name(prefix.name+s).exists() for s in
           ['_measurement.json','_progress.jsonl','_delivery.json','_supervision.json','_process.log']):
        raise ValueError('OUTPUT_PREFIX_ALREADY_EXISTS')
    with (root/'EXPERIMENT_LOG.md').open('a',encoding='utf-8') as stream:
        stream.write('\n'+local_date()+' PRELOG local native execution: '+
                     json.dumps(dict(output_prefix=prefix.relative_to(root).as_posix(),
                                     entity=str(args.entity) if args.entity else None,mesh=str(args.mesh) if args.mesh else None,stage=args.entity_stage,
                                     backend=args.native_backend,inherit_from=str(args.inherit_from)
                                     if args.inherit_from else None,
                                     max_memory_mb=args.max_memory_mb),ensure_ascii=False)+'\n')
    snapshot=root/'snapshot/native_model_delivery_run001/executions'
    snapshot.mkdir(parents=True,exist_ok=True)
    config=P/'configs/native_model_delivery_v001.json'
    execution_key=hashlib.sha256(str(prefix).encode('utf-8')).hexdigest()[:12]
    write_json(snapshot/(prefix.name+'_'+execution_key+'.json'),dict(config=file_reference(config,W),
               command=[sys.executable,'-B',str(Path(__file__).resolve())]+argv,
               backend=args.native_backend,stage=args.entity_stage,
               scale_mm_per_model_unit=args.scale_mm,max_memory_mb=args.max_memory_mb,
               source_files=[file_reference(P/name,W) for name in
                   ['scripts/measure_model_delivery.py','src/root_mesh_delivery_measurement.py','src/root_mesh_resolution.py','src/root_mesh_thinning_storage.py','src/root_mesh_body_index.py','src/root_exact_query_sections.py','src/root_branch_section_measurement.py','src/root_cubical_surface.py','src/root_closed_triangle_voxelization.py','src/root_native_entity_measurement.py',
                    'src/root_model_delivery_exports.py','src/root_reconstruction_measurement.py',
                    'src/root_cubical_graph_quotient.py','src/root_native_junction_fit.py',
                    'src/root_native_endpoint_retry.py','src/root_native_axis_order.py','src/root_model_parameters.py','src/root_local_network_policy.py']]))
    command=[sys.executable,'-B',str(Path(__file__).resolve())]+argv+['--_native-worker']
    started=time.monotonic();peak_ws=peak_private=0;reason=None
    process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                             text=True,encoding='utf-8',errors='replace',
                             env=dict(__import__('os').environ,PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',
                                      OPENBLAS_NUM_THREADS='6',OMP_NUM_THREADS='6'))
    log=prefix.with_name(prefix.name+'_process.log')
    def pump():
        with log.open('x',encoding='utf-8') as stream:
            for line in process.stdout:
                stream.write(line);stream.flush()
                print(line,end='',flush=True)
    reader=threading.Thread(target=pump,daemon=True);reader.start()
    try:
        while process.poll() is None:
            if args.cancel_file is not None and args.cancel_file.is_file():
                reason='USER_CANCELLED';process.terminate();break
            memory=_worker_memory(process)
            peak_ws=max(peak_ws,memory['peak_working_set'],memory['working_set'])
            peak_private=max(peak_private,memory['private_bytes'])
            if max(memory['working_set'],memory['private_bytes'])>args.max_memory_mb*1024**2:
                reason='WORKER_MEMORY_BUDGET_EXCEEDED';process.terminate();break
            time.sleep(.25)
        process.wait();reader.join(timeout=5)
    except KeyboardInterrupt:
        reason='USER_CANCELLED';process.terminate();process.wait();reader.join(timeout=5)
    outcome=dict(worker_exit_code=process.returncode,elapsed_seconds=time.monotonic()-started,
                 peak_working_set_bytes=peak_ws,peak_private_bytes=peak_private,
                 max_memory_mb=args.max_memory_mb,stop_reason=reason,
                 measurement_exists=prefix.with_name(prefix.name+'_measurement.json').is_file(),
                 delivery_exists=prefix.with_name(prefix.name+'_delivery.json').is_file())
    if outcome['measurement_exists']:
        outcome['measurement_status']=json.loads(prefix.with_name(prefix.name+'_measurement.json').read_text(encoding='utf-8'))['status']
    write_json(prefix.with_name(prefix.name+'_supervision.json'),outcome)
    with (root/'EXPERIMENT_LOG.md').open('a',encoding='utf-8') as stream:
        stream.write('\n'+local_date()+' RESULT local native execution: '+
                     json.dumps(dict(output_prefix=prefix.relative_to(root).as_posix(),**outcome))+'\n')
    if process.returncode!=0 or reason:
        raise RuntimeError(reason or 'NATIVE_WORKER_FAILED_SEE_PROCESS_LOG')
    if outcome.get('measurement_status') in ['NATIVE_MEASUREMENT_BLOCKED','MESH_MEASUREMENT_BLOCKED']:
        raise RuntimeError('NATIVE_MEASUREMENT_BLOCKED_SEE_REPORT')
    return json.loads(prefix.with_name(prefix.name+'_delivery.json').read_text(encoding='utf-8'))


def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    ap=argparse.ArgumentParser(description=__doc__)
    inputs=ap.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--mesh',type=Path);inputs.add_argument('--entity',type=Path)
    ap.add_argument('--output-prefix',type=Path,required=True)
    ap.add_argument('--coordinate-units');ap.add_argument('--grid-divisions',type=int,default=16)
    ap.add_argument('--execution-domain',choices=['local','registered_slurm_complete_source'],default='local');ap.add_argument('--reuse-measurement',type=Path)
    ap.add_argument('--entity-stage',choices=['regularized_negative_indices','source_negative_indices'],default='regularized_negative_indices')
    ap.add_argument('--mesh-workflow',action='store_true',help='supervised mesh target/challenge/original-surface workflow')
    ap.add_argument('--mesh-target-pitch',type=float)
    ap.add_argument('--mesh-target-reason',default='',help='input geometry/feature-scale basis, no GT parameters')
    ap.add_argument('--mesh-resume-target-from',type=Path,help='explicit diagnostic target checkpoint reuse; never called fresh')
    ap.add_argument('--scale-mm',type=float)
    ap.add_argument('--inherit-from',type=Path,help='M1 only: SHA-bound prior body/radius witness; reported as inherited')
    ap.add_argument('--resume-axis-from',type=Path,help='Reuse a bound prior automatic native axis and remeasure final sections')
    ap.add_argument('--native-backend',choices=['auto','existing_span','existing_guarded_lee'],default='auto')
    ap.add_argument('--max-memory-mb',type=int,default=8192)
    ap.add_argument('--cancel-file',type=Path,help=argparse.SUPPRESS)
    ap.add_argument('--_native-worker',action='store_true',help=argparse.SUPPRESS)
    a=ap.parse_args(argv)
    if a.scale_mm is not None and (not __import__('math').isfinite(a.scale_mm) or a.scale_mm<=0):
        ap.error('positive finite mm scale required')
    if not a.mesh_workflow and (a.mesh_target_pitch is not None or a.mesh_target_reason or a.mesh_resume_target_from is not None):
        ap.error('mesh target/cache options require --mesh-workflow')
    if a.mesh_workflow:
        if a.mesh is None or a.entity is not None or a.reuse_measurement or a.inherit_from or a.resume_axis_from:
            ap.error('mesh workflow requires only --mesh and a fresh measurement')
        if a.mesh_resume_target_from is not None:
            cache=a.mesh_resume_target_from.resolve()
            if not cache.is_relative_to(P/'experiments') or not cache.is_file():ap.error('mesh cache must be an existing project experiment measurement')
        if a._native_worker:
            from root_local_network_policy import disable_measurement_network
            disable_measurement_network()
            from root_mesh_delivery_measurement import run_mesh_worker
            return run_mesh_worker(a)
        return _supervise_native(a,argv)
    if a.entity is not None:
        if a.inherit_from and a.resume_axis_from:
            ap.error('--inherit-from and --resume-axis-from are mutually exclusive')
        if a.reuse_measurement is not None:
            ap.error('native entity reuse is explicit --inherit-from; do not relabel old mesh output')
        if a._native_worker:
            from root_local_network_policy import disable_measurement_network
            disable_measurement_network()
            from root_native_entity_measurement import run_native_worker
            return run_native_worker(a)
        return _supervise_native(a,argv)
    if a.inherit_from or a.resume_axis_from or a._native_worker:
        ap.error('native-only options require --entity')
    a.coordinate_units=a.coordinate_units or 'model_units_unspecified_not_mm'
    prefix=a.output_prefix.resolve();delivery=prefix.with_name(prefix.name+'_delivery.json');_output(delivery)
    root=next((p for p in prefix.parents if (p/'EXPERIMENT_LOG.md').is_file()),None);assert root is not None
    expected=file_sha(a.mesh)
    if a.reuse_measurement:
        measurement=a.reuse_measurement.resolve();r=json.loads(measurement.read_text(encoding='utf-8-sig'))
        if r['input_file_binding']['sha256']!=expected or r['input_file_binding']['bytes']!=a.mesh.stat().st_size:raise ValueError('REUSE_REQUIRES_IDENTICAL_INPUT_MODEL')
        if r['parameters']['grid_divisions']!=a.grid_divisions or r['coordinate_units']!=a.coordinate_units:raise ValueError('REUSE_REQUIRES_IDENTICAL_GLOBAL_PARAMETERS_AND_UNITS')
        if (r.get('automatic_generation_evidence',{}).get('integration_backend') or {}).get('name')!='fixed_existing_Source_span_current_neighborhood_Lee':raise ValueError('REUSE_BACKEND_MISMATCH_OLD_RESULT_PRESERVED')
    else:
        measurement=prefix.with_name(prefix.name+'_measurement.json');_,runtime=measure_integrated(['--mesh',str(a.mesh),'--output',str(measurement),'--coordinate-units',a.coordinate_units,'--grid-divisions',str(a.grid_divisions),'--execution-domain',a.execution_domain])
        with prefix.with_name(prefix.name+'_runtime.json').open('x',encoding='utf8') as f:json.dump(runtime,f,indent=2)
    data,out,html=write_delivery(measurement,prefix,root,model_path=a.mesh,case_label='固定跨度保护后端候选结果');assert file_sha(a.mesh)==expected
    print(json.dumps(dict(delivery=str(out),viewer=str(html),summary=data['summary']),ensure_ascii=False),flush=True)
    return data

if __name__=='__main__':main()
