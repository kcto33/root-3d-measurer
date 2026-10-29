"""Loopback-only native measurement workspace; jobs never invoke a shell."""
import hashlib
import json
import math
import os
import re
import secrets
import subprocess
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit

from root_model_delivery_exports import file_reference,write_json
from root_model_delivery_validation import validation_for
from root_model_delivery_downloads import display_validation,export_download,parameter_download
from root_model_delivery_reconstruction import reconstruction_comparison,reconstruction_csv

PROJECT=Path(__file__).resolve().parents[1]


def serve_workspace(experiment,port=0,open_browser=False,input_paths=()):
    experiment=Path(experiment).resolve()
    if not (experiment/'EXPERIMENT_LOG.md').is_file() or not (experiment/'ROOT_MEASURER_WORKSPACE.json').is_file():
        raise ValueError('ROOT_MEASURER_WORKSPACE_REQUIRED')
    base=experiment/'results/model_delivery_local'
    token=secrets.token_hex(24)
    entities={};jobs={};batches={};lock=threading.Lock()
    interaction_cache={};interaction_lock=threading.Lock()
    def register_entity(path):
        path=Path(path).resolve()
        if path.is_dir():path=path/'entity_prediction.json'
        if path.suffix.lower()=='.ply':
            if not path.is_file() or not 0<path.stat().st_size<=268435456:raise ValueError('MESH_SOURCE_BYTE_BUDGET')
            with path.open('rb') as stream:
                if stream.readline().strip()!=b'ply':raise ValueError('PLY_REQUIRED')
            key=hashlib.sha256(str(path).encode('utf-8')).hexdigest()[:16]
            entities[key]=dict(id=key,label=path.parent.name+' / '+path.stem,path=str(path),input_kind='triangle_ply',pitch=None,occupied_cells=None)
            return entities[key]
        data=json.loads(path.read_text(encoding='utf-8-sig'))
        if not data.get('fields') or 'all_source_negative_cells_accounted' not in data:
            raise ValueError('NATIVE_ENTITY_METADATA_REQUIRED')
        key=hashlib.sha256(str(path).encode('utf-8')).hexdigest()[:16]
        entities[key]=dict(id=key,label=data.get('root',path.parent.name),path=str(path),
                           input_kind='native_entity',
                           pitch=data.get('storage',{}).get('pitch'),
                           occupied_cells=data['all_source_negative_cells_accounted'])
        return entities[key]
    for path in input_paths:register_entity(path)
    def deliveries():
        found={}
        if base.is_dir():
            for path in sorted(base.glob('run*/*_delivery.json'),key=lambda p:p.stat().st_mtime,reverse=True):
                key=hashlib.sha256(path.relative_to(experiment).as_posix().encode()).hexdigest()[:16]
                try:
                    data=json.loads(path.read_text(encoding='utf-8'))
                    if data.get('schema_version')!='root-native-delivery-v001':continue
                    found[key]=dict(id=key,path=path,label=data['case_label'],
                                    provenance=data['provenance']['mode'],summary=data['summary'])
                except (ValueError,KeyError):continue
        return found
    def job_state(job):
        process=job['process'];prefix=job['prefix'];progress=prefix.with_name(prefix.name+'_progress.jsonl')
        event={}
        if progress.is_file():
            with progress.open('rb') as f:
                f.seek(max(0,progress.stat().st_size-8192));tail=f.read().splitlines()
            for line in reversed(tail):
                try:event=json.loads(line);break
                except ValueError:pass
        code=process.poll();supervision=prefix.with_name(prefix.name+'_supervision.json')
        outcome=json.loads(supervision.read_text(encoding='utf-8')) if supervision.is_file() else None
        measured=prefix.with_name(prefix.name+'_measurement.json')
        blocked=False
        if code is not None and measured.is_file():
            blocked=validation_for(prefix)['runtime_status'].startswith('BLOCKED')
        delivery=prefix.with_name(prefix.name+'_delivery.json')
        result_id=hashlib.sha256(delivery.relative_to(experiment).as_posix().encode()).hexdigest()[:16] if delivery.is_file() else None
        return dict(id=job['id'],entity_label=job['label'],running=code is None,result_id=result_id,
                    exit_code=code,progress=event,supervision=outcome,measurement_blocked=blocked,
                    elapsed_seconds=outcome['elapsed_seconds'] if code is not None and outcome else time.monotonic()-job['started'],
                    prefix=prefix.relative_to(experiment).as_posix())
    def result_routes(record):
        path=record['path'];prefix=path.with_name(path.name[:-len('_delivery.json')])
        data=json.loads(path.read_text(encoding='utf-8'))
        def bound(ref):
            p=(path.parent/ref['path']).resolve()
            if not p.is_relative_to(experiment):raise ValueError('RESULT_REFERENCE_OUTSIDE_EXPERIMENT')
            actual=file_reference(p)
            if any(actual[k]!=ref[k] for k in ['sha256','bytes']):raise ValueError('RESULT_REFERENCE_CHANGED')
            return p
        routes={'':(prefix.with_name(prefix.name+'_viewer.html'),'text/html; charset=utf-8'),
                'delivery.json':(path,'application/json'),
                'measurement.json':(bound(data['measurement']),'application/json'),
                'graph.json':(bound(data['exports']['graph']),'application/json'),
                'surface_positions.bin':(bound(data['surface']['buffers']['positions']),'application/octet-stream'),
                'surface_indices.bin':(bound(data['surface']['buffers']['indices']),'application/octet-stream')}
        for kind,ref in data['exports'].items():
            routes['export/'+kind]=(bound(ref),'application/json' if kind=='graph' else 'text/csv; charset=utf-8')
        return routes
    def batch_state(batch):
        folder=batch['folder'];progress=folder/'batch_progress.jsonl';event={}
        if progress.exists():
            for line in reversed(progress.read_text(encoding='utf-8').splitlines()):
                try:event=json.loads(line);break
                except ValueError:pass
        summary=folder/'batch_summary.json';code=batch['process'].poll()
        return dict(id=batch['id'],running=code is None,exit_code=code,total=batch['total'],progress=event,
                    summary=json.loads(summary.read_text(encoding='utf-8')) if summary.exists() else None,
                    elapsed_seconds=time.monotonic()-batch['started'])
    class Handler(BaseHTTPRequestHandler):
        def json_response(self,value,code=200):
            raw=json.dumps(value,ensure_ascii=False,allow_nan=False).encode('utf-8')
            self.send_response(code);self.send_header('Content-Type','application/json; charset=utf-8')
            self.send_header('Content-Length',str(len(raw)));self.send_header('Cache-Control','no-store')
            self.end_headers();self.wfile.write(raw)
        def send_file(self,path,kind,attachment=False,replace_token=False):
            if replace_token:
                raw=path.read_text(encoding='utf-8').replace('__LOCAL_TOKEN__',token).encode('utf-8')
            else:raw=None
            self.send_response(200);self.send_header('Content-Type',kind)
            self.send_header('Content-Security-Policy',"default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; connect-src 'self'; img-src 'self' data:; frame-src 'self'; object-src 'none'; base-uri 'self'")
            self.send_header('Content-Length',str(len(raw) if raw is not None else path.stat().st_size))
            self.send_header('Cache-Control','no-store')
            if attachment:self.send_header('Content-Disposition','attachment; filename="'+path.name+'"')
            self.end_headers()
            try:
                if raw is not None:self.wfile.write(raw)
                else:
                    with path.open('rb') as stream:
                        for block in iter(lambda:stream.read(1048576),b''):self.wfile.write(block)
            except (BrokenPipeError,ConnectionResetError,ConnectionAbortedError):pass
        def do_GET(self):
            path=urlsplit(self.path).path
            if path=='/':return self.send_file(PROJECT/'src/model_delivery_workspace.html','text/html; charset=utf-8',replace_token=True)
            if path=='/model_delivery_viewer.js':return self.send_file(PROJECT/'src/model_delivery_viewer.js','text/javascript; charset=utf-8')
            if path=='/model_delivery_interaction.js':return self.send_file(PROJECT/'src/model_delivery_interaction.js','text/javascript; charset=utf-8')
            if path=='/api/interaction-session':return self.json_response(dict(token=token))
            if path=='/api/state':
                records=deliveries()
                return self.json_response(dict(entities=list(entities.values()),
                    jobs=[job_state(j) for j in jobs.values()],
                    batches=[batch_state(b) for b in batches.values()],
                    deliveries=[{k:v for k,v in r.items() if k!='path'} for r in records.values()]))
            batch_match=re.fullmatch(r'/api/batches/([a-f0-9]{16})/(csv|json)',path)
            if batch_match and batch_match[1] in batches:
                target=batches[batch_match[1]]['folder']/('batch_summary.'+batch_match[2])
                if target.is_file():return self.send_file(target,'text/csv; charset=utf-8' if batch_match[2]=='csv' else 'application/json',attachment=True)
            match=re.fullmatch(r'/result/([a-f0-9]{16})/(.*)',path)
            if match:
                record=deliveries().get(match[1])
                if record and match[2]=='connected_validation.json':
                    prefix=record['path'].with_name(record['path'].name[:-len('_delivery.json')])
                    try:
                        from root_model_delivery_downloads import connected_validation_for
                        return self.json_response(connected_validation_for(prefix))
                    except (ValueError,OSError,KeyError) as exc:return self.json_response(dict(error=str(exc)),409)
                if record and match[2] in ['parameters.json','export/parameters','export/parameters_json']:
                    prefix=record['path'].with_name(record['path'].name[:-len('_delivery.json')])
                    try:raw,mime,filename=parameter_download(prefix,'csv' if match[2]=='export/parameters' else 'json')
                    except (ValueError,OSError,KeyError) as exc:return self.json_response(dict(error=str(exc)),409)
                    self.send_response(200);self.send_header('Content-Type',mime);self.send_header('Content-Length',str(len(raw)))
                    self.send_header('Cache-Control','no-store')
                    if match[2].startswith('export/'):self.send_header('Content-Disposition','attachment; filename="'+filename+'"')
                    self.end_headers();self.wfile.write(raw);return
                if record and match[2] in ['reconstruction_comparison.json','export/reconstruction']:
                    prefix=record['path'].with_name(record['path'].name[:-len('_delivery.json')])
                    try:
                        if match[2]=='reconstruction_comparison.json':return self.json_response(reconstruction_comparison(prefix))
                        raw=reconstruction_csv(prefix)
                    except (ValueError,OSError,KeyError) as exc:return self.json_response(dict(error=str(exc)),409)
                    self.send_response(200);self.send_header('Content-Type','text/csv; charset=utf-8')
                    self.send_header('Content-Length',str(len(raw)));self.send_header('Content-Disposition','attachment; filename="reconstruction_comparison.csv"')
                    self.end_headers();self.wfile.write(raw);return
                if record and match[2]=='validation.json':
                    prefix=record['path'].with_name(record['path'].name[:-len('_delivery.json')])
                    try:return self.json_response(display_validation(prefix))
                    except (ValueError,OSError,KeyError) as exc:return self.json_response(dict(error=str(exc)),409)
                if record and match[2] in ['export/paths','export/samples','export/graph','graph.json']:
                    prefix=record['path'].with_name(record['path'].name[:-len('_delivery.json')])
                    kind=match[2].split('/')[-1] if match[2]!='graph.json' else 'graph'
                    try:raw,mime,filename=export_download(prefix,kind)
                    except (ValueError,OSError,KeyError) as exc:return self.json_response(dict(error=str(exc)),409)
                    self.send_response(200);self.send_header('Content-Type',mime);self.send_header('Content-Length',str(len(raw)))
                    self.send_header('Cache-Control','no-store')
                    if match[2].startswith('export/'):self.send_header('Content-Disposition','attachment; filename="'+filename+'"')
                    self.end_headers();self.wfile.write(raw);return
                route=result_routes(record).get(match[2]) if record else None
                if route:return self.send_file(*route,attachment=match[2].startswith('export/'))
            self.send_error(404)
        def do_POST(self):
            if self.headers.get('X-Root-Local-Token')!=token:
                return self.json_response(dict(error='LOCAL_SESSION_TOKEN_REQUIRED'),403)
            length=int(self.headers.get('Content-Length','0'))
            if not 0<=length<=65536:return self.json_response(dict(error='REQUEST_SIZE'),400)
            try:
                data=json.loads(self.rfile.read(length) or b'{}')
                interaction=re.fullmatch(r'/result/([a-f0-9]{16})/interact',self.path)
                if interaction:
                    from root_model_interaction import ModelInteraction
                    from root_native_entity_measurement import check_entity_unchanged
                    record=deliveries().get(interaction[1])
                    if record is None:raise ValueError('RESULT_NOT_FOUND')
                    delivery=json.loads(record['path'].read_text(encoding='utf-8'))
                    measured=record['path'].parent/delivery['measurement']['path']
                    binding=file_reference(measured)
                    if binding['sha256']!=delivery['measurement']['sha256']:raise ValueError('MEASUREMENT_CHANGED')
                    with interaction_lock:
                        key=binding['sha256']
                        if key not in interaction_cache:
                            interaction_cache.clear()
                            interaction_cache[key]=ModelInteraction(json.loads(measured.read_text(encoding='utf-8')))
                        engine=interaction_cache[key]
                        if engine.report.get('entity_binding'):check_entity_unchanged(engine.report['entity_binding'])
                        if data.get('action')=='length':value=engine.length(data['intervals'])
                        elif data.get('action')=='radius':value=engine.radius(data)
                        else:raise ValueError('INTERACTION_ACTION_REQUIRED')
                        value['measurement']=binding
                        folder=record['path'].parent/'interactive_queries';folder.mkdir(exist_ok=True)
                        number=1
                        while (folder/f'{interaction[1]}_query{number:05d}.json').exists():number+=1
                        query=folder/f'{interaction[1]}_query{number:05d}.json'
                        write_json(query,dict(request=data,result=value,measurement=binding))
                        value['query_record']=file_reference(query,experiment)
                        with (experiment/'EXPERIMENT_LOG.md').open('a',encoding='utf-8') as f:
                            f.write('\n'+time.strftime('%Y-%m-%d')+' interactive '+str(data['action'])+': '+query.relative_to(experiment).as_posix()+'; original measurement/source bound; fresh query, no stored-radius interpolation.\n')
                    return self.json_response(value)
                if self.path=='/api/register-entity':
                    path=data.get('path')
                    if not isinstance(path,str) or not path.strip():raise ValueError('ENTITY_PATH_REQUIRED')
                    return self.json_response(dict(entity=register_entity(path.strip().strip('"'))))
                if self.path=='/api/select-entity':
                    startup=None
                    if os.name=='nt':
                        startup=subprocess.STARTUPINFO();startup.dwFlags|=subprocess.STARTF_USESHOWWINDOW
                        startup.wShowWindow=1
                    response=subprocess.run([sys.executable,'-B',str(PROJECT/'scripts/serve_model_delivery.py'),
                                             '--pick-entity'],capture_output=True,text=True,encoding='utf-8',
                                            startupinfo=startup,
                                            creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
                    if response.returncode:raise RuntimeError('ENTITY_FILE_PICKER_FAILED')
                    chosen=json.loads(response.stdout)['path']
                    return self.json_response(dict(entity=register_entity(chosen) if chosen else None))
                if self.path=='/api/batches':
                    with lock:
                        if any(j['process'].poll() is None for j in [*jobs.values(),*batches.values()]):raise ValueError('ONE_LOCAL_MEASUREMENT_AT_A_TIME')
                        inputs=data.get('inputs')
                        if not isinstance(inputs,list) or not 1<=len(inputs)<=1000 or not all(isinstance(p,str) and p.strip() for p in inputs):raise ValueError('EXPLICIT_NATIVE_INPUT_LIST_REQUIRED')
                        records=[register_entity(p.strip().strip('"')) for p in inputs]
                        if any(r['input_kind']!='native_entity' for r in records):raise ValueError('BATCH_NATIVE_ENTITIES_REQUIRED')
                        sources=[r['path'] for r in records]
                        if len(set(sources))!=len(sources):raise ValueError('DUPLICATE_BATCH_INPUT')
                        scale=data.get('scale_mm')
                        if scale is not None and (isinstance(scale,bool) or not math.isfinite(float(scale)) or float(scale)<=0):raise ValueError('POSITIVE_MM_SCALE_REQUIRED')
                        key=secrets.token_hex(8);registered=experiment/'snapshot/batch_requests';registered.mkdir(exist_ok=True)
                        manifest=registered/(key+'.json');cancel=registered/(key+'_cancel.json')
                        run=1
                        while (base/f'run{run:03d}').exists():run+=1
                        folder=base/f'run{run:03d}'
                        with (experiment/'EXPERIMENT_LOG.md').open('a',encoding='utf-8') as f:
                            f.write('\n'+time.strftime('%Y-%m-%d')+' PRELOG UI batch: '+manifest.relative_to(experiment).as_posix()+' inputs; '+folder.relative_to(experiment).as_posix()+' outputs; cancellation marker in snapshot/batch_requests; results/model_delivery_local/workspace_batch_'+key+'.log records launcher.\n')
                        write_json(manifest,dict(inputs=sources))
                        cmd=[sys.executable,'-B',str(PROJECT/'scripts/measure_model_batch.py'),'--manifest',str(manifest),'--output-dir',str(folder),'--cancel-file',str(cancel)]
                        if scale is not None:cmd+=['--scale-mm',str(float(scale))]
                        with (base/('workspace_batch_'+key+'.log')).open('x',encoding='utf-8') as log:
                            process=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,
                                env=dict(os.environ,PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1'),creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
                        batches[key]=dict(id=key,process=process,folder=folder,cancel=cancel,total=len(sources),started=time.monotonic())
                        return self.json_response(dict(batch_id=key))
                if self.path=='/api/jobs':
                    with lock:
                        if any(j['process'].poll() is None for j in [*jobs.values(),*batches.values()]):
                            raise ValueError('ONE_LOCAL_MEASUREMENT_AT_A_TIME')
                        entity=entities[data['entity_id']]
                        scale=data.get('scale_mm')
                        if scale is not None and (isinstance(scale,bool) or not math.isfinite(float(scale)) or float(scale)<=0):
                            raise ValueError('POSITIVE_MM_SCALE_REQUIRED')
                        target=data.get('mesh_target_pitch');reason=data.get('mesh_target_reason','')
                        if entity.get('input_kind')=='triangle_ply':
                            if target is not None and (isinstance(target,bool) or not math.isfinite(float(target)) or float(target)<=0):raise ValueError('POSITIVE_INPUT_DERIVED_TARGET_PITCH_REQUIRED')
                            if not isinstance(reason,str) or len(reason)>2000:raise ValueError('MESH_FEATURE_SCALE_REASON_REQUIRED')
                        run=1
                        while (base/f'run{run:03d}').exists():run+=1
                        folder=base/f'run{run:03d}';folder.mkdir(parents=True)
                        name=re.sub('[^a-zA-Z0-9_]+','_',entity['label'])[:60] or 'entity'
                        prefix=folder/(name+'_fresh');cancel=folder/'cancel.json'
                        cmd=[sys.executable,'-B',str(PROJECT/'scripts/measure_model_delivery.py'),
                             '--mesh' if entity.get('input_kind')=='triangle_ply' else '--entity',entity['path'],'--output-prefix',str(prefix),
                             '--native-backend','auto','--cancel-file',str(cancel)]
                        if entity.get('input_kind')=='triangle_ply':
                            cmd+=['--mesh-workflow','--mesh-target-reason',reason]
                            if target is not None:cmd+=['--mesh-target-pitch',str(float(target))]
                        if scale is not None:cmd+=['--scale-mm',str(float(scale))]
                        with (folder/'launcher.log').open('x',encoding='utf-8') as log:
                            process=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,
                                env=dict(os.environ,PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1'),
                                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
                        key=secrets.token_hex(8)
                        jobs[key]=dict(id=key,process=process,prefix=prefix,cancel=cancel,
                                       label=entity['label'],started=time.monotonic())
                        return self.json_response(dict(job_id=key))
                if self.path=='/api/cancel':
                    job=(batches if data.get('batch_id') else jobs)[data.get('batch_id') or data['job_id']]
                    if job['process'].poll() is None and not job['cancel'].exists():
                        write_json(job['cancel'],dict(reason='USER_CANCELLED_FROM_LOCAL_WORKSPACE'))
                    return self.json_response(dict(cancel_requested=True))
                self.send_error(404)
            except (ValueError,KeyError,RuntimeError,OSError) as exc:
                self.json_response(dict(error=str(exc)),400)
        def log_message(self,*_):pass
    server=ThreadingHTTPServer(('127.0.0.1',port),Handler)
    url='http://127.0.0.1:%d/'%server.server_address[1]
    print(url,flush=True)
    if open_browser:
        import webbrowser
        webbrowser.open(url)
    server.serve_forever()
