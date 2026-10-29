"""Measure an explicit native-input list sequentially; every failure stays listed."""
import argparse,csv,json,subprocess,sys,time
from pathlib import Path
P=Path(__file__).resolve().parents[1];W=Path.cwd()
sys.path.insert(0,str(P/'src'))
from root_model_delivery_exports import file_reference,write_json,local_date


def run(inputs,folder,scale=None,max_memory_mb=8192,cancel_file=None):
    folder=Path(folder).resolve()
    if folder.parent.name!='model_delivery_local' or not (folder.parents[2]/'ROOT_MEASURER_WORKSPACE.json').is_file():
        raise ValueError('EXPERIMENT_MODEL_DELIVERY_RUN_REQUIRED')
    folder.mkdir(exist_ok=False)
    experiment=folder.parents[2];records=[]
    with (experiment/'EXPERIMENT_LOG.md').open('a',encoding='utf-8') as f:
        f.write('\n'+local_date()+' PRELOG sequential native batch '+folder.relative_to(experiment).as_posix()+
                ': '+json.dumps(dict(inputs=[str(v) for v in inputs],scale=scale,memory_mb=max_memory_mb))+
                '; per-model supervision, measurements, exports and full batch denominator retained.\n')
    for i,path in enumerate(inputs):
        source=Path(path).resolve();prefix=folder/f'case{i:03d}'
        if cancel_file and Path(cancel_file).exists():
            records.append(dict(index=i,input=str(source),exit_code=None,prefix=prefix.relative_to(experiment).as_posix(),
                                complete_length=False,radius_valid_arc_fraction=None,status='CANCELLED_BEFORE_START'))
            continue
        command=[sys.executable,'-B',str(P/'scripts/measure_model_delivery.py'),'--entity',str(source),
                 '--output-prefix',str(prefix),'--max-memory-mb',str(max_memory_mb)]
        if cancel_file:command+=['--cancel-file',str(cancel_file)]
        if scale is not None:command+=['--scale-mm',str(scale)]
        event=dict(index=i,total=len(inputs),input=str(source),phase='measurement')
        print(json.dumps(event),flush=True)
        with (folder/'batch_progress.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(event)+'\n')
        with prefix.with_name(prefix.name+'_launcher.log').open('x',encoding='utf-8') as log:
            completed=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        measured=prefix.with_name(prefix.name+'_measurement.json');params=prefix.with_name(prefix.name+'_parameters.json')
        row=dict(index=i,input=str(source),exit_code=completed.returncode,prefix=prefix.relative_to(experiment).as_posix(),
                 complete_length=False,radius_valid_arc_fraction=None,status='FAILED')
        if params.is_file():
            value=json.loads(params.read_text(encoding='utf-8'))
            row.update(complete_length=value['complete_model_length_available'],
                       radius_valid_arc_fraction=value['coverage']['V_radius_out'],
                       connected_radius_valid_arc_fraction=value.get('connected_section_radius_model_units',{}).get('valid_arc_fraction'),
                       status='COMPLETE_AXIS_WITH_EXPLICIT_RADIUS_COVERAGE' if value['complete_model_length_available'] else 'INCOMPLETE_AXIS',
                       parameters=file_reference(params,W),measurement=file_reference(measured,W))
        records.append(row)
        print(json.dumps(row),flush=True)
    write_json(folder/'batch_summary.json',dict(schema='native-model-batch-v001',total=len(inputs),
        measurements=records,all_complete_lengths=all(r['complete_length'] for r in records),
        complete_radius_claim=False,missing_intervals_retained=True))
    with (folder/'batch_summary.csv').open('x',encoding='utf-8-sig',newline='') as f:
        fields=['index','input','exit_code','prefix','complete_length','radius_valid_arc_fraction','connected_radius_valid_arc_fraction','status']
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(records)
    with (experiment/'EXPERIMENT_LOG.md').open('a',encoding='utf-8') as f:
        f.write('\n'+local_date()+' RESULT batch '+folder.relative_to(experiment).as_posix()+
                ': '+str(sum(r['complete_length'] for r in records))+'/'+str(len(inputs))+' complete axes; per-model radius coverage and failures retained in batch_summary.json/CSV.\n')
    return records


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--manifest',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--scale-mm',type=float)
    ap.add_argument('--max-memory-mb',type=int,default=8192);ap.add_argument('--cancel-file',type=Path);a=ap.parse_args()
    value=json.loads(a.manifest.read_text(encoding='utf-8'));inputs=value['inputs']
    if not inputs or len(inputs)>1000:raise ValueError('BATCH_SIZE_1_TO_1000_REQUIRED')
    run([a.manifest.resolve().parent/path for path in inputs],a.output_dir,a.scale_mm,a.max_memory_mb,a.cancel_file)
