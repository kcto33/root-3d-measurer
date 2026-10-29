"""SHA-bound centerline body and area-radius measurement; no GT arguments."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
P=Path(__file__).resolve().parents[1]
W=P.parents[1]
sys.path.insert(0,str(P/'src'))
from root_reconstruction_measurement import ReconstructedSolid,load_reconstruction
from root_body_measurement_api import measure_centerline_body


def _sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()


def _json(path):
    def unique(items):
        out={}
        for key,value in items:
            if key in out:raise ValueError('DUPLICATE_REQUEST_OR_SOURCE_JSON_KEY')
            out[key]=value
        return out
    return json.loads(path.read_text(encoding='utf-8-sig'),object_pairs_hook=unique)


def _output(path):
    root=P/'experiments'
    if not path.is_relative_to(root):raise ValueError('REGISTERED_EXPERIMENT_OUTPUT_REQUIRED')
    parts=path.relative_to(root).parts
    if len(parts)<3 or parts[1] not in ('results','analysis','snapshot') or not (root/parts[0]/'EXPERIMENT_LOG.md').is_file():
        raise ValueError('REGISTERED_EXPERIMENT_OUTPUT_REQUIRED')
    if not path.parent.is_dir() or path.exists():raise ValueError('EXCLUSIVE_EXISTING_PARENT_OUTPUT_REQUIRED')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request',required=True);parser.add_argument('--request-sha256',required=True)
    parser.add_argument('--output',required=True);a=parser.parse_args(argv)
    req=Path(a.request).resolve();out=Path(a.output).resolve();_output(out)
    if _sha(req)!=a.request_sha256:raise ValueError('REQUEST_SHA_MISMATCH')
    request=_json(req)
    permitted={'coordinate_units','solid','kernel_graph','frozen_body','frozen_radius_profile','include_radius_profile'}
    if set(request)-permitted:raise ValueError('UNRECOGNIZED_REQUEST_FIELDS_NO_GT_ARGUMENTS')
    if not isinstance(request.get('include_radius_profile',False),bool):raise ValueError('BOOLEAN_RADIUS_REQUEST_REQUIRED')
    refs=[]
    def bound(ref):
        rel=Path(ref['workspace_relative_path'])
        if rel.is_absolute() or '..' in rel.parts:raise ValueError('WORKSPACE_RELATIVE_SOURCE_REQUIRED')
        path=(W/rel).resolve()
        if not path.is_relative_to(W):raise ValueError('WORKSPACE_SOURCE_REQUIRED')
        if _sha(path)!=ref['sha256'] or path.stat().st_size!=ref['bytes']:raise ValueError('SOURCE_SHA_OR_BYTES_MISMATCH')
        refs.append((ref,path));return path
    solidref=request['solid'];graphref=request['kernel_graph']
    sp=bound(solidref['input'])
    if solidref['kind']=='occupied_npz':
        with np.load(sp,allow_pickle=False) as z:
            if set(z.files)!={'cells','pitch'} or z['pitch'].shape!=():raise ValueError('CELLS_AND_SCALAR_PITCH_NPZ_REQUIRED')
            solid=ReconstructedSolid(z['cells'],float(z['pitch']),dict(source_sha256=solidref['input']['sha256']))
    elif solidref['kind']=='sparse_reconstruction':
        meta=_json(sp);crops=solidref['source_crops']
        if len(crops)!=len(meta['fields']):raise ValueError('FULL_SOURCE_CROP_BINDING_REQUIRED')
        caches={f['cache']['path']:(bound(r),r['sha256']) for f,r in zip(meta['fields'],crops)}
        solid=load_reconstruction(sp,solidref['input']['sha256'],caches)
    else:raise ValueError('SUPPORTED_SOLID_INPUT_KIND_REQUIRED')
    graph=_json(bound(graphref['input']))
    if graphref['kind']=='source_id_graph':
        with np.load(bound(graphref['source_skeleton_nodes']),allow_pickle=False) as z:cells=z['global_integer_cells']
        graph=dict(points=graph['graph_points'],adjacency=graph['graph_adjacency'],source_skeleton_cells=cells.tolist(),representation=graph['representation'])
    elif graphref['kind']!='tube_graph':raise ValueError('SUPPORTED_KERNEL_GRAPH_KIND_REQUIRED')
    cached_body=_json(bound(request['frozen_body'])) if 'frozen_body' in request else None
    cached_radius=_json(bound(request['frozen_radius_profile'])) if 'frozen_radius_profile' in request else None
    want=request.get('include_radius_profile',False)
    bodypath=out.with_name(out.stem+'_body.json');radiuspath=out.with_name(out.stem+'_radius.json')
    for path,needed in [(bodypath,cached_body is None),(radiuspath,want and cached_radius is None)]:
        if needed:_output(path)
    result=measure_centerline_body(solid,coordinate_units=request['coordinate_units'],kernel_graph=graph,
        frozen_body=cached_body,include_radius_profile=want,frozen_radius_profile=cached_radius)
    if _sha(req)!=a.request_sha256:raise ValueError('REQUEST_CHANGED_DURING_MEASUREMENT')
    for ref,path in refs:
        if _sha(path)!=ref['sha256']:raise ValueError('SOURCE_CHANGED_DURING_MEASUREMENT')
    def emit(path,value):
        with path.open('x',encoding='utf8') as f:json.dump(value,f,indent=2,ensure_ascii=False,allow_nan=False);f.write('\n')
        return dict(workspace_relative_path=path.relative_to(W).as_posix(),sha256=_sha(path),bytes=path.stat().st_size)
    bodyref=request.get('frozen_body')
    if cached_body is None and result['body_geometry'] is not None:bodyref=emit(bodypath,result['body_geometry'])
    radiusref=request.get('frozen_radius_profile')
    if want and cached_radius is None and result['area_equivalent_radius_profile'] is not None:radiusref=emit(radiuspath,result['area_equivalent_radius_profile'])
    profile=result['area_equivalent_radius_profile'];source=result['source_solid_measurement']
    receipt=dict(schema_version=result['schema_version'],request_sha256=a.request_sha256,request_reference=str(req.relative_to(W)) if req.is_relative_to(W) else None,
        units=result['units'],metrics=result['metrics'],source_geometric_volume=source['metrics']['occupied_geometric_volume'],
        source_graph_nodes=source['source_graph_ledger']['source_nodes'],source_graph_edges=source['source_graph_ledger']['source_edges'],
        body_geometry_reference=bodyref,area_equivalent_radius_profile_reference=radiusref,
        radius_definition=None if profile is None else profile['radius_definition'],radius_coverage=None if profile is None else profile['new_coverage'],
        radius_status_counts=None if profile is None else profile['status_counts'],
        cached_body_reused=result['frozen_body_reused'],source_geometry_modified=False,
        GT_reference_or_endpoint_span_accepted=False,physical_measurements_admitted=False,source_SHA_before_after=True)
    emit(out,receipt);return receipt


if __name__=='__main__':main()
