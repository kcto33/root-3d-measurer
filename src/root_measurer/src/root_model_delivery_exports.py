"""Native entity delivery tables; graph length and coverage have separate denominators."""
import csv
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime,timedelta,timezone

import numpy as np

from root_graph_measurement_ledger import graph_measurement_ledger


def local_date():
    return datetime.now(timezone(timedelta(hours=8))).strftime('%Y-%m-%d')


def file_reference(path, relative_to=None):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return dict(path=str(path) if relative_to is None else
                __import__('os').path.relpath(path, relative_to).replace('\\', '/'),
                sha256=digest.hexdigest(), bytes=path.stat().st_size)


def write_json(path, data):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(data, stream, ensure_ascii=False, allow_nan=False, separators=(',', ':'))
        stream.write('\n')


def measurement_tables(report):
    """Return path/sample rows, using each edge exactly once even at junctions."""
    graph = report.get('candidate_centerline')
    if graph is None:
        return [], [], dict(C_axis_ref=None, C_radius_ref=None, V_radius_out=None,
                            reference_status='REFERENCE_NOT_USED', output_total_length=0.0)
    points = graph['points']
    ledger = graph_measurement_ledger(points, graph['adjacency'],
                                      coordinate_units=report['coordinate_units'])
    profile = report.get('radius_curve') or {}
    sampling = profile.get('sampling', 'edge_midpoint')
    rows = profile.get('rows', [])
    by_id = {r['node_id'] if sampling == 'node_halfedge' else
             tuple(sorted(r['node_ids'])): r for r in rows}
    scale = report.get('scale_mm_per_model_unit')
    paths, samples = [], []
    valid_weight, complete_segment_weight = 0.0, 0.0
    endpoint_status={r['node_id']:r['status'] for r in report.get('geometric_endpoints',[])}
    for branch in ledger['branches']:
        ids = branch['source_node_ids']
        arc = [0.0]
        for a, b in zip(ids, ids[1:]):
            arc.append(arc[-1] + math.dist(points[a], points[b]))
        closed = branch['closed']
        path_id = f"path{branch['branch_id']:05d}"
        path = dict(path_id=path_id, component_id=branch['component_id'], closed=closed,
                    node_ids=ids, length=branch['source_length'],
                    length_mm=None if scale is None else scale * branch['source_length'],
                    start_node=ids[0], end_node=ids[-1],
                    start_status='junction' if len(graph['adjacency'][ids[0]]) >= 3 else
                    'closed' if closed else endpoint_status.get(ids[0],'endpoint_unresolved'),
                    end_status='junction' if len(graph['adjacency'][ids[-1]]) >= 3 else
                    'closed' if closed else endpoint_status.get(ids[-1],'endpoint_unresolved'), samples=[])
        for at, (a, b) in enumerate(zip(ids, ids[1:])):
            edge_length = arc[at+1] - arc[at]
            stations = [(a, arc[at], edge_length/2), (b, arc[at+1], edge_length/2)] \
                if sampling == 'node_halfedge' else [(None, (arc[at]+arc[at+1])/2, edge_length)]
            edge_radii = []
            for node, location, weight in stations:
                row = by_id.get(node if node is not None else tuple(sorted((a, b))), {})
                radius = row.get('radius')
                available = radius is not None and row.get('available', True)
                if radius is not None and (not math.isfinite(radius) or radius < 0):
                    raise ValueError('FINITE_NONNEGATIVE_MEASURED_RADIUS_REQUIRED')
                if available:
                    valid_weight += weight
                edge_radii.append(available)
                position = row.get('position', points[node] if node is not None else
                                   [(x+y)/2 for x, y in zip(points[a], points[b])])
                sample = dict(path_id=path_id, component_id=branch['component_id'],
                              node_id=node, edge_nodes=[a, b], arc_position=location,
                              arc_weight=weight, position=position, radius=radius,
                              radius_mm=None if radius is None or scale is None else radius*scale,
                              diameter=None if radius is None else 2*radius,
                              diameter_mm=None if radius is None or scale is None else 2*radius*scale,
                              cross_section_area=None if radius is None else math.pi*radius**2,
                              status=row.get('status', 'NOT_MEASURED'), available=available,
                              normal=row.get('normal'), sampling=sampling)
                if report.get('entity_binding'):
                    from root_section_semantics import connected_section
                    section=connected_section(row)
                    sample.update(connected_section_radius=section['radius'],connected_section_area=section['area'],
                                  connected_section_classification=section['classification'],
                                  connected_section_radius_mm=None if section['radius'] is None or scale is None else section['radius']*scale)
                samples.append(sample)
                path['samples'].append(sample)
            if all(edge_radii):
                complete_segment_weight += edge_length
        paths.append(path)
        from root_model_parameters import radius_statistics
        stats=radius_statistics(path['samples'])
        path.update(radius_mean=stats['mean'],radius_min=stats['minimum'],radius_max=stats['maximum'],
                    valid_radius_arc_fraction=stats['valid_arc_fraction'],
                    radius_mean_mm=None if scale is None or stats['mean'] is None else scale*stats['mean'])
    total = ledger['measurement_graph_length']
    assert math.isclose(math.fsum(x['arc_weight'] for x in samples), total,
                        rel_tol=1e-12, abs_tol=1e-12)
    assert math.isclose(math.fsum(x['length'] for x in paths), total,
                        rel_tol=1e-12, abs_tol=1e-12)
    def fraction(value):
        if not total:return None
        ratio=value/total
        if not -1e-12<=ratio<=1+1e-12:raise AssertionError('COVERAGE_WEIGHT_OUTSIDE_DENOMINATOR')
        return min(1.,max(0.,ratio))
    coverage = dict(C_axis_ref=None, C_radius_ref=None,
                    reference_status='REFERENCE_NOT_USED', output_total_length=total,
                    output_valid_radius_arc_weight=valid_weight,
                    V_radius_out=fraction(valid_weight),
                    complete_segment_radius_fraction=fraction(complete_segment_weight),
                    sample_count=len(rows), available_sample_count=sum(
                        r.get('radius') is not None and r.get('available', True) for r in rows),
                    sampling=sampling,
                    definition='sample-supported output arc weight, not continuous or reference coverage')
    return paths, samples, coverage


def export_tables(report, prefix):
    prefix = Path(prefix)
    paths, samples, coverage = measurement_tables(report)
    refs = {}
    path_fields = ['path_id', 'component_id', 'closed', 'start_node', 'end_node',
                   'length', 'length_mm', 'coordinate_units', 'start_status', 'end_status',
                   'radius_mean','radius_mean_mm','radius_min','radius_max','valid_radius_arc_fraction',
                   'provenance_mode']
    station_fields = ['path_id', 'component_id', 'node_id', 'edge_start', 'edge_end',
                      'x', 'y', 'z', 'arc_position', 'arc_weight', 'radius', 'radius_mm',
                      'diameter','diameter_mm','cross_section_area','connected_section_radius',
                      'connected_section_radius_mm','connected_section_area','connected_section_classification',
                      'coordinate_units', 'available', 'status', 'sampling', 'provenance_mode']
    for kind, fields, values in [('paths', path_fields, paths), ('samples', station_fields, samples)]:
        target = prefix.with_name(prefix.name + '_' + kind + '.csv')
        with target.open('x', encoding='utf-8-sig', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for row in values:
                value = {k: row.get(k) for k in fields}
                if kind == 'samples':
                    value.update(zip(('x', 'y', 'z'), row['position']))
                    value.update(edge_start=row['edge_nodes'][0], edge_end=row['edge_nodes'][1])
                value.update(coordinate_units=report['coordinate_units'],
                             provenance_mode=report['provenance']['mode'])
                writer.writerow(value)
        refs[kind] = file_reference(target, prefix.parent)
    graph_path = prefix.with_name(prefix.name + '_graph.json')
    write_json(graph_path, dict(schema_version='root-delivery-graph-v001',
                               candidate_centerline=report.get('candidate_centerline'),
                               paths=[{k:v for k,v in p.items() if k != 'samples'} for p in paths],
                               coordinate_units=report['coordinate_units'],
                               scale_mm_per_model_unit=report.get('scale_mm_per_model_unit'),
                 source=report.get('entity_binding') or dict(input_kind='triangle_ply',file=report['input_file_binding']), coverage=coverage,
                               endpoint_status=(report.get('length') or {}).get('endpoint_status','endpoint_unresolved')))
    refs['graph'] = file_reference(graph_path, prefix.parent)
    return refs, paths, coverage


def native_surface_buffers(solid, cache_dir):
    """Display-only float32 boundary; cube-union measurement coordinates stay unchanged."""
    from root_reconstruction_measurement import boundary_faces
    cache = Path(cache_dir)
    key = hashlib.sha256(solid.cells.tobytes()+np.float64(solid.pitch).tobytes()).hexdigest()
    folder = cache / key
    manifest = folder / 'surface.json'
    if manifest.is_file():
        data = json.loads(manifest.read_text(encoding='utf-8'))
        for ref in data['buffers'].values():
            actual = file_reference(folder/ref['path'])
            if actual['sha256'] != ref['sha256'] or actual['bytes'] != ref['bytes']:
                raise ValueError('SURFACE_CACHE_CHANGED')
        return folder, data
    folder.mkdir(parents=True, exist_ok=False)
    centers, axes = boundary_faces(solid)
    vertices = np.repeat(centers[:,None,:], 4, axis=1)
    offsets = np.asarray([[-.5,-.5],[.5,-.5],[.5,.5],[-.5,.5]])
    for axis in range(3):
        selected = np.flatnonzero(axes == axis)
        other = [k for k in range(3) if k != axis]
        for j, dimension in enumerate(other):
            vertices[selected,:,dimension] += offsets[:,j]
    positions = (vertices.reshape(-1,3)*solid.pitch).astype('<f4')
    indices = (np.arange(len(centers),dtype='<u4')[:,None]*4 +
               np.asarray([0,1,2,0,2,3],dtype='<u4')).reshape(-1)
    buffers = {}
    for name, array in [('positions', positions), ('indices', indices)]:
        target = folder/(name+'.bin')
        with target.open('xb') as stream:
            array.tofile(stream)
        buffers[name] = file_reference(target, folder)
    data = dict(schema_version='cube-union-display-surface-v001', quads=len(centers),
                vertices=len(positions), triangles=len(indices)//3,
                buffers=buffers, source_cells_sha256=hashlib.sha256(solid.cells.tobytes()).hexdigest(),
                pitch=solid.pitch, measurement_geometry_modified=False,
                display_float32_only=True, all_boundary_faces_displayed=True)
    write_json(manifest, data)
    return folder, data


def write_native_delivery(report_path, prefix, experiment_root, solid):
    prefix, experiment_root = Path(prefix), Path(experiment_root)
    report = json.loads(Path(report_path).read_text(encoding='utf-8'))
    exports, paths, coverage = export_tables(report, prefix)
    surface_folder, surface = native_surface_buffers(
        solid, experiment_root/'results/model_delivery_local/geometry')
    surface = dict(surface, buffers={k:file_reference(surface_folder/v['path'], prefix.parent)
                                   for k,v in surface['buffers'].items()})
    from root_model_parameters import write_parameters
    parameters=write_parameters(report,prefix,solid=solid,boundary_face_count=surface['quads'])
    measurement_ref = file_reference(report_path, prefix.parent)
    result = dict(schema_version='root-native-delivery-v001', measurement=measurement_ref,
                  case_label=report.get('case_label', '原生实体测量'),
                  summary=dict(status=report['status'], candidate_length=report.get('length'),
                               paths=len(paths), coverage=coverage,
                               automatic_nodes=len((report.get('candidate_centerline') or {}).get('points', []))),
                  exports=exports, surface=surface, parameters=parameters,provenance=report['provenance'],
                  functional_acceptance='NOT_YET_PASSED', precision_acceptance='NOT_YET_PASSED')
    delivery = prefix.with_name(prefix.name+'_delivery.json')
    write_json(delivery, result)
    viewer = prefix.with_name(prefix.name+'_viewer.html')
    with viewer.open('x',encoding='utf-8') as stream:
        stream.write(Path(__file__).with_name('model_delivery_viewer.html').read_text(encoding='utf-8'))
    return result, delivery, viewer
