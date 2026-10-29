"""Create immutable synthetic native inputs; scoring axes never enter predictor metadata."""
import argparse
import hashlib
import json
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
PROJECT = Path(__file__).resolve().parents[1]
WORKSPACE = Path.cwd()
sys.path.insert(0, str(PROJECT / 'src'))
from root_model_delivery_exports import write_json, file_reference
from root_graph_measurement_ledger import graph_measurement_ledger

def voxelize_capsules(points, edges, radii, pitch):
    chunks = []
    for (a, b), radius in zip(edges, radii):
        p, q = (points[a], points[b])
        d = q - p
        length = float(d @ d)
        if length <= 0 or radius <= 0:
            raise ValueError('POSITIVE_SEGMENTS_AND_RADII_REQUIRED')
        lo = np.floor((np.minimum(p, q) - radius) / pitch - 0.5).astype(int)
        hi = np.ceil((np.maximum(p, q) + radius) / pitch - 0.5).astype(int)
        grid = np.stack(np.meshgrid(*(np.arange(x, y + 1) for x, y in zip(lo, hi)), indexing='ij'), axis=-1).reshape(-1, 3)
        v = (grid + 0.5) * pitch - p
        t = np.clip(v @ d / length, 0, 1)
        distance = np.einsum('ij,ij->i', v - t[:, None] * d, v - t[:, None] * d)
        chunks.append(grid[distance <= radius * radius])
    cells = np.unique(np.concatenate(chunks), axis=0)
    if not len(cells):
        raise ValueError('EMPTY_GENERATED_MODEL')
    return cells.astype(np.int64)

def reference_ray_exit(cells, h, point, direction):
    direction = np.asarray(direction, float)
    direction /= np.linalg.norm(direction)
    lo = cells * h
    hi = (cells + 1) * h
    enter = np.zeros(len(cells))
    leave = np.full(len(cells), np.inf)
    valid = np.ones(len(cells), bool)
    for a in range(3):
        if direction[a] == 0:
            valid &= (lo[:, a] <= point[a]) & (point[a] <= hi[:, a])
        else:
            x = (lo[:, a] - point[a]) / direction[a]
            y = (hi[:, a] - point[a]) / direction[a]
            enter = np.maximum(enter, np.minimum(x, y))
            leave = np.minimum(leave, np.maximum(x, y))
    valid &= (leave >= enter) & (leave >= 0)
    end = 0.0
    for a, b in sorted(zip(enter[valid], leave[valid])):
        if a > end + 1e-12:
            break
        end = max(end, float(b))
    return None if end == 0 else (np.asarray(point) + end * direction).tolist()

def construction(case):
    source = None
    if case == 'straight':
        points = np.array([[0, 0, 0], [0, 0, 1.0]], float)
        edges = [[0, 1]]
        radii = [0.1]
        pitch = 0.02
    elif case in ('curve', 'loop'):
        angles = np.linspace(0, math.pi / 2, 129) if case == 'curve' else np.linspace(0, 2 * math.pi, 256, endpoint=False)
        points = np.c_[np.cos(angles), np.sin(angles), np.zeros(len(angles))]
        edges = [[i, i + 1] for i in range(len(points) - 1)]
        if case == 'loop':
            edges.append([len(points) - 1, 0])
        radii = [0.1] * len(edges)
        pitch = 0.02 if case == 'curve' else 0.025
    elif case == 'branching_y':
        points = np.array([[0, 0, 0], [0, 0, -1], [math.sqrt(3) / 2, 0, 0.5], [-math.sqrt(3) / 2, 0, 0.5]])
        edges = [[0, 1], [0, 2], [0, 3]]
        radii = [0.1] * 3
        pitch = 0.02
    else:
        raise ValueError('UNKNOWN_DEMO_CASE')
    return (points, edges, list(map(float, radii)), pitch, source)

def build(case, base):
    folder = base / case
    if folder.exists():
        raise ValueError('ACCEPTANCE_CASE_ALREADY_EXISTS')
    points, edges, radii, h, source = construction(case)
    cells = voxelize_capsules(points, edges, radii, h)
    lo = cells.min(axis=0) - 1
    shape = cells.max(axis=0) - lo + 2
    folder.mkdir(parents=True)
    crop = folder / 'entity_crop0000.npz'
    local = cells - lo
    np.savez_compressed(crop, source_negative_indices=local, regularized_negative_indices=local, pitch=np.array(h), origin=(lo + 0.5) * h, shape=shape.astype(np.int64))
    ref = file_reference(crop, folder)
    metadata = dict(schema_version='generated-native-entity-v001', root=case, status='COMPLETE_GENERATED_FINITE_ENTITY', fields=[dict(cache=ref)], duplicate_occupied_cells=False, all_source_negative_cells_accounted=len(cells), storage=dict(pitch=h, occupied_voxels=len(cells), global_occupied_shape=(cells.max(axis=0) - cells.min(axis=0) + 1).tolist()), geometry='center-sampled union of closed geometric capsules; no observation confidence claim', no_closing_or_postprocessing=True)
    write_json(folder / 'entity_prediction.json', metadata)
    adjacency = [[] for _ in points]
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    final_points = points.copy()
    ends = []
    for i, row in enumerate(adjacency):
        if len(row) != 1:
            continue
        direction = points[i] - points[row[0]]
        end = reference_ray_exit(cells, h, points[i], direction)
        ends.append(dict(node_id=i, construction_axis_endpoint=points[i].tolist(), finite_entity_endpoint=end, status='RESOLVED' if end is not None else 'REFERENCE_ENDPOINT_UNRESOLVED'))
        if end is not None:
            final_points[i] = end
    ledger = graph_measurement_ledger(final_points, adjacency, coordinate_units='normalized_model_units')
    write_json(folder / 'reference_graph.json', dict(points=final_points.tolist(), adjacency=adjacency, construction_points=points.tolist(), edges=edges, construction_segment_radii=radii, pitch=h, endpoints=ends, source=source, max_radius=max(radii), reference_length=ledger['measurement_graph_length'], endpoint_definition='construction axis extended to first boundary of the exact generated cube union', radius_definition='area-equivalent section of generated cube union at reference position and tangent', predictor_must_not_read_this_file=True))
    from root_model_delivery_exports import local_date
    manifest = dict(case=case, created_date=local_date(), input=file_reference(folder / 'entity_prediction.json', WORKSPACE), crop=file_reference(crop, WORKSPACE), scorer_only=file_reference(folder / 'reference_graph.json', WORKSPACE), source=source, generator=file_reference(Path(__file__).resolve(), WORKSPACE), occupancy_cells=len(cells), pitch=h, center_sampling_declared=True, construction_discretization_is_separate_from_measurement_error=True, split='independent' if case in ('sparse_curved_seed5207', 'sparse_curved_seed5209', 'sparse_curved_seed5210', 'native_arbor_seed2602', 'native_arbor_seed2603', 'native_arbor_seed2604', 'native_arbor_seed2605', 'native_arbor_seed2606', 'native_arbor_seed2607') else 'development', prediction_opened=False)
    write_json(folder / 'dataset_manifest.json', manifest)
    with (folder / 'README.md').open('x', encoding='utf-8') as f:
        f.write('# ' + case + '\n\nGenerated ' + local_date() + ' from the saved construction, using center-sampled closed capsule unions. Native indices, pitch, origin and shape define the actual finite input. Source and delivered stages are identical; no closing is applied. The entity metadata binds only its NPZ. reference_graph.json is scorer-only and must not be passed to the predictor. Voxelization/finite-geometry differences from a continuous generating tube are reported separately. Generated synthetic demonstration geometry; not an independently held-out acceptance dataset. SHA-256, source and generator version are recorded in dataset_manifest.json.\n')
    return manifest
