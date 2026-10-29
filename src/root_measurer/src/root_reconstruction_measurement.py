"""Read reconstructed cubical solids without inventing observation evidence.

The adapter consumes prediction metadata and its SHA-bound sparse crop files.
Its quantities describe the reconstructed solid; they are not biological
centerlines, physical root tips, or an observation-confidence qualification.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import math
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def _indices(values):
    a = np.asarray(values)
    if a.ndim != 2 or a.shape[1] != 3 or not np.issubdtype(a.dtype, np.integer):
        raise ValueError('INTEGER_CELL_INDICES_REQUIRED')
    if a.size and np.max(np.abs(a.astype(object))) > 2**50:
        raise ValueError('CELL_COORDINATE_RANGE')
    a = a.astype(np.int64, copy=True)
    if len(np.unique(a, axis=0)) != len(a):
        raise ValueError('DUPLICATE_OCCUPIED_CELL')
    return a


@dataclass(frozen=True)
class ReconstructedSolid:
    cells: np.ndarray
    pitch: float
    provenance: dict

    def __post_init__(self):
        a = _indices(self.cells)
        if not math.isfinite(self.pitch) or self.pitch <= 0:
            raise ValueError('FINITE_POSITIVE_PITCH_REQUIRED')
        a.setflags(write=False)
        object.__setattr__(self, 'cells', a)


def load_reconstruction(metadata_path, expected_sha256, cache_paths,
                        stage='regularized_negative_indices', lattice_atol=1e-14):
    """cache_paths maps each registered remote path to (local path, SHA).

    No reference geometry, centerline, radius or organ-label argument exists.
    The sparse NPZ omits observation weights; their absence remains explicit.
    """
    if stage not in ('source_negative_indices', 'regularized_negative_indices'):
        raise ValueError('UNSUPPORTED_RECONSTRUCTION_STAGE')
    metadata_path = Path(metadata_path)
    if file_sha(metadata_path) != expected_sha256:
        raise ValueError('PREDICTION_SHA_MISMATCH')
    metadata = json.loads(metadata_path.read_text(encoding='utf-8-sig'))
    if not metadata.get('fields') or metadata.get('duplicate_occupied_cells'):
        raise ValueError('INCOMPLETE_OR_DUPLICATED_RECONSTRUCTION')
    cells, receipts = [], []
    pitch = None
    source_count = 0
    for field in metadata['fields']:
        ref = field['cache']
        if ref['path'] not in cache_paths:
            raise ValueError('UNREGISTERED_RECONSTRUCTION_CROP')
        path, registered_sha = cache_paths[ref['path']]
        if registered_sha != ref['sha256'] or file_sha(path) != registered_sha:
            raise ValueError('CROP_SHA_MISMATCH')
        with np.load(path, allow_pickle=False) as z:
            if set(z.files) != {'source_negative_indices', 'regularized_negative_indices',
                                'pitch', 'origin', 'shape'}:
                raise ValueError('UNSUPPORTED_CROP_SCHEMA')
            h = float(z['pitch'])
            if not math.isfinite(h) or h <= 0 or (pitch is not None and h != pitch):
                raise ValueError('INCONSISTENT_CROP_PITCH')
            origin = np.asarray(z['origin'], float)
            shape = np.asarray(z['shape'])
            if origin.shape != (3,) or not np.isfinite(origin).all():
                raise ValueError('INVALID_CROP_ORIGIN')
            if shape.shape != (3,) or not np.issubdtype(shape.dtype, np.integer) or (shape <= 0).any():
                raise ValueError('INVALID_CROP_SHAPE')
            offset = np.rint(origin / h - .5)
            if (np.abs(offset) > 2**50).any() or not np.allclose(origin, (offset + .5) * h,
                                                                            atol=lattice_atol, rtol=0):
                raise ValueError('NONCANONICAL_HALF_CELL_LATTICE')
            source, regularized = (_indices(z[k]) for k in
                                  ('source_negative_indices', 'regularized_negative_indices'))
            for a in (source, regularized):
                if ((a < 0) | (a >= shape)).any():
                    raise ValueError('OCCUPIED_INDEX_OUTSIDE_CROP')
            # The declared one-face stage adds cells; never silently drop source material.
            if not set(map(tuple, source)).issubset(set(map(tuple, regularized))):
                raise ValueError('SOURCE_MATERIAL_REMOVED')
            selected = regularized if stage.startswith('regularized') else source
            cells.append(selected + offset.astype(np.int64))
            source_count += len(source)
            receipts.append(dict(path=str(path), sha256=registered_sha,
                                 cells=len(selected), shape=shape.tolist()))
            pitch = h
        if file_sha(path) != registered_sha:
            raise ValueError('CROP_CHANGED_DURING_READ')
    if source_count != metadata['all_source_negative_cells_accounted']:
        raise ValueError('INCOMPLETE_SOURCE_CELL_ACCOUNTING')
    if file_sha(metadata_path) != expected_sha256:
        raise ValueError('PREDICTION_CHANGED_DURING_READ')
    return ReconstructedSolid(np.concatenate(cells), pitch, dict(
        metadata_sha256=expected_sha256, crops=receipts, stage=stage,
        observed_weights_available=False, known_empty_exterior_available=False,
        unknown_filled=False, coordinates='canonical integer cell, center=(index+.5)*pitch',
        units='normalized_model_coordinates', geometry_modified=False,
        biological_admission=False))


def neighbor_graph(cells):
    """26-neighbor graph; each undirected edge stored exactly once."""
    cells = _indices(cells)
    if not len(cells):
        return np.empty((0, 2), np.int64), np.empty(0, np.int64)
    pairs = cKDTree(cells).query_pairs(1, p=np.inf, eps=0, output_type='ndarray')
    pairs = pairs.astype(np.int64, copy=False)
    graph = coo_matrix((np.ones(len(pairs), np.uint8), (pairs[:, 0], pairs[:, 1])),
                       shape=(len(cells), len(cells))).tocsr()
    _, labels = connected_components(graph, directed=False)
    return pairs, labels


def solid_summary(solid):
    from root_volume_topology import cubical_euler
    counts = cubical_euler(solid.cells)
    _, labels = neighbor_graph(solid.cells)
    sizes = sorted(np.bincount(labels).tolist(), reverse=True) if len(labels) else []
    return dict(counts, foreground_components=len(sizes), component_occupied_counts=sizes,
                volume=len(solid.cells) * solid.pitch**3, pitch=solid.pitch,
                foreground_connectivity=26, boundary='union of closed cubical cells',
                enclosed_cavities=None, tunnel_dimension=None,
                full_homology_recomputed=False, physical_topology_qualified=False,
                observation_confidence_available=False, biological_admission=False)


def boundary_faces(solid):
    """Face rectangles of the exact cube union, in grid coordinates."""
    cells = solid.cells
    if not len(cells):
        return np.empty((0, 3)), np.empty(0, np.int8)
    lo = cells.min(axis=0) - 1
    dims = cells.max(axis=0) - lo + 2
    if math.prod(map(int, dims)) > np.iinfo(np.int64).max:
        raise ValueError('CELL_KEY_OVERFLOW')
    encode = lambda a: ((a[:, 0] - lo[0]) * dims[1] + a[:, 1] - lo[1]) * dims[2] + a[:, 2] - lo[2]
    codes = np.sort(encode(cells))
    centers, axes = [], []
    for axis in range(3):
        for sign in (-1, 1):
            delta = np.zeros(3, np.int64)
            delta[axis] = sign
            key = encode(cells + delta)
            found = np.searchsorted(codes, key)
            present = (found < len(codes)) & (codes[np.minimum(found, len(codes)-1)] == key)
            face = cells[~present].astype(float) + .5
            face[:, axis] += .5 * sign
            centers.append(face)
            axes.append(np.full(len(face), axis, np.int8))
    return np.concatenate(centers), np.concatenate(axes)


def cubical_radii(solid, query_cells):
    """Exact distance from occupied cell centers to cube-union boundary faces.

    Unlike background-center EDT, this is a boundary distance. Neither is
    automatically an anatomical cross-section radius or an observed-boundary proof.
    """
    q = _indices(query_cells)
    if not set(map(tuple, q)).issubset(set(map(tuple, solid.cells))):
        raise ValueError('RADIUS_QUERY_NOT_IN_SOLID')
    if not len(q):
        return np.empty(0)
    centers, axes = boundary_faces(solid)
    tree = cKDTree(centers)
    points = q.astype(float) + .5
    near_distance, nearest = tree.query(points, eps=0)
    results = []
    for point, upper, nearest_id in zip(points, near_distance, nearest):
        # A face point lies at most sqrt(.5) from its center. This ball
        # includes every rectangle capable of beating the current upper bound.
        ids = tree.query_ball_point(point, np.nextafter(upper + math.sqrt(.5), np.inf), eps=0)
        delta = np.abs(centers[ids] - point)
        residual = np.maximum(delta - .5, 0)
        row = np.arange(len(ids))
        residual[row, axes[ids]] = delta[row, axes[ids]]
        results.append(math.sqrt(float(np.min(np.sum(residual**2, axis=1)))))
    return np.asarray(results) * solid.pitch


def cubical_ray_exit(solid, point, direction):
    """First exit of a supplied ray from the geometric cube union.

    This diagnostic accepts a direction, not a physical-tip model. Missing
    occupancy is NOT asserted to be observed empty. No extension is adopted.
    """
    point = np.asarray(point, float) / solid.pitch
    direction = np.asarray(direction, float)
    if point.shape != (3,) or direction.shape != (3,) or not np.isfinite(point).all() or not np.isfinite(direction).all():
        raise ValueError('FINITE_RAY_REQUIRED')
    norm = np.linalg.norm(direction)
    if not np.isfinite(norm) or norm == 0:
        raise ValueError('NONZERO_RAY_DIRECTION_REQUIRED')
    direction = direction / norm
    occupied = set(map(tuple, solid.cells))
    index = np.floor(point).astype(np.int64)
    if tuple(index) not in occupied:
        raise ValueError('RAY_START_NOT_IN_OCCUPIED_CELL')
    step = np.sign(direction).astype(np.int64)
    active = direction != 0
    entry = last_entry = 0.
    for visited in range(len(occupied) + 1):
        exits = np.full(3, np.inf)
        planes = index + (step > 0)
        exits[active] = (planes[active] - point[active]) / direction[active]
        leave = float(exits.min())
        if leave < entry - 1e-10:
            raise ValueError('NONMONOTONIC_GEOMETRIC_RAY')
        tied = np.abs(exits - leave) <= 16 * np.finfo(float).eps * max(1., abs(leave))
        next_index = index + step * tied
        if tuple(next_index) not in occupied:
            return dict(status='GEOMETRIC_CUBICAL_BOUNDARY_NOT_PHYSICAL_TIP',
                        boundary_distance=leave * solid.pitch,
                        boundary_point=((point + leave * direction) * solid.pitch).tolist(),
                        last_occupied_cell_interval=[entry * solid.pitch, leave * solid.pitch],
                        occupied_cells_traversed=visited + 1,
                        observation_available=False, physical_tip_position=None,
                        adopted_extension=0., unknown_filled=False)
        last_entry, entry, index = entry, leave, next_index
    raise RuntimeError('GEOMETRIC_RAY_ITERATION_GUARD')


def skeleton_diagnostic(solid, max_dense_cells=350000000):
    """Geometric Lee/graph diagnostic. Never returns an admitted root length."""
    from skimage.morphology import skeletonize
    from root_lattice_graph import reduce_lattice_graph
    if not len(solid.cells):
        return dict(status='EMPTY', biological_length=None), np.empty((0, 3), np.int64)
    lo = solid.cells.min(axis=0) - 1
    shape = solid.cells.max(axis=0) - lo + 2
    dense_count = math.prod(map(int, shape))
    if dense_count > max_dense_cells:
        raise ValueError('MEASUREMENT_DENSE_GUARD')
    occupied = np.zeros(tuple(shape), bool)
    occupied[tuple((solid.cells - lo).T)] = True
    skel = skeletonize(occupied, method='lee')
    nodes = np.argwhere(skel) + lo
    del occupied, skel
    pairs, labels = neighbor_graph(nodes)
    adjacency = [set() for _ in nodes]
    for a, b in pairs:
        adjacency[a].add(int(b)); adjacency[b].add(int(a))
    points = (nodes + .5) * solid.pitch
    if len(nodes):
        reduced, reduction = reduce_lattice_graph(points, adjacency, solid.pitch)
    else:
        reduced, reduction = [], {}
    reduced_pairs = np.asarray([(i, j) for i, ns in enumerate(reduced) for j in ns if i < j], np.int64).reshape(-1, 2)
    length = lambda edges: float(np.linalg.norm(points[edges[:, 0]] - points[edges[:, 1]], axis=1).sum())
    before_components = len(np.unique(neighbor_graph(solid.cells)[1]))
    after_components = len(np.unique(labels))
    degree = np.asarray([len(ns) for ns in reduced])
    report = dict(status='GEOMETRIC_DIAGNOSTIC_NOT_LENGTH_QUALIFICATION' if before_components == after_components else
                  'FAILED_SKELETON_COMPONENT_PRESERVATION', dense_cells=dense_count, nodes=len(nodes),
                  occupied_components=before_components, skeleton_components=after_components,
                  isolated_skeleton_nodes=int(np.count_nonzero(degree == 0)), endpoints=int(np.count_nonzero(degree == 1)),
                  raw_unique_edges=len(pairs), raw_unique_edge_length=length(pairs),
                  reduced_unique_edges=len(reduced_pairs), reduced_unique_edge_length=length(reduced_pairs),
                  reduced_graph_cycle_rank=len(reduced_pairs) - len(nodes) + after_components,
                  reduction=reduction, biological_length=None, physical_endpoints_qualified=False,
                  unknown_boundary_support_available=False, no_entity_geometry_changed=True,
                  no_component_pruning=True, units='normalized_model_coordinates')
    return report, nodes
