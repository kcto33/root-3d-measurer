"""Closed cube-union ray intervals; geometric diagnostics, never tip extension."""
import numpy as np


def cube_union_exit(cells, pitch, point, direction):
    cells = np.asarray(cells)
    p = np.asarray(point, float) / pitch
    d = np.asarray(direction, float)
    if cells.ndim != 2 or cells.shape[1] != 3 or cells.dtype.kind not in 'iu':
        raise ValueError('INTEGER_CELLS_REQUIRED')
    if not len(cells) or not np.isfinite(pitch) or pitch <= 0:
        raise ValueError('NONEMPTY_FINITE_POSITIVE_PITCH_REQUIRED')
    if p.shape != (3,) or d.shape != (3,) or not np.isfinite(p).all() or not np.isfinite(d).all() or not np.linalg.norm(d):
        raise ValueError('FINITE_NONZERO_RAY_REQUIRED')
    d = d / np.linalg.norm(d)
    lo = np.full(len(cells), -np.inf); hi = np.full(len(cells), np.inf)
    valid = np.ones(len(cells), bool)
    for a in range(3):
        if d[a] == 0:
            valid &= (cells[:, a] <= p[a]) & (p[a] <= cells[:, a]+1)
        else:
            t0 = (cells[:, a]-p[a])/d[a]; t1 = (cells[:, a]+1-p[a])/d[a]
            lo = np.maximum(lo, np.minimum(t0, t1)); hi = np.minimum(hi, np.maximum(t0, t1))
    # Cancellation error is measured in lattice coordinates, including origin.
    tol = 256*np.finfo(float).eps*max(1., float(np.max(np.abs(p))), float(np.max(np.abs(cells))))
    valid &= (hi >= lo-tol) & (hi >= -tol)
    spans = sorted(zip(lo[valid].tolist(), hi[valid].tolist()))
    if not any(a <= tol and b >= -tol for a, b in spans):
        raise ValueError('RAY_START_NOT_IN_CLOSED_SOLID')
    end = 0.; started = False; contacts = 0
    for a, b in spans:
        if not started:
            if a > tol: break
            if b >= -tol: started = True; end = max(0., b)
        elif a > end+tol:
            break
        else:
            if abs(a-end) <= tol and b > end+tol: contacts += 1
            end = max(end, b)
    boundary = p+end*d
    forward = np.floor(boundary).astype(np.int64)
    for a in range(3):
        nearest = round(boundary[a])
        if abs(boundary[a]-nearest) <= tol and d[a] != 0:
            forward[a] = nearest if d[a] > 0 else nearest-1
    if np.any(np.all(cells == forward, axis=1)):
        raise ValueError('NUMERIC_FORWARD_CELL_STILL_OCCUPIED')
    return dict(status='GEOMETRIC_CLOSED_CUBE_UNION_EXIT_NOT_PHYSICAL_TIP',
                boundary_distance=float(end*pitch), boundary_point=(boundary*pitch).tolist(),
                forward_cell=forward.tolist(), zero_width_ray_contacts=contacts,
                contact_ambiguity=bool(contacts), roundoff_tolerance_lattice=float(tol),
                physical_tip_position=None, adopted_extension=0., unknown_filled=False)
