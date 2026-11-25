"""Ordered independent plane queries; each spawned worker indexes its solid once."""
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import numpy as np

_engine = _solid = _slabs = _kind = None


def _initialize(solid, kind):
    global _engine, _solid, _slabs, _kind
    from root_native_entity_measurement import NativeSectionIndex
    if kind == 'body':
        from root_branch_section_measurement import owned_section_class
        kernel = owned_section_class()
    else:
        from root_cubical_disk_section import interior_boundary_section_class
        kernel = interior_boundary_section_class()
    _solid = solid
    _engine = NativeSectionIndex(solid, kernel)
    _slabs = {}
    _kind = kind


def section_query(engine, point, normal, slabs, kind):
    from root_native_entity_measurement import axis_aligned_section
    if kind=='body':
        # The legacy body kernel rejects grid-face and cut-cell-boundary queries.
        # The general axis fast path joins those squares, so use the unchanged
        # kernel near lattice boundaries rather than broadening body admission.
        p=np.asarray(point,float)/engine.solid.pitch
        tol=256*np.finfo(float).eps*max(1.,float(abs(p).max()),engine.maximum)
        if np.any(abs(p-np.rint(p))<=tol):return engine.section(point,normal)
    result=axis_aligned_section(engine.solid,point,normal,slabs)
    return engine.section(point,normal) if result is None else result


def _query(query):
    point, normal = query
    return section_query(_engine,point,normal,_slabs,_kind)


def sections(solid, queries, workers, kind='disk'):
    """Executor.map preserves source order, including missing/failed sections."""
    if not queries:
        return []
    if workers == 1:
        # Keep serial callers independent of module-level worker state.
        from root_native_entity_measurement import NativeSectionIndex, axis_aligned_section
        if kind == 'body':
            from root_branch_section_measurement import owned_section_class
            kernel = owned_section_class()
        else:
            from root_cubical_disk_section import interior_boundary_section_class
            kernel = interior_boundary_section_class()
        engine = NativeSectionIndex(solid, kernel); slabs = {}
        return [section_query(engine,p,n,slabs,kind) for p,n in queries]
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn'),
                             initializer=_initialize, initargs=(solid, kind)) as executor:
        return list(executor.map(_query, queries, chunksize=max(1, min(32, len(queries)//(workers*4)))))
