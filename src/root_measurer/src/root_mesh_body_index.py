"""Complete segment AABB index over original cubes; same support predicate."""
import itertools
import math
from types import FunctionType,SimpleNamespace
import numpy as np
import root_section_centroid as centroid
import root_branch_section_measurement as associated


class SegmentSupportIndex:
    def __init__(self,cells,pitch,lookup=None,max_box_cells=4096):
        self.cells=cells;self.pitch=pitch;self.lookup=lookup if lookup is not None else {tuple(map(int,c)):i for i,c in enumerate(cells)}
        self.max_box_cells=max_box_cells;self.queries=self.fallbacks=self.selected_cells=0

    def covered(self,cells,h,a,b):
        self.queries+=1
        if cells is not self.cells or h!=self.pitch:
            self.fallbacks+=1;return centroid._covered(cells,h,a,b)
        # These are exactly the original predicate's arithmetic operations.
        # One-cell padding retains closed faces and any rounding-boundary cells.
        p=np.asarray(a)/h;d=(np.asarray(b)-a)/h
        low=np.floor(np.minimum(p,p+d)).astype(np.int64)-1
        high=np.floor(np.maximum(p,p+d)).astype(np.int64)+1
        if math.prod(map(int,high-low+1))>self.max_box_cells:
            self.fallbacks+=1;return centroid._covered(cells,h,a,b)
        ids=[self.lookup[key] for key in itertools.product(*(range(int(x),int(y)+1) for x,y in zip(low,high))) if key in self.lookup]
        selected=cells[sorted(ids)];self.selected_cells+=len(selected)
        return centroid._covered(selected,h,a,b)


def indexed_mesh_body(solid,points,adjacency,progress=None):
    from root_native_entity_measurement import NativeSectionIndex
    index=NativeSectionIndex(solid,associated.owned_section_class())
    support=SegmentSupportIndex(solid.cells,solid.pitch,lookup=index.occupied)
    class Engine:
        def __init__(self,source):
            if source is not solid:raise ValueError('SAME_SOURCE_REQUIRED_FOR_MESH_INDEX')
        def section(self,point,normal):return index.section(point,normal)
    namespace=SimpleNamespace(**vars(centroid));namespace._covered=support.covered
    base=associated.branch_section_centroid_measurement
    fn=FunctionType(base.__code__,dict(base.__globals__,original=namespace,owned_section_class=lambda:Engine),argdefs=base.__defaults__)
    fn.__kwdefaults__=base.__kwdefaults__
    result=fn(solid,points,adjacency,progress=progress)
    result['segment_support_index']=dict(queries=support.queries,full_scan_fallbacks=support.fallbacks,
        selected_source_cells=support.selected_cells,original_predicate_unchanged=True,
        complete_segment_candidate_AABB=True,original_section_centroid_codeobject=True)
    return result
