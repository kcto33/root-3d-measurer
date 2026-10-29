"""Read native observation states without confusing absence with empty material."""
from dataclasses import dataclass
import math
import numpy as np
from root_reconstruction_measurement import ReconstructedSolid, _indices

KEY = np.dtype([('x','<i8'),('y','<i8'),('z','<i8')])


def _keys(cells):
    return np.ascontiguousarray(cells,dtype=np.int64).view(KEY).reshape(-1)


@dataclass(frozen=True)
class ObservedMeasurementField:
    indices: np.ndarray
    tsdf: np.ndarray
    weights: np.ndarray
    pitch: float
    frames: int
    provenance: dict

    def __post_init__(self):
        cells=_indices(self.indices)
        if len(cells)>350000000:
            raise ValueError('OBSERVED_FIELD_CELL_GUARD')
        if not math.isfinite(self.pitch) or self.pitch<=0:
            raise ValueError('FINITE_POSITIVE_PITCH_REQUIRED')
        if isinstance(self.frames,bool) or not isinstance(self.frames,int) or self.frames<=0:
            raise ValueError('REGISTERED_POSITIVE_FRAME_COUNT_REQUIRED')
        values=np.asarray(self.tsdf);weights=np.asarray(self.weights)
        if values.shape!=(len(cells),) or weights.shape!=values.shape:
            raise ValueError('OBSERVED_FIELD_ARRAY_SHAPES')
        if not np.issubdtype(values.dtype,np.number) or np.iscomplexobj(values) or not np.isfinite(values).all():
            raise ValueError('FINITE_REAL_TSDF_REQUIRED')
        if not np.issubdtype(weights.dtype,np.number) or np.iscomplexobj(weights) or not np.isfinite(weights).all():
            raise ValueError('FINITE_REAL_WEIGHTS_REQUIRED')
        if np.any((values < -1)|(values > 1)):
            raise ValueError('NATIVE_TSDF_RANGE')
        if np.any((weights<=0)|(weights>self.frames)|(weights!=np.floor(weights))):
            raise ValueError('POSITIVE_NATIVE_INTEGRATION_COUNTS_REQUIRED')
        values=np.array(values,copy=True);weights=np.array(weights,copy=True)
        for a in (cells,values,weights):a.setflags(write=False)
        object.__setattr__(self,'indices',cells);object.__setattr__(self,'tsdf',values)
        object.__setattr__(self,'weights',weights)

    def summary(self):
        negative=int(np.count_nonzero(self.tsdf<0))
        return dict(observed_cells=len(self.indices),known_negative_cells=negative,
            known_nonnegative_cells=len(self.indices)-negative,frames=self.frames,pitch=self.pitch,
            known_negative_geometric_volume=negative*self.pitch**3,
            all_observed_weights_retained=True,unknown_filled=0,
            unqueried_unknown_cell_count=None,semantic_material_volume=None,
            physical_volume=None,physical_units_qualified=False)

    def known_negative_domain(self):
        return ReconstructedSolid(self.indices[self.tsdf<0],self.pitch,dict(
            source=self.provenance,stage='NATIVE_KNOWN_NEGATIVE_OBSERVATION_DOMAIN',
            observed_weights_available=True,unknown_filled=False,
            outside_negative_domain_not_assumed_observed_empty=True,
            semantic_R16_entity=False,biological_admission=False,
            units='declared_scene_coordinate_units_not_physical_mm'))

    def query(self,cells):
        # Repeated query IDs are allowed; source IDs remain unique.
        query=np.asarray(cells)
        if query.ndim!=2 or query.shape[1]!=3 or not np.issubdtype(query.dtype,np.integer):
            raise ValueError('INTEGER_QUERY_INDICES_REQUIRED')
        if query.size and np.max(np.abs(query.astype(object)))>2**50:
            raise ValueError('QUERY_COORDINATE_RANGE')
        source=_keys(self.indices);order=np.argsort(source);keys=source[order]
        wanted=_keys(query);at=np.searchsorted(keys,wanted)
        found=at<len(keys)
        if len(keys):found &= keys[np.minimum(at,len(keys)-1)]==wanted
        result=[]
        for i,okay in enumerate(found):
            if not okay:
                result.append(dict(index=query[i].tolist(),state='UNKNOWN',tsdf=None,weight=None))
            else:
                row=order[at[i]];value=float(self.tsdf[row])
                result.append(dict(index=query[i].tolist(),state='KNOWN_NEGATIVE' if value<0 else 'KNOWN_NONNEGATIVE',
                                   tsdf=value,weight=int(self.weights[row])))
        return result
