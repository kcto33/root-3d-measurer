"""Fixed existing Source-span thinning policy for automatic integrated delivery."""
from types import FunctionType
import numpy as np
from scipy import ndimage
import root_tube_measurement as tube
import root_automatic_model_measurement as automatic
from root_guarded_thinning import _input,CONNECT26
from root_span_guarded_thinning import span_anchored_component

def span_backend(volume,method='lee',*,max_dense_cells=350000000,receipt=None):
    if method!='lee':raise ValueError('ONLY_EXISTING_SOURCE_SPAN_LEE_SUPPORTED')
    source=_input(volume,max_dense_cells);labels,count=ndimage.label(source,structure=CONNECT26);out=np.zeros_like(source);rows=[]
    for component,box in enumerate(ndimage.find_objects(labels),1):
        part=labels[box]==component;result,audit=span_anchored_component(part,max_dense_cells=max_dense_cells);out[box]|=result
        rows.append(dict(component_id=component,source_voxels=int(part.sum()),derived_voxels=int(result.sum()),existing_span_audit=audit))
    assert not np.any(out&~source) and ndimage.label(out,structure=CONNECT26)[1]==count
    if receipt is not None:receipt.update(status='EXISTING_SOURCE_SPAN_CURRENT_NEIGHBORHOOD_BACKEND',source_components=count,derived_components=count,components=rows,source_material_unchanged=True,compiled_Lee_reexecuted=False,missing_compiled_components=[],physical_endpoint_inference=False)
    return out

def measure_span_tube(solid,*,max_dense_cells=64000000):
    fn=FunctionType(tube.measure_tube_solid.__code__,dict(tube.measure_tube_solid.__globals__,guarded_skeletonize=span_backend),tube.measure_tube_solid.__name__,tube.measure_tube_solid.__defaults__);fn.__kwdefaults__=tube.measure_tube_solid.__kwdefaults__
    assert fn.__code__ is tube.measure_tube_solid.__code__
    r=fn(solid,max_dense_cells=max_dense_cells);r['integration_binding']=dict(original_tube_pipeline_codeobject_exact=True,only_existing_span_thinning_backend_bound=True,original_module_globals_modified=False);return r

def measure_model_span(mesh,**kwargs):
    fn=FunctionType(automatic.measure_model.__code__,dict(automatic.measure_model.__globals__,measure_tube_solid=measure_span_tube),automatic.measure_model.__name__,automatic.measure_model.__defaults__);fn.__kwdefaults__=automatic.measure_model.__kwdefaults__
    assert fn.__code__ is automatic.measure_model.__code__
    r=fn(mesh,**kwargs);r['integration_backend']=dict(name='existing_Source_span_current_neighborhood_Lee',automatic_model_pipeline_codeobject_exact=True,original_module_globals_modified=False);return r
