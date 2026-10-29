"""Experimental source-derived span anchors, not physical root-tip inference.

EDT is distance to background sample centers, not exact cubical-boundary radius.
This alternative remains separate from the qualified nonempty-only guard.
"""
import ast
import copy
import hashlib
import inspect
import numpy as np
from scipy import ndimage
import root_guarded_thinning as original


def plateau_span_anchors(volume):
    source=original._input(volume,350000000)
    if not source.any():return np.zeros_like(source),dict(anchor_count=0)
    distance=ndimage.distance_transform_edt(np.pad(source,1))[1:-1,1:-1,1:-1]
    plateau=source & (distance==float(distance.max()))
    points=np.argwhere(plateau)
    center=points.mean(axis=0);centered=points-center
    values,vectors=np.linalg.eigh(centered.T@centered)
    tied=int(np.count_nonzero(values==values[-1]))
    if tied>1:
        direction=np.eye(3)[int(np.argmax(np.ptp(points,axis=0)))]
    else:direction=vectors[:,-1]
    if direction[int(np.argmax(np.abs(direction)))]<0:direction=-direction
    scores=centered@direction
    chosen=[]
    for extreme in [float(scores.min()),float(scores.max())]:
        candidates=points[scores==extreme]
        distances=np.sum((candidates-center)**2,axis=1)
        candidates=candidates[distances==distances.min()]
        chosen.append(tuple(map(int,sorted(map(tuple,candidates))[0])))
    anchor_mask=np.zeros_like(source)
    for point in chosen:anchor_mask[point]=True
    assert not np.any(anchor_mask&~source)
    return anchor_mask,dict(anchor_count=int(anchor_mask.sum()),anchor_indices=[list(p) for p in sorted(set(chosen))],
        maximal_EDT_plateau_samples=len(points),leading_PCA_eigenvalue=float(values[-1]),
        exact_repeated_leading_eigenvalue=tied>1,anchor_axis=direction.tolist(),
        source_EDT_definition='distance to background sample centers on padded source grid',
        physical_endpoints_inferred=False,geometric_centerline_qualified=False)


def span_anchored_component(volume,max_dense_cells=350000000):
    source=original._input(volume,max_dense_cells)
    anchors,receipt=plateau_span_anchors(source)
    tree=ast.parse(inspect.getsource(original.sequential_lee_thinning));before=ast.dump(tree,include_attributes=False)
    count=[0]
    class Bind(ast.NodeTransformer):
        def visit_If(self,node):
            self.generic_visit(node)
            if isinstance(node.test,ast.Call) and isinstance(node.test.func,ast.Name) and node.test.func.id=='deletable_current':
                guard=ast.parse('not anchor_mask[a-1,b-1,c-1]',mode='eval').body
                node.test=ast.BoolOp(op=ast.And(),values=[node.test,guard]);count[0]+=1
            return node
    bound=Bind().visit(tree);ast.fix_missing_locations(bound);assert count[0]==2
    class Restore(ast.NodeTransformer):
        def visit_If(self,node):
            if (isinstance(node.test,ast.BoolOp) and isinstance(node.test.op,ast.And)
                and isinstance(node.test.values[0],ast.Call) and isinstance(node.test.values[0].func,ast.Name)
                and node.test.values[0].func.id=='deletable_current'):
                node.test=node.test.values[0]
            return self.generic_visit(node)
    assert ast.dump(Restore().visit(copy.deepcopy(bound)),include_attributes=False)==before
    namespace=dict(vars(original),anchor_mask=anchors)
    exec(compile(bound,'experimental_source_span_anchor_predicate_binding','exec'),namespace)
    result,audit=namespace['sequential_lee_thinning'](source,max_dense_cells=max_dense_cells)
    assert np.all(result[anchors])
    receipt.update(original_sequential_AST_restored_exact=True,anchor_predicate_insertions=2,
        original_nonempty_audit=audit,method_adopted=False)
    return result,receipt


def experimental_span_skeletonize(volume,method='lee',*,max_dense_cells=350000000,receipt=None):
    if method!='lee':raise ValueError('ONLY_LEE_EXPERIMENT_SUPPORTED')
    source=original._input(volume,max_dense_cells);source_hash=hashlib.sha256(source.tobytes()).hexdigest()
    result=original.compiled_lee(source,method='lee')
    labels,count=ndimage.label(source,structure=original.CONNECT26)
    present=set(map(int,np.unique(labels[result])))-{0}
    missing=sorted(set(range(1,count+1))-present);repairs=[]
    if missing:
        result=result.copy();boxes=ndimage.find_objects(labels)
        for component in missing:
            box=boxes[component-1];part=labels[box]==component
            recovered,audit=span_anchored_component(part,max_dense_cells=max_dense_cells)
            result[box]|=recovered
            repairs.append(dict(source_component=component,source_voxels=int(part.sum()),
                recovered_derived_voxels=int(recovered.sum()),audit=audit))
    assert not np.any(result&~source) and hashlib.sha256(source.tobytes()).hexdigest()==source_hash
    after=ndimage.label(result,structure=original.CONNECT26)[1];assert count==after
    audit=dict(status='EXPERIMENTAL_SPAN_ANCHORED_FALLBACK_NOT_ADOPTED',source_components=int(count),
        derived_components=int(after),missing_compiled_components=missing,fallback_components=repairs,
        original_compiled_output_unchanged=not missing,source_material_unchanged=True,
        occupied_material_added=0,components_joined=0,genuine_cycles_forced_to_tree=False,
        universal_topological_or_biological_qualification=False,physical_endpoints_inferred=False,
        method_adopted=False)
    if receipt is not None:receipt.update(audit)
    return result
