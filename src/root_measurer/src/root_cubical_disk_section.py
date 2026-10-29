"""Resolve internal cut-cell query boundaries without weakening disk topology."""
import ast
import copy
from functools import lru_cache
import hashlib
from pathlib import Path
import numpy as np
from root_section_centroid import KERNEL_SHA
from root_branch_section_measurement import _receipt


def _on_outer_boundary(boundary_graph,point,normal,tol):
    def vertex(key):
        p=np.array(key[1:4],float)
        if key[0]:p[key[4]]+=np.dot(normal,point-p)/normal[key[4]]
        return p
    for a,neighbors in boundary_graph.items():
        start=vertex(a)
        for b in neighbors:
            if a>=b:continue
            end=vertex(b);d=end-start;denominator=float(d@d)
            t=np.clip(float((point-start)@d)/denominator,0,1) if denominator else 0
            if np.linalg.norm(point-start-t*d)<=tol:return True
    return False


@lru_cache(maxsize=1)
def interior_boundary_section_class():
    path=Path(__file__).with_name('root_cubical_sections.py');raw=path.read_bytes()
    assert hashlib.sha256(raw).hexdigest()==KERNEL_SHA
    tree=ast.parse(raw);reference=ast.dump(tree,include_attributes=False)
    fun=next(f for cls in tree.body if isinstance(cls,ast.ClassDef) and cls.name=='CubicalSectionMeasurer'
             for f in cls.body if isinstance(f,ast.FunctionDef) and f.name=='section')
    removed=[(i,n) for i,n in enumerate(fun.body) if isinstance(n,ast.If) and 'QUERY_ON_CELL_SECTION_BOUNDARY' in ast.dump(n)]
    assert len(removed)==1;old_index,old_gate=removed[0];del fun.body[old_index]
    additions=[]
    at=next(i for i,n in enumerate(fun.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='selected' for t in n.targets))
    gate=ast.parse("if len({find(int(i)) for i in hits}) != 1:\n    return dict(base,status='QUERY_ON_AMBIGUOUS_SECTION_COMPONENT_CONTACT',candidate_cells=len(candidates),polygons=len(polygons))").body[0]
    fun.body.insert(at,gate);additions.append(ast.dump(gate,include_attributes=False))
    at=next(i for i,n in enumerate(fun.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='visited' for t in n.targets))
    gate=ast.parse("if _on_outer_boundary(boundary_graph,point,normal,tol):\n    return dict(result,status='QUERY_ON_SECTION_OUTER_BOUNDARY')").body[0]
    fun.body.insert(at,gate);additions.append(ast.dump(gate,include_attributes=False))
    at=next(i for i,n in enumerate(fun.body) if isinstance(n,ast.Return) and 'CLOSED_DISK_CUBICAL_SECTION' in ast.dump(n))
    nodes=ast.parse("result['section_centroid_evidence']=_receipt(polygons,ids,point,normal,basis,self.solid.pitch)\nresult['query_boundary_resolution']=dict(hit_cut_polygons=len(hits),internal_partition_boundary_resolved=bool(len(hits)>1 or any(on_boundary[int(i)] for i in hits)),outer_union_boundary_excluded=True)").body
    for offset,node in enumerate(nodes):fun.body.insert(at+offset,node);additions.append(ast.dump(node,include_attributes=False))
    restored=copy.deepcopy(tree)
    rf=next(f for cls in restored.body if isinstance(cls,ast.ClassDef) and cls.name=='CubicalSectionMeasurer'
            for f in cls.body if isinstance(f,ast.FunctionDef) and f.name=='section')
    rf.body=[n for n in rf.body if ast.dump(n,include_attributes=False) not in additions]
    rf.body.insert(old_index,copy.deepcopy(old_gate))
    assert ast.dump(restored,include_attributes=False)==reference
    namespace={'_on_outer_boundary':_on_outer_boundary,'_receipt':_receipt}
    exec(compile(ast.fix_missing_locations(tree),str(path),'exec'),namespace)
    return namespace['CubicalSectionMeasurer']
