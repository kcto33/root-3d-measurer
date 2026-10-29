"""Source-only oriented-cut reference with all disconnected components retained."""
import ast
import copy
from collections import defaultdict
from functools import lru_cache
import hashlib
from pathlib import Path
import root_star_section_reference as original


def _components(segments):
    incident=defaultdict(list)
    for i,(a,b) in enumerate(segments):incident[a].append(i);incident[b].append(i)
    if any(len(rows)!=2 for rows in incident.values()):return None
    unseen=set(range(len(segments)));parts=[]
    while unseen:
        start=min(unseen);unseen.remove(start);stack=[start];ids=[]
        while stack:
            i=stack.pop();ids.append(i)
            for point in segments[i]:
                for j in incident[point]:
                    if j in unseen:unseen.remove(j);stack.append(j)
        parts.append([segments[i] for i in sorted(ids)])
    assert sum(map(len,parts))==len(segments)
    return parts


@lru_cache(maxsize=1)
def _reference_kind():
    path=Path(original.__file__);raw=path.read_bytes();tree=ast.parse(raw);before=ast.dump(tree,include_attributes=False)
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='StarSectionReference')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='section')
    chain=next(n for n in fn.body if isinstance(n,ast.For) and ast.unparse(n.target)=='(group, (begin, end))')
    at=next(i for i,n in enumerate(chain.body) if isinstance(n,ast.Assign) and any(isinstance(x,ast.Name) and x.id=='balance' for x in n.targets))
    old=copy.deepcopy(chain.body[at:]);new=copy.deepcopy(old)
    class Group(ast.NodeTransformer):
        def visit_Name(self,node):
            if node.id=='segments':return ast.copy_location(ast.Name('component_segments',node.ctx),node)
            return node
        def visit_keyword(self,node):
            self.generic_visit(node)
            if node.arg=='chain_id':node.value=ast.List(elts=[ast.Name('group',ast.Load()),ast.Name('component_id',ast.Load())],ctx=ast.Load())
            return node
    new=[Group().visit(n) for n in new]
    addition=ast.parse("components=_components(segments)\nif components is None:return dict(base,status='REFERENCE_MISSING_NONMANIFOLD_CUT_COMPONENTS')\nfor component_id,component_segments in enumerate(components):\n    pass").body
    addition[-1].body=new;chain.body[at:]=addition
    restored=copy.deepcopy(tree)
    rc=next(n for n in restored.body if isinstance(n,ast.ClassDef) and n.name==cls.name)
    rf=next(n for n in rc.body if isinstance(n,ast.FunctionDef) and n.name==fn.name)
    rl=next(n for n in rf.body if isinstance(n,ast.For) and ast.unparse(n.target)=='(group, (begin, end))')
    rl.body[at:]=old
    assert ast.dump(restored,include_attributes=False)==before
    namespace=dict(vars(original),_components=_components)
    exec(compile(ast.fix_missing_locations(tree),str(path),'exec'),namespace)
    return namespace[cls.name],hashlib.sha256(raw).hexdigest()


class DisconnectedCutSectionReference:
    def __init__(self,mesh,face_ranges=None):
        kind,self.original_source_sha256=_reference_kind();self.engine=kind(mesh,face_ranges);self.mesh=mesh

    def section(self,point,normal):
        result=self.engine.section(point,normal)
        result.update(reference_grouping='all exact endpoint-connected oriented-cut components within each source chain',
            all_oriented_cut_components_retained=True,ordered_contour_path_not_reconstructed=True,
            original_star_reference_source_sha256=self.original_source_sha256,
            grouping_only_inverse_AST_exact=True,measured_radius_or_GT_not_accepted=True)
        return result
