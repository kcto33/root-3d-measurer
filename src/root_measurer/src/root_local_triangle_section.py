"""Query-containing isolated triangle section disk; no overlapping-area union."""
import ast
import copy
from fractions import Fraction
from functools import lru_cache
import hashlib
from pathlib import Path
from root_triangle_plane_section import _intersects,_orient


def _contains(poly,q):
    inside=False
    for a,b in zip(poly,poly[1:]+poly[:1]):
        if _orient(a,b,q)==0 and all(min(a[i],b[i])<=q[i]<=max(a[i],b[i]) for i in range(2)):return True
        if (a[1]>q[1])!=(b[1]>q[1]):
            x=a[0]+(q[1]-a[1])*(b[0]-a[0])/(b[1]-a[1])
            if x>q[0]:inside=not inside
    return inside


def _selection(self,loops,p,n):
    axis=max(range(3),key=lambda i:abs(n[i]));other=[i for i in range(3) if i!=axis]
    polygons=[[tuple(v[k] for k in other) for v in path] for path in loops];q=tuple(p[k] for k in other)
    receipt=dict(all_contours=len(loops),unselected_contours_retained=True,overlap_or_holes_filled=False,
        selection_rule='unique simple disk containing actual query; all contours mutually disjoint and nonnested')
    self.local_selection_receipt=receipt
    for poly in polygons:
        for i in range(len(poly)):
            for j in range(i+1,len(poly)):
                if j==i+1 or (i==0 and j==len(poly)-1):continue
                if _intersects(poly[i],poly[(i+1)%len(poly)],poly[j],poly[(j+1)%len(poly)]):
                    return dict(path=None,status='MISSING_SELF_INTERSECTING_LOCAL_CONTOUR')
    for i,a in enumerate(polygons):
        for b in polygons[i+1:]:
            if any(_intersects(u,v,x,y) for u,v in zip(a,a[1:]+a[:1]) for x,y in zip(b,b[1:]+b[:1])):
                return dict(path=None,status='MISSING_CROSSING_OR_TOUCHING_SECTION_CONTOURS')
            if _contains(a,b[0]) or _contains(b,a[0]):
                return dict(path=None,status='MISSING_NESTED_SECTION_CONTOURS_NOT_UNIONED')
    selected=[i for i,poly in enumerate(polygons) if _contains(poly,q)]
    receipt.update(query_containing_contours=selected,unselected_contours=len(loops)-len(selected))
    if len(selected)!=1:return dict(path=None,status='MISSING_UNIQUE_QUERY_CONTAINING_DISK')
    receipt['selected_contour_id']=selected[0]
    return dict(path=loops[selected[0]],status='COMPLETE_LOCAL_QUERY_CONTAINING_DISK')


@lru_cache(maxsize=1)
def local_triangle_section_class():
    path=Path(__file__).with_name('root_triangle_plane_section.py');raw=path.read_bytes();tree=ast.parse(raw)
    original=ast.dump(tree,include_attributes=False)
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='TrianglePlaneSectionMeasurer')
    section=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='section')
    i=next(i for i,n in enumerate(section.body) if isinstance(n,ast.If) and ast.unparse(n.test)=='len(loops) != 1')
    old_if=copy.deepcopy(section.body[i]);old_next=copy.deepcopy(section.body[i+1])
    assert ast.unparse(old_next).startswith('path = loops[0]')
    replacement=ast.parse("selection=self._select_local_disk(loops,p,n)\nif selection['path'] is None:return dict(base,status=selection['status'],contours=len(loops))\npath=selection['path']").body
    # The original next statement is a semicolon-expanded assignment followed
    # by an independent axis assignment. Only path=loops[0] is replaced.
    section.body[i:i+2]=replacement
    restored=copy.deepcopy(tree)
    rc=next(n for n in restored.body if isinstance(n,ast.ClassDef) and n.name==cls.name)
    rf=next(n for n in rc.body if isinstance(n,ast.FunctionDef) and n.name=='section')
    rf.body[i:i+len(replacement)]=[old_if,old_next]
    assert ast.dump(restored,include_attributes=False)==original
    namespace={};exec(compile(ast.fix_missing_locations(tree),str(path),'exec'),namespace)
    kind=namespace['TrianglePlaneSectionMeasurer'];kind._select_local_disk=_selection
    kind.original_source_sha256=hashlib.sha256(raw).hexdigest();kind.inverse_AST_exact=True
    return kind


class LocalTriangleSectionMeasurer:
    def __init__(self,mesh):self.kernel=local_triangle_section_class()(mesh)
    def section(self,point,normal):
        self.kernel.local_selection_receipt=None
        result=self.kernel.section(point,normal)
        result['local_contour_selection']=self.kernel.local_selection_receipt
        result['whole_surface_material_union_admitted']=False
        result['original_triangle_section_source_SHA']=self.kernel.original_source_sha256
        result['inverse_AST_original_kernel_exact']=True
        if result['radius'] is not None:
            result['radius_definition']='sqrt(area/pi) of the isolated simple section disk containing the actual query'
        return result
