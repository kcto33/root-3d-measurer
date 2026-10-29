"""Fresh whole-source sections at the exact requested rational point/normal.

The former binary64 query provider remains unchanged. Only query encoding and
complete-plane indexing differ here; strict local disk refusals are retained.
"""
import ast,builtins,copy,hashlib,math,os
from fractions import Fraction as F
from functools import lru_cache
from types import FunctionType,MethodType
import numpy as np
import root_local_triangle_section as local
from root_triangle_measurement_input import TriangleInput
from root_fragment_triangle_section import _boundary
from root_graph_supported_length import _q
from root_graph_source_length import q,pack

def _bounded(values):
    if any(max(abs(x.numerator).bit_length(),x.denominator.bit_length())>131072 for x in values):
        raise ValueError('UNCHANGED_INTERNAL_131072_BIT_GUARD')
    return values

def _exact_query(point,normal):
    if len(point)!=3 or len(normal)!=3:raise ValueError('EXACT_THREE_COORDINATES_REQUIRED')
    def decode(x):return _q(q(x)) if isinstance(x,(tuple,list)) else _q(x)
    p,n=tuple(map(decode,point)),tuple(map(decode,normal))
    if not any(n):raise ValueError('EXACT_NONZERO_PLANE_NORMAL_REQUIRED')
    return p,n

@lru_cache(maxsize=1)
def _exact_indexed_kind():
    evidence={}
    def compile_index(tree,filename,mode):
        tree=copy.deepcopy(tree);before=ast.dump(tree,include_attributes=False)
        cls=next(x for x in tree.body if isinstance(x,ast.ClassDef) and x.name=='TrianglePlaneSectionMeasurer')
        section=next(x for x in cls.body if isinstance(x,ast.FunctionDef) and x.name=='section')
        prefix=copy.deepcopy(section.body[:5])
        assert ast.unparse(prefix[0]).startswith('p = np.asarray(point, float)') and ast.unparse(prefix[4]).startswith('n = tuple(')
        section.body[:5]=ast.parse('p,n=_exact_query(point,normal)').body
        distance=next(x for x in section.body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='distance' for t in x.targets))
        loop=next(x for x in section.body if isinstance(x,ast.For) and ast.unparse(x.target)=='(face_id, face)')
        od,oi=copy.deepcopy(distance.value),copy.deepcopy(loop.iter)
        distance.value=ast.parse('self._exact_indexed_distances',mode='eval').body
        loop.iter=ast.parse('zip(self._exact_indexed_face_ids,self.faces)',mode='eval').body
        imported=ast.parse('from root_exact_query_sections import _exact_query').body[0]
        tree.body.insert(0,imported);executable=copy.deepcopy(tree)
        tree.body.pop(0);distance.value=od;loop.iter=oi;section.body[:1]=prefix
        assert ast.dump(tree,include_attributes=False)==before
        evidence.update(query_prefix_and_distance_ID_iteration_inverse_AST_exact=True,old_local_contour_rules_unchanged=True)
        return builtins.compile(ast.fix_missing_locations(executable),filename,mode)
    factory=local.local_triangle_section_class.__wrapped__
    bound=FunctionType(factory.__code__,dict(factory.__globals__,compile=compile_index),factory.__name__,factory.__defaults__)
    kind=bound();kind.exact_query_binding=evidence
    return kind

class ExactQuerySectionMeasurer:
    def __init__(self,mesh,*,execution_domain='local',max_queries=5,max_candidate_faces=100000,max_work=10000000):
        if not isinstance(mesh,TriangleInput):raise ValueError('QUALIFIED_ORIGINAL_TRIANGLE_INPUT_REQUIRED')
        nv,nf=len(mesh.vertices),len(mesh.faces)
        if not 0<nv<=1200000 or not 0<nf<=2000000:raise ValueError('COMPLETE_ORIGINAL_SOURCE_SIZE_GUARD')
        if execution_domain=='local':
            if nv>3000 or nf>1000:raise ValueError('UNCHANGED_LOCAL_SOURCE_GUARDS')
        elif execution_domain!='registered_slurm_complete_source' or not os.environ.get('SLURM_JOB_ID','').isdigit():raise ValueError('REGISTERED_SLURM_BEFORE_COMPLETE_SOURCE_INDEX')
        if type(max_queries)is not int or not 1<=max_queries<=5 or type(max_candidate_faces)is not int or not 1<=max_candidate_faces<=100000 or type(max_work)is not int or not 0<=max_work<=10000000:
            raise ValueError('QUERY_CANDIDATE_WORK_GUARDS_ONLY_TIGHTENED')
        self.mesh=mesh;self.max_queries=max_queries;self.max_candidates=max_candidate_faces;self.max_work=max_work;self.queries=0
        self.before=hashlib.sha256(mesh.vertices.tobytes()+mesh.faces.tobytes()).hexdigest()
        self.power=max(float(x).as_integer_ratio()[1].bit_length()-1 for x in mesh.vertices.flat)
        self.integer_vertices=[]
        for row in mesh.vertices:
            self.integer_vertices.append(tuple(a<<(self.power-(d.bit_length()-1)) for a,d in (float(x).as_integer_ratio() for x in row)))
        kind=_exact_indexed_kind();self.kernel=kind(mesh);self.selected_path=None;self.pair_bound=0
        def select(instance,loops,p,n):
            self.selected_path=None;total=sum(map(len,loops));self.pair_bound=total*(total-1)//2
            for loop in loops:
                for point in loop:_bounded(point)
            if any(len(loop)>256 for loop in loops):return dict(path=None,status='MISSING_UNCHANGED_256_CYCLE_GUARD')
            # Local selection tests all contour pairs; the retained original
            # section kernel checks the selected contour once again.
            if self.index_work+2*self.pair_bound>self.max_work:return dict(path=None,status='MISSING_EXACT_SECTION_WORK_GUARD')
            result=local._selection(instance,loops,p,n);self.selected_path=result['path'];return result
        self.kernel._select_local_disk=MethodType(select,self.kernel)

    def _unchanged(self):
        if hashlib.sha256(self.mesh.vertices.tobytes()+self.mesh.faces.tobytes()).hexdigest()!=self.before:raise ValueError('ORIGINAL_SOURCE_ARRAYS_CHANGED')

    def section(self,point,normal):
        p,n=_exact_query(point,normal);self._unchanged()
        if self.queries>=self.max_queries:raise ValueError('REGISTERED_QUERY_COUNT_EXHAUSTED')
        self.queries+=1;self.selected_path=None;self.pair_bound=0;nv,nf=len(self.mesh.vertices),len(self.mesh.faces)
        base=dict(area=None,radius=None,centroid=None,source_sha256=self.mesh.source_sha256,source_arrays_sha256=self.before,
            exact_query_point=list(map(pack,p)),exact_query_normal=list(map(pack,n)),query_was_rounded_to_binary64=False,
            old_binary64_query_provider_modified=False,Source_geometry_or_IDs_modified=False,Source_region_or_physical_radius_admitted=False)
        if 3*nv+nf>self.max_work:return dict(base,status='MISSING_EXACT_SECTION_WORK_GUARD')
        nd=math.lcm(*(x.denominator for x in n));pd=math.lcm(*(x.denominator for x in p));coeff=[int(x*nd) for x in n]
        origin=[int(x*pd) for x in p];scale=1<<self.power;offset=sum(a*b for a,b in zip(coeff,origin))*scale
        denominator=nd*pd*scale
        if max(abs(offset).bit_length(),denominator.bit_length(),*(abs(x).bit_length() for x in coeff))>131072:raise ValueError('UNCHANGED_INTERNAL_131072_BIT_GUARD')
        def distance(v):return sum(a*x*pd for a,x in zip(coeff,v))-offset
        signs=np.fromiter(((d>0)-(d<0) for d in map(distance,self.integer_vertices)),dtype=np.int8,count=nv)
        fs=signs[self.mesh.faces];mask=(fs.min(axis=1)<=0)&(fs.max(axis=1)>=0);ids=np.flatnonzero(mask)
        self.index_work=3*nv+nf+3*len(ids)
        evidence=dict(source_vertices=nv,source_faces=nf,candidate_faces=len(ids),strictly_one_side_excluded_faces=int((~mask).sum()),
            original_candidate_face_ids_first1000=ids[:1000].tolist(),candidate_ID_details_complete=len(ids)<=1000,
            complete_candidate_face_ID_array_SHA256=hashlib.sha256(ids.astype('<i8').tobytes()).hexdigest(),
            plane_sign_array_sha256=hashlib.sha256(signs.tobytes()).hexdigest(),all_original_faces_classified=True,
            all_coplanar_tangent_and_crossing_faces_retained=True,exact_integer_rational_plane_signs=True,query_index=self.queries,
            counted_index_work=self.index_work,binding=_exact_indexed_kind().exact_query_binding.copy())
        def missing(status):
            self._unchanged();return dict(base,status=status,index_evidence=evidence)
        if len(ids)>self.max_candidates:return missing('MISSING_EXACT_QUERY_FACE_GUARD')
        if self.index_work>self.max_work:return missing('MISSING_EXACT_SECTION_WORK_GUARD')
        distances=[None]*nv
        for i in np.unique(self.mesh.faces[ids]):distances[int(i)]=F(distance(self.integer_vertices[int(i)]),denominator)
        k=self.kernel;k._exact_indexed_distances=distances;k._exact_indexed_face_ids=list(map(int,ids));k.faces=self.mesh.faces[ids];k.local_selection_receipt=None
        result=k.section(p,n);self._unchanged()
        result.update(base,area=result.get('area'),radius=result.get('radius'),index_evidence=evidence,
            geometry='exact rational query plane and original binary64 triangle geometry; only final displays/area-radius square roots rounded',
            local_contour_selection=k.local_selection_receipt,centroid_missing_reason=result['status'],
            contour_pair_comparison_upper_bound=2*self.pair_bound,work_count_is_not_CPU_instruction_count=True)
        if result['radius'] is None:return result
        _bounded(tuple(F(int(v['numerator']),int(v['denominator'])) for v in result['exact_area_vector']))
        path=self.selected_path;assert path is not None
        drop=max(range(3),key=lambda i:abs(n[i]));keep=[i for i in range(3) if i!=drop];xy=[tuple(v[i] for i in keep) for v in path]
        weights=[];S=F(0);T=[F(0),F(0)]
        for a,b in zip(xy,xy[1:]+xy[:1]):
            w=_bounded((a[0]*b[1]-a[1]*b[0],))[0];weights.append(w);S=_bounded((S+w,))[0]
            for i in range(2):T[i]=_bounded((T[i]+(a[i]+b[i])*w,))[0]
        assert S
        c=_bounded(tuple(T[i]/(3*S) for i in range(2)))
        inside=local._contains(xy,c) and not _boundary(xy,c);result['centroid_strictly_inside_selected_disk']=inside
        if not inside:result['centroid_missing_reason']='AREA_CENTROID_OUTSIDE_OR_ON_NONCONVEX_DISK';return result
        center=[F(0)]*3
        for i,value in zip(keep,c):center[i]=value
        center[drop]=p[drop]-sum(n[i]*(center[i]-p[i]) for i in keep)/n[drop]
        _bounded(center)
        assert sum(n[i]*(center[i]-p[i]) for i in range(3))==0
        result.update(centroid=list(map(float,center)),exact_centroid=list(map(pack,center)),centroid_missing_reason=None,
            centroid_definition='exact polygon moments on this exact rational query plane; no graph movement',
            final_area_and_radius_sqrt_rounded=True,fresh_complete_original_Source_plane_query=True)
        return result
