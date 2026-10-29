"""Reject section components with detected foreign source-branch ownership."""
import ast
import copy
from functools import lru_cache
import hashlib
import math
from pathlib import Path
import types
import numpy as np
from scipy.spatial import cKDTree
from root_graph_measurement_ledger import graph_measurement_ledger
import root_section_centroid as original


class BranchOwnership:
    def __init__(self,points,adjacency):
        self.points=np.asarray(points,float)
        self.ledger=graph_measurement_ledger(points,adjacency,coordinate_units='source_graph_units')
        edges=[];owners=[];self.node_branches=[set() for unused in self.points]
        for branch in self.ledger['branches']:
            owner=branch['branch_id']
            for node in branch['source_node_ids']:self.node_branches[node].add(owner)
            for edge in branch['source_edge_ids']:edges.append(edge);owners.append(owner)
        e=np.array(edges,int).reshape(-1,2);self.owners=np.array(owners,int)
        self.a=self.points[e[:,0]];self.b=self.points[e[:,1]];self.d=self.b-self.a
        self.squared=np.einsum('ij,ij->i',self.d,self.d)
        mid=(self.a+self.b)/2;self.tree=cKDTree(mid) if len(mid) else None
        self.max_half=max((math.dist(a,b)/2 for a,b in zip(self.a,self.b)),default=0.)

    def _distances(self,q,ids):
        ids=np.asarray(ids,int);v=q-self.a[ids];denom=self.squared[ids]
        t=np.zeros(len(ids));positive=denom>0
        t[positive]=np.clip(np.einsum('ij,ij->i',v[positive],self.d[ids][positive])/denom[positive],0,1)
        delta=v-t[:,None]*self.d[ids]
        return np.array([math.hypot(*row) for row in delta])

    def nearest_owners(self,q):
        q=np.asarray(q,float)
        if q.shape!=(3,) or not np.isfinite(q).all():raise ValueError('FINITE_OWNERSHIP_QUERY')
        if self.tree is None:return set()
        initial=int(self.tree.query(q)[1]);upper=float(self._distances(q,[initial])[0])
        pad=64*np.finfo(float).eps*max(1.,math.hypot(*q),upper,self.max_half)
        ids=self.tree.query_ball_point(q,upper+self.max_half+pad)
        distances=self._distances(q,ids);minimum=float(distances.min())
        return set(map(int,self.owners[np.asarray(ids)[distances<=minimum+pad]]))


def _witnesses(polygons,ids,point,normal,pitch):
    points=[]
    for i in sorted(ids):
        vertices=[]
        for key in polygons[i]:
            v=np.array(key[1:4],float)
            if key[0]:v[key[4]]+=np.dot(normal,point-v)/normal[key[4]]
            vertices.append(v)
        points.append(np.mean(vertices,axis=0)*pitch)
    return points


def _receipt(polygons,ids,point,normal,basis,pitch):
    return dict(original._moments(polygons,ids,point,normal,basis,pitch),
        _interior_polygon_witnesses=[q.tolist() for q in _witnesses(polygons,ids,point,normal,pitch)])


@lru_cache(maxsize=1)
def owned_section_class():
    path=Path(__file__).with_name('root_cubical_sections.py');raw=path.read_bytes()
    assert hashlib.sha256(raw).hexdigest()==original.KERNEL_SHA
    tree=ast.parse(raw);reference=ast.dump(tree,include_attributes=False);matched=[]
    for cls in tree.body:
        if isinstance(cls,ast.ClassDef) and cls.name=='CubicalSectionMeasurer':
            for fun in cls.body:
                if isinstance(fun,ast.FunctionDef) and fun.name=='section':
                    for index,node in enumerate(fun.body):
                        if isinstance(node,ast.Return) and 'CLOSED_DISK_CUBICAL_SECTION' in ast.dump(node):
                            fun.body.insert(index,ast.parse("result['section_centroid_evidence']=_receipt(polygons,ids,point,normal,basis,self.solid.pitch)").body[0])
                            matched.append(index);break
    assert len(matched)==1
    restored=copy.deepcopy(tree)
    for cls in restored.body:
        if isinstance(cls,ast.ClassDef) and cls.name=='CubicalSectionMeasurer':
            for fun in cls.body:
                if isinstance(fun,ast.FunctionDef) and fun.name=='section':del fun.body[matched[0]]
    assert ast.dump(restored,include_attributes=False)==reference
    namespace={'_receipt':_receipt};exec(compile(ast.fix_missing_locations(tree),str(path),'exec'),namespace)
    return namespace['CubicalSectionMeasurer']


def branch_section_centroid_measurement(solid,points,adjacency,*,progress=None):
    ownership=BranchOwnership(points,adjacency);events=[];p=np.asarray(points,float)
    class Engine:
        def __init__(self,source):self.kernel=owned_section_class()(source)
        def section(self,point,normal):
            result=self.kernel.section(point,normal)
            evidence=result.get('section_centroid_evidence')
            if evidence is None:return result
            witnesses=evidence.pop('_interior_polygon_witnesses')
            matches=np.flatnonzero(np.all(p==point,axis=1))
            if len(matches)!=1:
                result.pop('section_centroid_evidence');return dict(result,status='MISSING_UNIQUE_SOURCE_QUERY_ID',radius=None)
            node=int(matches[0]);own=ownership.node_branches[node]
            foreign=[]
            for q in witnesses:
                nearest=ownership.nearest_owners(q)
                if not nearest or nearest!=own:foreign.append(dict(point=q,nearest_branches=sorted(nearest)))
            events.append(dict(node_id=node,source_branch_ids=sorted(own),witnesses=len(witnesses),foreign_witnesses=foreign))
            if foreign:
                result.pop('section_centroid_evidence')
                return dict(result,status='MISSING_MIXED_BRANCH_SECTION',radius=None)
            return result
    namespace=dict(vars(original));namespace['centroid_section_class']=lambda:Engine
    method=types.FunctionType(original.section_centroid_measurement.__code__,namespace,
        name='section_centroid_measurement',argdefs=original.section_centroid_measurement.__defaults__)
    method.__kwdefaults__=original.section_centroid_measurement.__kwdefaults__
    assert method.__code__ is original.section_centroid_measurement.__code__
    result=method(solid,points,adjacency,progress=progress)
    result['branch_section_association']=dict(events=events,queries=len(events),
        rejected_queries=sum(bool(e['foreign_witnesses']) for e in events),
        exact_original_centroid_function_codeobject=True,
        no_detected_foreign_polygon_witness_not_entire_polygon_uniqueness_proof=True,
        missing_radius_coverage_must_be_reported=True)
    result['status']='COMPLETE_EXPERIMENTAL_BRANCH_ASSOCIATED_SECTION_GRAPH'
    return result
