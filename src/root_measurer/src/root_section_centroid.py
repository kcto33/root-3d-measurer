"""One-pass, endpoint-fixed section-centroid geometry estimator; no GT input."""
import ast
from collections import Counter
import copy
from functools import lru_cache
import hashlib
import math
from pathlib import Path
import numpy as np
from root_graph_measurement_ledger import graph_measurement_ledger

KERNEL_SHA='399def4141b43d5d46e08373606c8d5bcfdf3da4ec8b6c341998265884f607ce'


def _moments(polygons,ids,point,normal,basis,pitch):
    def vertex(key):
        value=np.array(key[1:4],float)
        if key[0]:
            axis=key[4]
            value[axis]+=np.dot(normal,point-value)/normal[axis]
        return value
    areas=[];moments=[]
    for i in sorted(ids):
        xy=np.array([(vertex(key)-point)@basis for key in polygons[i]])
        nxt=np.roll(xy,-1,axis=0);cross=xy[:,0]*nxt[:,1]-nxt[:,0]*xy[:,1]
        area=math.fsum(map(float,cross))/2
        assert area>0
        # Integral x dA = sum((x_i+x_j)*cross_ij)/6.
        moment=np.array([math.fsum(map(float,(xy[:,a]+nxt[:,a])*cross))/6 for a in range(2)])
        areas.append(area);moments.append(moment)
    total=math.fsum(areas)
    center=np.array([math.fsum(float(m[a]) for m in moments)/total for a in range(2)])
    return dict(centroid=((point+basis@center)*pitch).tolist(),area_moment_sum=total*pitch**2,
        selected_cut_polygons=len(ids),definition='area centroid of unchanged selected closed-disk plane section')


@lru_cache(maxsize=1)
def centroid_section_class():
    path=Path(__file__).with_name('root_cubical_sections.py');raw=path.read_bytes()
    assert hashlib.sha256(raw).hexdigest()==KERNEL_SHA
    original=ast.parse(raw);tree=copy.deepcopy(original);matches=[]
    for cls in tree.body:
        if isinstance(cls,ast.ClassDef) and cls.name=='CubicalSectionMeasurer':
            for fun in cls.body:
                if isinstance(fun,ast.FunctionDef) and fun.name=='section':
                    for index,node in enumerate(fun.body):
                        if isinstance(node,ast.Return) and 'CLOSED_DISK_CUBICAL_SECTION' in ast.dump(node):
                            fun.body.insert(index,ast.parse("result['section_centroid_evidence']=_moments(polygons,ids,point,normal,basis,self.solid.pitch)").body[0])
                            matches.append((cls.name,fun.name,index));break
    assert len(matches)==1
    restored=copy.deepcopy(tree)
    for cls in restored.body:
        if isinstance(cls,ast.ClassDef) and cls.name==matches[0][0]:
            for fun in cls.body:
                if isinstance(fun,ast.FunctionDef) and fun.name==matches[0][1]:del fun.body[matches[0][2]]
    assert ast.dump(restored,include_attributes=False)==ast.dump(original,include_attributes=False)
    namespace={'_moments':_moments}
    exec(compile(ast.fix_missing_locations(tree),str(path),'exec'),namespace)
    return namespace['CubicalSectionMeasurer']


def _normal(path,at,pitch,closed):
    lengths=np.array([math.dist(a,b) for a,b in zip(path,path[1:])]);arc=np.r_[0,np.cumsum(lengths)]
    if np.any(lengths<=0) or not len(lengths):return None
    position=arc[at];extent=arc[-1]
    if closed:
        # Use at most one traversal; no repeated-loop fitting.
        width=min(10*pitch,extent)
        locations=position+np.linspace(-width/2,width/2,21)
        locations=np.mod(locations,extent)
    else:
        low=max(0.,position-5*pitch);high=min(extent,position+5*pitch)
        if at==0:high=min(extent,10*pitch)
        if at==len(path)-1:low=max(0.,extent-10*pitch)
        locations=np.linspace(low,high,21)
    q=np.column_stack([np.interp(locations,arc,path[:,axis]) for axis in range(3)])-path[at]
    derivative=np.array([math.fsum(2*(i-10)/77*float(q[i,a]) for i in range(21)) for a in range(3)])
    size=math.hypot(*derivative)
    return derivative/size if size>0 and math.isfinite(size) else None


def _covered(cells,h,a,b):
    p=a/h;d=(b-a)/h
    lower=np.minimum(p,p+d);upper=np.maximum(p,p+d)
    keep=np.all((cells<=upper)&(cells+1>=lower),axis=1);boxes=cells[keep]
    lo=np.zeros(len(boxes));hi=np.ones(len(boxes));valid=np.ones(len(boxes),bool)
    for axis in range(3):
        if d[axis]==0:valid&=(boxes[:,axis]<=p[axis])&(p[axis]<=boxes[:,axis]+1)
        else:
            one=(boxes[:,axis]-p[axis])/d[axis];two=(boxes[:,axis]+1-p[axis])/d[axis]
            lo=np.maximum(lo,np.minimum(one,two));hi=np.minimum(hi,np.maximum(one,two))
    # Numerical slab union, not an exact rational or physical-connectivity proof.
    tol=256*np.finfo(float).eps*max(1.,float(np.max(np.abs(p))),float(np.max(np.abs(p+d))))
    valid&=hi>=lo-tol;end=0.
    for entry,exit in sorted(zip(lo[valid],hi[valid])):
        if entry>end+tol:return False
        end=max(end,float(exit))
        if end>=1-tol:return True
    return False


def section_centroid_measurement(solid,points,adjacency,*,progress=None):
    """Keep ID graph and endpoints; derive transverse section centroids once."""
    p=np.asarray(points,float);base=graph_measurement_ledger(p,adjacency,coordinate_units='caller_solid_units')
    h=solid.pitch;engine=centroid_section_class()(solid)
    normals={};ambiguous=set()
    for branch in base['branches']:
        ids=branch['source_node_ids'];path=p[ids]
        for at,node in enumerate(ids):
            if len(adjacency[node])>=3:continue
            direction=_normal(path,at,h,branch['closed'])
            if direction is None:continue
            if node in normals and not np.allclose(abs(float(normals[node]@direction)),1.,atol=1e-12):ambiguous.add(node)
            else:normals[node]=direction
    positions=p.copy();radii=[None]*len(p);statuses=[];centroids=[]
    movable={i for i,ns in enumerate(adjacency) if len(ns)==2}
    occupied=set(map(tuple,solid.cells))
    for i,point in enumerate(p):
        direction=normals.get(i)
        if direction is None or i in ambiguous:
            statuses.append('MISSING_UNIQUE_BRANCH_NORMAL');centroids.append(None);continue
        result=engine.section(point,direction);statuses.append(result['status'])
        evidence=result.get('section_centroid_evidence');centroids.append(evidence)
        if evidence is not None:
            radii[i]=result['radius'];candidate=np.array(evidence['centroid'])
            if i in movable and tuple(np.floor(candidate/h).astype(np.int64)) in occupied:positions[i]=candidate
        if progress is not None and (i+1)%100==0:progress(dict(nodes_processed=i+1,total_nodes=len(p)))
    rollback=set();passes=0
    while True:
        bad=[]
        for a,b in base['source_edge_ids']:
            if np.array_equal(positions[a],p[a]) and np.array_equal(positions[b],p[b]):continue
            if not _covered(solid.cells,h,positions[a],positions[b]):bad.extend([a,b])
        changed={i for i in bad if not np.array_equal(positions[i],p[i])}
        if not changed:break
        for i in changed:positions[i]=p[i]
        rollback|=changed;passes+=1
        assert passes<=len(p)
    ledger=graph_measurement_ledger(positions,adjacency,coordinate_units='caller_solid_units',radii=radii)
    assert ledger['source_edge_ids']==base['source_edge_ids']
    pinned=[i for i,ns in enumerate(adjacency) if len(ns)!=2]
    assert np.array_equal(positions[pinned],p[pinned])
    return dict(status='COMPLETE_EXPERIMENTAL_SECTION_CENTROID_GRAPH',source_graph_ledger=base,
        centroid_graph_ledger=ledger,source_points=p.tolist(),points=positions.tolist(),
        adjacency=[sorted(map(int,ns)) for ns in adjacency],area_equivalent_radii=radii,
        section_status_counts=dict(Counter(statuses)),section_statuses=statuses,section_centroids=centroids,
        endpoint_and_junction_coordinates_exact=True,all_source_IDs_and_edges_preserved=True,
        source_solid_modified=False,GT_axis_or_radius_used=False,passes=1,
        numerical_edge_support_rollback_node_ids=sorted(rollback),rollback_passes=passes,
        newly_shifted_edges_with_numerical_closed_cube_support=True,
        moved_nodes=int(np.count_nonzero(np.any(positions!=p,axis=1))),
        radius_definition='area-equivalent radius of selected source cube-union plane component',
        physical_tip_extension=0.,production_metric_adopted=False,biological_measurements_admitted=False)
