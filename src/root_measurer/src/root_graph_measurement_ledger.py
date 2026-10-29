"""Source-ID branch and coverage ledger; no smoothing or biological inference."""
import math
import numpy as np


def graph_measurement_ledger(points, adjacency, *, coordinate_units, radii=None, display_points=None):
    p=np.asarray(points,dtype=float)
    if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():raise ValueError('FINITE_N_BY_3_POINTS_REQUIRED')
    if not isinstance(coordinate_units,str) or not coordinate_units.strip():raise ValueError('EXPLICIT_UNITS_REQUIRED')
    if len(adjacency)!=len(p):raise ValueError('ADJACENCY_SIZE')
    graph=[]
    for i,neighbors in enumerate(adjacency):
        row=list(neighbors)
        if any(isinstance(j,(bool,np.bool_)) or not isinstance(j,(int,np.integer)) for j in row):raise ValueError('INTEGER_IDS_REQUIRED')
        row=list(map(int,row))
        if len(set(row))!=len(row):raise ValueError('DUPLICATE_ADJACENCY_ID')
        if i in row or any(j<0 or j>=len(p) for j in row):raise ValueError('INVALID_EDGE_ID')
        graph.append(set(row))
    if any(i not in graph[j] for i,ns in enumerate(graph) for j in ns):raise ValueError('ASYMMETRIC_GRAPH')
    edges=sorted((i,j) for i,ns in enumerate(graph) for j in ns if i<j)
    def distance(a,b):
        with np.errstate(over='ignore',invalid='ignore'):delta=b-a
        result=math.hypot(*delta)
        if not math.isfinite(result):raise ValueError('GEOMETRIC_LENGTH_FLOAT_RANGE')
        return result
    lengths={edge:distance(p[edge[0]],p[edge[1]]) for edge in edges}
    def total(values):
        try:result=math.fsum(values)
        except OverflowError as error:raise ValueError('TOTAL_LENGTH_FLOAT_RANGE') from error
        if not math.isfinite(result):raise ValueError('TOTAL_LENGTH_FLOAT_RANGE')
        return result
    source_length=total(lengths.values())
    labels={};components=[]
    for start in range(len(p)):
        if start in labels:continue
        stack=[start];nodes=[];label=len(components);labels[start]=label
        while stack:
            i=stack.pop();nodes.append(i)
            for j in sorted(graph[i],reverse=True):
                if j not in labels:labels[j]=label;stack.append(j)
        components.append(sorted(nodes))
    used=set();paths=[]
    def walk(start,following):
        path=[start];previous=start;current=following
        while True:
            edge=tuple(sorted((previous,current)))
            if edge in used:raise AssertionError('SOURCE_EDGE_REUSED')
            used.add(edge);path.append(current)
            if current==start or len(graph[current])!=2:break
            following=next(iter(graph[current]-{previous}))
            previous,current=current,following
        return path
    for start,ns in enumerate(graph):
        if len(ns)==2:continue
        for following in sorted(ns):
            if tuple(sorted((start,following))) not in used:paths.append(walk(start,following))
    for start,ns in enumerate(graph):
        for following in sorted(ns):
            if tuple(sorted((start,following))) not in used:paths.append(walk(start,following))
    assert used==set(edges) and sum(len(path)-1 for path in paths)==len(edges)
    rs=[None]*len(p) if radii is None else list(radii)
    if len(rs)!=len(p):raise ValueError('RADIUS_SIZE')
    for i,r in enumerate(rs):
        if r is None:continue
        if isinstance(r,(bool,np.bool_)) or not math.isfinite(float(r)) or float(r)<0:raise ValueError('FINITE_NONNEGATIVE_RADIUS_OR_NONE')
        rs[i]=float(r)
    half_supported=total(lengths[i,j]*(int(rs[i] is not None)+int(rs[j] is not None))/2 for i,j in edges)
    complete_supported=total(lengths[i,j] for i,j in edges if rs[i] is not None and rs[j] is not None)
    branches=[]
    for i,path in enumerate(paths):
        branch_edges=[tuple(sorted(edge)) for edge in zip(path,path[1:])]
        branches.append(dict(branch_id=i,source_node_ids=path,component_id=labels[path[0]],
            closed=path[0]==path[-1],source_edge_ids=[list(edge) for edge in branch_edges],
            source_length=total(lengths[edge] for edge in branch_edges),radii=[rs[j] for j in path]))
    display_length=None;bound=None
    if display_points is not None:
        q=np.asarray(display_points,dtype=float)
        if q.shape!=p.shape or not np.isfinite(q).all():raise ValueError('DISPLAY_POINTS_EXACT_ID_SHAPE')
        display_length=total(distance(q[i],q[j]) for i,j in edges)
        shifts=[distance(a,b) for a,b in zip(p,q)]
        bound=total(shifts[i]+shifts[j] for i,j in edges)
        slack=32*np.finfo(float).eps*max(source_length,display_length,bound)
        if abs(display_length-source_length)>bound+slack:raise AssertionError('DISPLAY_TRIANGLE_BOUND')
    return dict(status='COMPLETE_SOURCE_ID_GRAPH_MEASUREMENT_LEDGER',coordinate_units=coordinate_units,
        source_nodes=len(p),source_edges=len(edges),source_components=len(components),
        source_cycle_rank=len(edges)-len(p)+len(components),component_node_ids=components,
        source_edge_ids=[list(e) for e in edges],branches=branches,
        measurement_graph_length=source_length,serialized_branch_length=total(b['source_length'] for b in branches),
        all_source_edges_owned_exactly_once=True,display_curve_length=display_length,
        display_total_length_displacement_bound=bound,display_never_overwrites_measurement=True,
        endpoint_node_ids=[i for i,ns in enumerate(graph) if len(ns)==1],
        junction_node_ids=[i for i,ns in enumerate(graph) if len(ns)>=3],
        isolated_node_ids=[i for i,ns in enumerate(graph) if not ns],
        zero_length_source_edge_ids=[list(e) for e,v in lengths.items() if v==0],
        radius_definition='caller-supplied geometric radius at original ID; no biological interpretation',
        missing_radius_node_ids=[i for i,r in enumerate(rs) if r is None],
        radius_node_halfedge_supported_length=half_supported,radius_complete_segment_supported_length=complete_supported,
        radius_node_arc_coverage=None if not source_length else half_supported/source_length,
        radius_complete_segment_arc_coverage=None if not source_length else complete_supported/source_length,
        physical_total_root_length=None,physical_radius_qualified=False,
        coordinates_or_graph_modified=False,coincident_IDs_merged=False,unknown_filled=False)
