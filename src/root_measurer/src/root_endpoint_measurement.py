"""ID-preserving geometric endpoint measurement with explicit qualification."""
import math
import numpy as np
from root_terminal_tangent_model import terminal_model_diagnostic


def measure_graph_endpoints(points, adjacency, pitch, *, coordinate_units):
    p=np.asarray(points,dtype=float)
    if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():raise ValueError('FINITE_GRAPH_POINTS_REQUIRED')
    if not isinstance(coordinate_units,str) or not coordinate_units.strip():raise ValueError('EXPLICIT_SOURCE_UNIT_LABEL_REQUIRED')
    if isinstance(pitch,bool) or not math.isfinite(pitch) or pitch<=0:raise ValueError('POSITIVE_FINITE_PITCH_REQUIRED')
    if len(adjacency)!=len(p):raise ValueError('GRAPH_NODE_COUNT')
    graph=[]
    for i,neighbors in enumerate(adjacency):
        if any(isinstance(j,(bool,np.bool_)) or not isinstance(j,(int,np.integer)) for j in neighbors):raise ValueError('INTEGER_NODE_IDS_REQUIRED')
        row=set(map(int,neighbors))
        if i in row or any(j<0 or j>=len(p) for j in row):raise ValueError('INVALID_GRAPH_EDGE_ID')
        graph.append(row)
    if any(i not in graph[j] for i,ns in enumerate(graph) for j in ns):raise ValueError('SYMMETRIC_SOURCE_GRAPH_REQUIRED')
    edges=[(i,j) for i,ns in enumerate(graph) for j in ns if i<j]
    lengths=[math.hypot(*(p[j]-p[i])) for i,j in edges]
    if not all(math.isfinite(v) for v in lengths):raise ValueError('GRAPH_EDGE_LENGTH_FLOAT_RANGE')
    remaining=set(range(len(p)));components=[]
    while remaining:
        stack=[min(remaining)];component=set()
        while stack:
            i=stack.pop()
            if i in component:continue
            component.add(i);stack.extend(graph[i]-component)
        remaining-=component;components.append(sorted(component))
    endpoint_ids=[i for i,ns in enumerate(graph) if len(ns)==1];rows=[]
    for leaf in endpoint_ids:
        path=[leaf];previous=None;current=leaf
        while True:
            choices=graph[current]-({previous} if previous is not None else set())
            if not choices:break
            if len(choices)!=1:raise AssertionError('TERMINAL_CHAIN_BRANCH_SELECTION')
            following=next(iter(choices))
            if following in path:raise AssertionError('TERMINAL_CHAIN_CYCLE_REVISIT')
            path.append(following);previous,current=current,following
            if len(graph[current])!=2:break
        try:
            diagnostic=terminal_model_diagnostic(p[path],pitch);failure=None
        except ValueError as error:
            diagnostic=None;failure=str(error)
        rows.append(dict(node_id=leaf,source_chain_node_ids=path,position=p[leaf].tolist(),
            position_definition='untyped degree1 source graph vertex',
            inward_geometric_direction=None if diagnostic is None else diagnostic['linear']['direction'],
            model_sensitivity=diagnostic,numerical_missing_reason=failure,
            physical_root_tip=None,biological_tip_identity_qualified=False,
            statistical_position_confidence_interval=None,physical_tangent_error_bound=None,
            endpoint_moved=False,extension=0.))
    return dict(status='COMPLETE_ID_PRESERVING_GEOMETRIC_ENDPOINT_BUNDLE',coordinate_units=coordinate_units,
        pitch=float(pitch),source_nodes=len(p),source_edges=len(edges),source_components=len(components),
        source_cycle_rank=len(edges)-len(p)+len(components),source_graph_length=math.fsum(lengths),
        source_graph_length_definition='sum of original unique undirected geometric edges',
        physical_total_root_length=None,display_smoothing_applied=False,source_graph_rewritten=False,
        graph_endpoint_count=len(rows),endpoints=rows,all_degree1_node_ids_once=True,
        isolated_node_ids=[i for i,ns in enumerate(graph) if not ns],
        isolated_points=[p[i].tolist() for i,ns in enumerate(graph) if not ns],
        geometrically_coincident_ids_not_merged=True,
        numerical_direction_available_count=sum(r['inward_geometric_direction'] is not None for r in rows),
        model_selected_by_fit_residual=False,physical_measurements_admitted=False)
