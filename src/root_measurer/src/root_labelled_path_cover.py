"""Exact complete labelled one-dimensional measurement paths, never graph repair.

Older geometry branches bundle XYZ/radii and lose ID/isolated provenance. This
provider only references original node/edge IDs and proves their full cover.
"""
import hashlib,os


def labelled_path_cover(adjacency,*,expected_nodes,expected_edges,execution_domain='local',max_work=10000000):
    if (type(expected_nodes)is not int or not 0<expected_nodes<=1200000 or type(expected_edges)is not int or
            not 0<=expected_edges<=2000000 or len(adjacency)!=expected_nodes):
        raise ValueError('COMPLETE_ORIGINAL_GRAPH_DENOMINATORS_REQUIRED')
    if execution_domain=='local':
        if expected_nodes>256 or expected_edges>1000:raise ValueError('UNCHANGED_LOCAL_256_NODE_1000_EDGE_GUARDS')
    elif execution_domain!='registered_slurm_complete_graph' or not os.environ.get('SLURM_JOB_ID','').isdigit():
        raise ValueError('REAL_REGISTERED_SLURM_REQUIRED_FOR_COMPLETE_GRAPH')
    if type(max_work)is not int or not 0<=max_work<=10000000:raise ValueError('WORK_GUARD_ONLY_TIGHTENED')
    work=0
    def use(n=1):
        nonlocal work
        work+=n
        if work>max_work:raise ValueError('NO_PARTIAL_COVER_AT_WORK_GUARD')
    adj=[]
    for i,neighbors in enumerate(adjacency):
        use();row=list(neighbors);use(len(row))
        if any(type(j)is not int or not 0<=j<expected_nodes or j==i for j in row) or len(set(row))!=len(row):
            raise ValueError('UNIQUE_INTEGER_ORIGINAL_ADJACENCY_REQUIRED')
        adj.append(set(row))
    for i,row in enumerate(adj):
        for j in row:
            use()
            if i not in adj[j]:raise ValueError('COMPLETE_SYMMETRIC_ADJACENCY_REQUIRED')
    edges=[(i,j) for i,row in enumerate(adj) for j in sorted(row) if i<j]
    if len(edges)!=expected_edges:raise ValueError('COMPLETE_EDGE_DENOMINATOR_REQUIRED')
    ids={edge:i for i,edge in enumerate(edges)};degree=list(map(len,adj));labels=[-1]*expected_nodes;components=0
    for seed in range(expected_nodes):
        if labels[seed]>=0:continue
        labels[seed]=components;todo=[seed]
        while todo:
            node=todo.pop();use()
            for neighbor in adj[node]:
                use()
                if labels[neighbor]<0:labels[neighbor]=components;todo.append(neighbor)
        components+=1
    used=set();paths=[]
    def walk(start,first):
        node_ids=[start,first];edge_ids=[];a,b=start,first
        while True:
            use();key=tuple(sorted((a,b)));eid=ids[key]
            if eid in used:raise ValueError('DUPLICATE_PATH_EDGE_BEFORE_COMPLETE_COVER')
            used.add(eid);edge_ids.append(eid)
            if b==start or degree[b]!=2:break
            nxt=next(j for j in adj[b] if j!=a);a,b=b,nxt;node_ids.append(b)
        closed=node_ids[0]==node_ids[-1]
        kind='pure_cycle' if closed and degree[start]==2 else 'junction_return_loop' if closed else 'open_chain'
        paths.append(dict(path_id='labelled_path%06d'%len(paths),original_node_ids=node_ids,original_edge_ids=edge_ids,
            component_id=labels[start],closed=closed,kind=kind,edge_count=len(edge_ids),
            endpoint_original_node_ids=[node_ids[0],node_ids[-1]],
            centering_admitted=None,section_radius_profile=None,biological_branch_ownership=None))
    for start in range(expected_nodes):
        if degree[start]==2:continue
        for first in sorted(adj[start]):
            use()
            if ids[tuple(sorted((start,first)))] not in used:walk(start,first)
    for start in range(expected_nodes):
        for first in sorted(adj[start]):
            use()
            if ids[tuple(sorted((start,first)))] not in used:walk(start,first)
    isolated=[i for i,d in enumerate(degree) if d==0]
    owned=[i for path in paths for i in path['original_edge_ids']]
    represented=set(isolated)
    for path in paths:
        represented.update(path['original_node_ids'])
        internal=path['original_node_ids'][1:-1]
        if any(degree[i]!=2 for i in internal):raise ValueError('INTERNAL_JUNCTION_MUST_NOT_BE_SKIPPED')
        if path['kind']=='open_chain' and any(degree[i]==2 for i in path['endpoint_original_node_ids']):
            raise ValueError('NONMAXIMAL_OPEN_CHAIN')
    if len(owned)!=len(set(owned)):
        raise ValueError('DUPLICATE_EDGE_IN_PATH_COVER')
    if set(owned)!=set(range(expected_edges)) or represented!=set(range(expected_nodes)):
        raise ValueError('ALL_ORIGINAL_EDGES_AND_NODES_REQUIRED')
    pure=sum(p['kind']=='pure_cycle' for p in paths);critical=sum(d!=2 for d in degree)
    quotient_vertices=critical+pure;rank=expected_edges-expected_nodes+components
    if len(paths)-quotient_vertices+components!=rank:raise ValueError('QUOTIENT_CYCLE_RANK_MUST_BE_PRESERVED')
    digest=hashlib.sha256(''.join('%d:%d\n'%e for e in edges).encode()).hexdigest()
    return dict(available=True,status='COMPLETE_LABELLED_MAXIMAL_PATH_CYCLE_COVER',nodes=expected_nodes,edges=expected_edges,
        paths=paths,path_count=len(paths),isolated_original_node_ids=isolated,components=components,
        critical_original_node_ids=[i for i,d in enumerate(degree) if d!=2],pure_cycle_count=pure,
        graph_cycle_rank=rank,quotient_cycle_rank=rank,quotient_vertex_count=quotient_vertices,
        canonical_edge_order_sha256=digest,all_original_edges_owned_exactly_once=True,all_original_nodes_represented=True,
        original_graph_modified=False,coordinates_or_old_radii_read=False,complete_centerline_admitted=None,
        fixed_logical_work=work,work_is_not_CPU_instruction_count=True,physical_length=None,actual_GT_accuracy=None)
