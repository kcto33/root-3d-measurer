"""Graph representative of cube-nerve H1; retain tunnels and all unfilled edges.

A deterministic spanning forest supplies cycle coordinates. Filled cube
intersection triangle boundaries are reduced over GF(2); only their pivot
chords are removed. The remaining graph cycle space maps isomorphically onto
H1 of the cube nerve. H2/cavities are reported separately, not filled away.
This is a graph interpretation of an existing skeleton, not a thinning method.
"""
import numpy as np
from root_reconstruction_measurement import neighbor_graph
from root_measurement_graph_topology import skeleton_voxel_h1
from root_volume_topology import cubical_euler


def quotient_graph(cells,pitch):
    cells=np.asarray(cells,np.int64)
    if cells.ndim!=2 or cells.shape[1]!=3 or pitch<=0:raise ValueError('CUBICAL_GRAPH_INPUT')
    pairs,labels=neighbor_graph(cells)
    edges=sorted(tuple(map(int,e)) for e in pairs)
    adj=[set() for _ in cells]
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    parent=list(range(len(cells)))
    def find(a):
        while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
        return a
    forest=set()
    for a,b in sorted(edges,key=lambda e:(int(np.sum((cells[e[0]]-cells[e[1]])**2)),e)):
        u,v=find(a),find(b)
        if u!=v:parent[u]=v;forest.add((a,b))
    chords=[e for e in edges if e not in forest];chord_id={e:i for i,e in enumerate(chords)}
    pivots={};filled_edges=set();triangles=0
    for a,ns in enumerate(adj):
        for b in sorted(n for n in ns if n>a):
            for c in sorted(n for n in ns&adj[b] if n>b):
                triangle=[(a,b),(a,c),(b,c)];filled_edges.update(triangle);triangles+=1
                column=0
                for edge in triangle:
                    if edge in chord_id:column^=1<<chord_id[edge]
                while column:
                    pivot=column.bit_length()-1
                    if pivot in pivots:column^=pivots[pivot]
                    else:pivots[pivot]=column;break
    removed={chords[i] for i in pivots};out=[set() for _ in cells]
    for a,b in edges:
        if (a,b) not in removed:out[a].add(b);out[b].add(a)
    # A newly exposed leaf internal to filled 2D junction support is not a
    # newly discovered arm. Every original end and every unfilled edge stays.
    ports={i for i,ns in enumerate(adj) if len(ns)<=1}
    unfilled=set(edges)-filled_edges
    ports.update(i for edge in unfilled for i in edge)
    pruned=[];queue=[i for i,ns in enumerate(out) if len(ns)==1 and i not in ports]
    while queue:
        a=queue.pop()
        if len(out[a])!=1 or a in ports:continue
        b=next(iter(out[a]));edge=tuple(sorted((a,b)))
        assert edge in filled_edges
        out[a].clear();out[b].remove(a);pruned.append(a)
        if len(out[b])==1 and b not in ports:queue.append(b)
    keep={i for i,ns in enumerate(out) if ns or not adj[i]}
    for label in np.unique(labels):
        members=np.flatnonzero(labels==label)
        if not any(int(i) in keep for i in members):keep.add(int(members[0]))
    ids=sorted(keep);mapping={old:new for new,old in enumerate(ids)}
    graph=[sorted(mapping[j] for j in out[i]) for i in ids]
    homology=skeleton_voxel_h1(cells)
    components=len(np.unique(labels));rank=sum(map(len,graph))//2-len(ids)+components
    assert rank==homology['skeleton_voxel_H1']
    assert all(b in out[a] for a,b in unfilled)
    h2=cubical_euler(cells)['euler']-components+rank
    audit=dict(status='PASS_H1_GRAPH_REPRESENTATIVE',source_nodes=len(cells),graph_nodes=len(ids),
        graph_component_count=components,graph_cycle_rank=rank,skeleton_voxel_H1=rank,skeleton_voxel_H2=h2,
        filled_triangle_count=triangles,boundary_rank=len(pivots),pivot_chords_removed=len(removed),
        passive_filled_region_leaves_removed=len(pruned),unfilled_edges_preserved=True,
        original_terminal_nodes_preserved=True,source_material_modified=False,
        source_voxel_ids=ids,all_material_topology_represented_by_graph=h2==0,
        cavities_filled=0,genuine_tunnels_forced_to_tree=False,biological_axis_qualified=False)
    return (cells[ids]+.5)*pitch,graph,audit
