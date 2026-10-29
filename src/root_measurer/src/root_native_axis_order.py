"""Fixed compiled-Lee axis-order adaptation, without altering material.

The full fixed endpoint rules remain authoritative. Old/new port matching is
a diagnostic of changed estimators, not a reference truth for root identity.
"""
import math
import numpy as np
from scipy.optimize import linear_sum_assignment
from root_graph_measurement_ledger import graph_measurement_ledger

def original_frame_kernel(solid,order,max_dense_cells):
    from root_reconstruction_measurement import ReconstructedSolid
    from root_tube_measurement import measure_tube_solid
    if sorted(order)!=[0,1,2]:raise ValueError('AXIS_PERMUTATION_REQUIRED')
    inverse=np.argsort(order)
    source=ReconstructedSolid(solid.cells[:,order].copy(),solid.pitch,dict(algorithm_axis_order=list(order)))
    kernel=measure_tube_solid(source,max_dense_cells=max_dense_cells)
    kernel['points']=np.asarray(kernel['points'])[:,inverse].tolist()
    kernel['source_skeleton_cells']=np.asarray(kernel['source_skeleton_cells'])[:,inverse].tolist()
    from root_sparse_solid_topology import sparse_solid_topology
    material=sparse_solid_topology(solid.cells,max_domain_cells=max_dense_cells)
    derived=sparse_solid_topology(np.asarray(kernel['source_skeleton_cells']),max_domain_cells=max_dense_cells)
    keys=['foreground_components','tunnel_dimension','enclosed_empty_components','euler']
    if [material[k] for k in keys]!=[derived[k] for k in keys]:
        raise ValueError('AXIS_ORDER_SOURCE_TOPOLOGY_NOT_PRESERVED')
    kernel['thinning']['full_source_topology_check']=dict(source=material,derived=derived,
        B0_B1_B2_Euler_equal=True,source_material_changed=False)
    for item in kernel['representation'].get('corrected_faces',[]):
        item['coordinates']=np.asarray(item['coordinates'])[inverse].tolist()
    kernel['algorithm_axis_order']=dict(order=list(order),inverse=inverse.tolist(),source_geometry_changed=False)
    return kernel


def port_connectivity_gate(before,after,pitch):
    def describe(graph):
        p=np.asarray(graph['points']);adj=graph['adjacency']
        ledger=graph_measurement_ledger(p,adj,coordinate_units='model_units')
        ports=[i for i,a in enumerate(adj) if len(a)!=2]
        return p,adj,ledger,ports
    p,a,l,ports=describe(before);q,b,m,other=describe(after)
    checks=dict(components=l['source_components']==m['source_components'],
        cycles=l['source_cycle_rank']==m['source_cycle_rank'],
        paths=len(l['branches'])==len(m['branches']),
        port_degrees=sorted(len(a[i]) for i in ports)==sorted(len(b[j]) for j in other))
    receipt=dict(checks=checks,maximum_port_displacement_pitches=None,
        port_distance_limit_pitches=4,source_axis_nodes_removed=False)
    if not all(checks.values()):return dict(receipt,qualified=False)
    costs=np.array([[math.dist(p[i],q[j]) if len(a[i])==len(b[j]) else 1e100 for j in other] for i in ports])
    mapping={ports[int(i)]:other[int(j)] for i,j in zip(*linear_sum_assignment(costs))} if ports else {}
    maximum=max((math.dist(p[i],q[j]) for i,j in mapping.items()),default=0.)
    checks['local_port_bijection']=maximum<=4*pitch
    def links(ledger,convert):
        result=[]
        for branch in ledger['branches']:
            if branch['closed']:result.append(('ring',));continue
            ids=branch['source_node_ids'];ends=sorted([convert(ids[0]),convert(ids[-1])])
            result.append(('path',ends[0],ends[1]))
        return sorted(result)
    checks['same_port_connections']=links(l,lambda i:mapping[i])==links(m,lambda i:i)
    receipt.update(maximum_port_displacement_pitches=maximum/pitch,port_mapping=mapping)
    return dict(receipt,qualified=all(checks.values()))


def endpoint_candidate_admission(primary,candidate):
    """Never trade fewer failures/shorter length for a complete endpoint set."""
    from root_model_delivery_validation import complete_length_state
    def signature(result):
        graph=result.get('candidate_centerline')
        if not graph or not graph.get('points'):return None
        ledger=graph_measurement_ledger(graph['points'],graph['adjacency'],coordinate_units='model_units')
        return dict(components=ledger['source_components'],cycles=ledger['source_cycle_rank'],
            paths=len(ledger['branches']),ports=sorted(len(a) for a in graph['adjacency'] if len(a)!=2))
    a=signature(primary);b=signature(candidate)
    original=primary['entity_binding'];new=candidate['entity_binding']
    state=complete_length_state(candidate)
    proof=candidate.get('stages',{}).get('skeleton',{}).get('thinning',{}).get('full_source_topology_check',{})
    checks=dict(same_input=all(original.get(k)==new.get(k) for k in ['cells_sha256','pitch','stage']),
        full_source_topology=proof.get('B0_B1_B2_Euler_equal') is True,
        complete_endpoint_set=b is not None and not state['unresolved_geometric_endpoint_ids'],
        graph_signature_preserved=a is not None and a==b,
        pipeline_completed=not any(f.get('code')=='NATIVE_PIPELINE_BLOCKED' for f in candidate.get('failures',[])))
    return dict(qualified=all(checks.values()),checks=checks,original_signature=a,candidate_signature=b,
        candidate_length_status=state,old_port_identity_not_assumed=True,GT_used=False)
