"""Fixed mesh target and adversarial resolution checks, separate from accuracy."""
import math
import numpy as np
from scipy.spatial import cKDTree
from root_graph_measurement_ledger import graph_measurement_ledger


def resolution_schedule(target_pitch,reason):
    if target_pitch is None or not isinstance(reason,str) or not reason.strip():
        raise ValueError('RESOLUTION_UNRESOLVED_INPUT_FEATURE_SCALE_REQUIRED')
    if isinstance(target_pitch,bool) or not math.isfinite(target_pitch) or target_pitch<=0:
        raise ValueError('POSITIVE_INPUT_DERIVED_TARGET_PITCH_REQUIRED')
    return [dict(name='coarse',pitch=2*target_pitch,origin=[0.,0.,0.]),
            dict(name='target',pitch=target_pitch,origin=[0.,0.,0.]),
            dict(name='finer_challenge',pitch=target_pitch/2,origin=[0.,0.,0.]),
            dict(name='half_grid_phase',pitch=target_pitch,origin=[target_pitch/2]*3)]


def graph_evidence(graph,pitch):
    if len(graph.get('points',[]))==0:
        return dict(components=0,cycles=0,endpoints=0,paths=0,isolated=0,length=0.,witnesses=np.empty((0,3)))
    p=np.asarray(graph['points'],float);adj=graph['adjacency']
    ledger=graph_measurement_ledger(p,adj,coordinate_units='model_units')
    witnesses=[]
    for a,b in ledger['source_edge_ids']:
        count=max(1,math.ceil(math.dist(p[a],p[b])/(pitch/2)))
        witnesses.extend(p[a]+t*(p[b]-p[a]) for t in np.linspace(0,1,count+1))
    return dict(components=ledger['source_components'],cycles=ledger['source_cycle_rank'],
                endpoints=len(ledger['endpoint_node_ids']),paths=len(ledger['branches']),
                isolated=len(ledger['isolated_node_ids']),length=ledger['measurement_graph_length'],
                witnesses=np.asarray(witnesses,float).reshape(-1,3))


def compare_graphs(a,b,pitch):
    first=graph_evidence(a,pitch);second=graph_evidence(b,pitch)
    topology=all(first[k]==second[k] for k in ['components','cycles','endpoints','paths','isolated'])
    length=abs(first['length']-second['length'])/max(first['length'],second['length'],1e-300)
    distances=[]
    for x,y in [(first['witnesses'],second['witnesses']),(second['witnesses'],first['witnesses'])]:
        distances.append(float(cKDTree(y).query(x)[0].max()) if len(x) and len(y) else None)
    coverage=all(d is not None and d<=math.sqrt(3)*pitch for d in distances)
    return dict(topology_and_path_counts_equal=topology,relative_length_change=length,
        numerical_stability_signal_only=length<=.02,bidirectional_graph_sample_max_distance=distances,
        sample_spacing=pitch/2,graph_sample_coverage_pass=coverage,
        passed=topology and length<=.02 and coverage,
        continuous_branch_completeness_or_precision_proved=False)


def qualify_resolution(graphs,target_pitch):
    required=['coarse','target','finer_challenge','half_grid_phase']
    missing=[k for k in required if k not in graphs]
    if missing:return dict(status='RESOLUTION_CHALLENGES_INCOMPLETE',passed=False,missing=missing)
    comparisons={name:compare_graphs(graphs['target'],graphs[name],target_pitch)
                 for name in ['coarse','finer_challenge','half_grid_phase']}
    # Coarse instability is allowed only because the declared target was reached.
    # It cannot be substituted for the mandatory finer and phase challenges.
    passed=all(comparisons[k]['passed'] for k in ['finer_challenge','half_grid_phase'])
    return dict(status='PASSED_FIXED_CHALLENGES_PENDING_OVERLAY' if passed else 'RESOLUTION_CHALLENGE_FAILED',
        passed=passed,comparisons=comparisons,target_pitch=target_pitch,
        automatic_precision_acceptance=False,original_surface_overlay_required=True)
