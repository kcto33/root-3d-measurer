"""One bounded subvoxel fairing step, with native geometry retained explicitly.

Connectivity, branch endpoints and native radius queries are retained. Fitted
coordinates are geometric estimates with half-voxel displacement uncertainty.
"""
from collections import defaultdict
import numpy as np

def fit_polylines(decoded):
    if decoded['status']!='OK':return decoded
    pitch=float(decoded['audit']['pitch'])
    if not np.isfinite(pitch) or pitch<=0:raise ValueError('POLYLINE_PITCH')
    key=lambda p:tuple(np.round(p,12));adj=defaultdict(set)
    for b in decoded['branches']:
        for a,z in zip(b['nodes'],b['nodes'][1:]):adj[key(a)].add(key(z));adj[key(z)].add(key(a))
    branches=[];mapping={};maximum=0.;source_length=0.;fitted_length=0.
    for source in decoded['branches']:
        p=np.asarray(source['nodes'],float);q=p.copy()
        if len(p)>2:
            closed=key(p[0])==key(p[-1]) and len(adj[key(p[0])])==2
            if closed:
                unique=p[:-1];v=.25*np.roll(unique,1,axis=0)+.5*unique+.25*np.roll(unique,-1,axis=0)
                q=np.r_[v,v[:1]]
            else:q[1:-1]=.25*p[:-2]+.5*p[1:-1]+.25*p[2:]
            shift=q-p;norm=np.linalg.norm(shift,axis=1);factor=np.minimum(1.,.5*pitch/np.maximum(norm,1e-30));q=p+shift*factor[:,None]
            maximum=max(maximum,float(np.linalg.norm(q-p,axis=1).max()))
        for old,new in zip(p,q):
            k=key(old);v=key(new)
            if k in mapping and mapping[k]!=v:raise ValueError('INCONSISTENT_SHARED_POLYLINE_NODE')
            mapping[k]=v
        b=dict(source);b['native_nodes']=source['nodes'];b['native_trace_length']=source['arc_length'];b['nodes']=q.tolist()
        b['arc_length']=float(np.linalg.norm(np.diff(q,axis=0),axis=1).sum());source_length+=source['arc_length'];fitted_length+=b['arc_length'];branches.append(b)
    if len(set(mapping.values()))!=len(mapping):raise ValueError('FITTED_NODE_COLLISION_NO_JOIN_ALLOWED')
    # Source connectivity is identified by the one-to-one native -> fitted map.
    original_edges={tuple(sorted((a,z))) for a,ns in adj.items() for z in ns}
    fitted_edges={tuple(sorted((mapping[a],mapping[z]))) for a,z in original_edges}
    if len(original_edges)!=len(fitted_edges):raise ValueError('FITTED_CONNECTIVITY_LOSS')
    if maximum>.5*pitch+1e-10:raise AssertionError('FITTED_DISPLACEMENT_BOUND')
    audit=dict(decoded['audit']);audit['polyline_geometry']={'method':'one symmetric three-point quarter-half-quarter fairing with halfvoxel cap',
        'native_trace_length':source_length,'fitted_curve_length':fitted_length,'maximum_displacement':maximum,'displacement_cap':.5*pitch,
        'native_nodes_retained':True,'native_radius_values_unchanged':True,'radius_query_locations':'native samples associated by index; not exact distance at shifted fitted points',
        'radius_location_uncertainty_bound':.5*pitch,'physical_edges':len(original_edges),'new_edges':0,'components_joined':0,
        'junction_and_openend_positions_unchanged':True,'iterations':1,'biological_accuracy_qualified':False}
    result=dict(decoded);result.update(branches=branches,total_length=fitted_length,audit=audit,native_adjacency_audit=audit)
    return result
