"""Replace a failed long terminal route using a supported interior anchor."""
import math
import numpy as np
from root_graph_measurement_ledger import graph_measurement_ledger
from root_section_centroid import _covered


def retry_short_boundary(solid,points,adj,records,ownership):
    """Admit a local first boundary of the fixed graph-end direction.

    Does not assert a pure section, physical root tip or sub-pitch accuracy.
    All source/ownership/direction challenges are retained in the receipt.
    """
    from root_reconstruction_measurement import cubical_ray_exit
    h=solid.pitch;ledger=graph_measurement_ledger(points,adj,coordinate_units='model_units')
    for record in records:
        if record['status']=='geometric_surface_estimate':continue
        node=record['node_id']
        branch=next((b for b in ledger['branches'] if node in [b['source_node_ids'][0],b['source_node_ids'][-1]]),None)
        if branch is None or branch['closed'] or not 0<branch['source_length']<12*h:continue
        own=ownership.node_branches[node] if node<len(ownership.node_branches) else set()
        if len(own)!=1 or len(adj[node])!=1:continue
        neighbor=adj[node][0];q=np.asarray(points[node],float);direction=q-np.asarray(points[neighbor],float)
        norm=math.hypot(*direction)
        if norm<=0:continue
        direction/=norm;axis=np.eye(3)[int(np.argmin(np.abs(direction)))]
        u=np.cross(direction,axis);u/=math.hypot(*u);v=np.cross(direction,u)
        directions=[direction]+[math.cos(math.pi/18)*direction+sign*math.sin(math.pi/18)*a for a in [u,v] for sign in [-1,1]]
        probes=[];refusal=None
        for at,n in enumerate(directions):
            try:ray=cubical_ray_exit(solid,q,n)
            except ValueError as exc:
                reason=str(exc)
                probes.append(dict(direction=n.tolist(),ray=None,refusals=[reason]))
                if refusal is None:refusal=reason
                continue
            end=np.asarray(ray['boundary_point'],float)
            failures=[]
            if ray['boundary_distance']>(2 if at==0 else 4)*h:
                failures.append('SHORT_ENDPOINT_EXIT_NOT_LOCAL')
            supported=bool(_covered(solid.cells,h,q,end))
            if not supported:failures.append('SHORT_ENDPOINT_EXIT_OUTSIDE_SOURCE')
            owners=[sorted(ownership.nearest_owners(q+t*(end-q))) for t in np.linspace(0,.99,11)]
            if any(set(o)!=own for o in owners):
                failures.append('SHORT_ENDPOINT_FOREIGN_BRANCH_CHALLENGE')
            if failures and refusal is None:refusal=failures[0]
            probes.append(dict(direction=n.tolist(),ray=ray,owners=owners,
                entire_segment_source_supported=supported,refusals=failures))
        spread=(max(math.dist(a['ray']['boundary_point'],b['ray']['boundary_point'])
                    for a in probes for b in probes) if all(p['ray'] is not None for p in probes) else None)
        if refusal is None and spread>h:refusal='SHORT_ENDPOINT_DIRECTION_CHALLENGE_UNSTABLE'
        record['short_boundary_evidence']=dict(probes=probes,refusal=refusal,challenge_angle_degrees=10,
            boundary_spread=spread,native_pitch=h,stability_not_accuracy_certificate=True,
            all_five_fixed_directions_evaluated=True,
            spread_pass=spread is not None and spread<=h,
            thresholds=dict(central_exit_pitches=2,challenge_exit_pitches=4,spread_pitches=1,
                            branch_max_pitches=12,owner_witnesses_per_ray=11))
        if refusal:continue
        end=probes[0]['ray']['boundary_point']
        if not _covered(solid.cells,h,np.asarray(points[neighbor]),np.asarray(end)):
            record['short_boundary_evidence']['refusal']='SHORT_ENDPOINT_CHANGED_EDGE_OUTSIDE_SOURCE';continue
        points[node]=end
        record.update(status='geometric_surface_estimate',initial_failure=record.get('reason'),
            position=end,direction=direction.tolist(),extension=probes[0]['ray']['boundary_distance'],
            endpoint_method='source_owned_short_graph_end_boundary_challenge',
            ray_status='FIRST_FINITE_BOUNDARY_OF_FIXED_GRAPH_END_DIRECTION',
            region='finite delivered cube union',physical_tip_qualified=False,
            pure_terminal_section_qualified=False,original_entity_unchanged=True,
            retry_status='RECOVERED_BY_LOCAL_SOURCE_AND_DIRECTION_CHALLENGE')
        record.pop('reason',None)
    return points,adj,records


def retry_from_interior(solid,graph,records,ownership):
    from root_native_entity_measurement import geometric_endpoints
    points=[list(p) for p in graph['points']];adj=[list(a) for a in graph['adjacency']];h=solid.pitch
    before=graph_measurement_ledger(points,adj,coordinate_units='model_units')
    for record in records:
        if record['status']=='geometric_surface_estimate':continue
        node=record['node_id']
        branch=next((b for b in before['branches'] if node in [b['source_node_ids'][0],b['source_node_ids'][-1]]),None)
        if branch is None or branch['closed']:continue
        ids=branch['source_node_ids'];ids=ids if ids[0]==node else ids[::-1]
        path=np.asarray([points[i] for i in ids]);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))]
        if arc[-1]<12*h:
            record['retry_status']='INSUFFICIENT_TERMINAL_PORT_SPAN'
            record['interior_retry_evidence']=dict(branch_length=float(arc[-1]),native_pitch=h,
                minimum_branch_pitches=12,candidate_count=0,attempts=[],
                status='INELIGIBLE_UNDER_FIXED_PORT_SPAN_RULE')
            continue
        candidates=[k for k in range(1,len(ids)-1) if 3*h<=arc[k]<=min(8*h,arc[-1]/4)]
        if len(candidates)>12:candidates=candidates[:12]
        receipt=dict(branch_length=float(arc[-1]),native_pitch=h,minimum_branch_pitches=12,
                     candidate_count=len(candidates),attempts=[],status='NO_QUALIFIED_ANCHOR')
        record['interior_retry_evidence']=receipt
        for k in candidates:
            trimmed=[list(a) for a in adj]
            for a,b in zip(ids[:k],ids[1:k+1]):trimmed[a].remove(b);trimmed[b].remove(a)
            proposal,attempts=geometric_endpoints(solid,dict(points=points,adjacency=trimmed),
                _ownership=ownership,_only_endpoints={ids[k]},_retry=False)
            success=next((r for r in attempts if r['node_id']==ids[k] and r['status']=='geometric_surface_estimate'),None)
            receipt['attempts'].append(dict(anchor_node=ids[k],arc=float(arc[k]),
                source_trace=next((r for r in attempts if r['node_id']==ids[k]),None)))
            if success is None:continue
            target=ids[k+1];chain=[ids[k]];previous=None
            while chain[-1]!=target:
                following=[n for n in proposal['adjacency'][chain[-1]] if n!=previous]
                if len(following)!=1:break
                previous=chain[-1];chain.append(following[0])
                if len(chain)>len(proposal['points']):raise AssertionError('TERMINAL_RETRACE_CYCLE')
            if chain[-1]!=target:continue
            q=np.array([proposal['points'][i] for i in chain]);s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(q,axis=0),axis=1))]
            if s[-1]<=0 or np.any(np.diff(s)<=0):continue
            if not all(_covered(solid.cells,h,a,b) for a,b in zip(q,q[1:])):continue
            boundaries=np.linspace(0,float(s[-1]),k+2)
            for at,i in zip(boundaries,ids[:k+2]):points[i]=[float(np.interp(at,s,q[:,a])) for a in range(3)]
            for t,(a,b) in enumerate(zip(ids[:k+1],ids[1:k+2])):
                adj[a].remove(b);adj[b].remove(a);route=[a]
                for j in np.flatnonzero((s>boundaries[t]+1e-13*h)&(s<boundaries[t+1]-1e-13*h)):
                    new=len(points);points.append(q[j].tolist());adj.append([]);route.append(new)
                route.append(b)
                for u,v in zip(route,route[1:]):adj[u].append(v);adj[v].append(u)
            initial=record.get('reason');original=record['original_position']
            record.update(success,node_id=node,original_position=original,initial_failure=initial,
                          endpoint_method='source_owned_interior_anchor_retrace',retry_anchor_node=ids[k],
                          retry_status='RECOVERED_WITH_SUPPORTED_TERMINAL_ROUTE',
                          extension=float(s[-1]-arc[k+1]),new_terminal_length=float(s[-1]),
                          replaced_terminal_length=float(arc[k+1]))
            record.pop('reason',None)
            record['interior_retry_evidence']=dict(receipt,status='QUALIFIED_ANCHOR_RECOVERED')
            break
        if record['status']!='geometric_surface_estimate':record.setdefault('retry_status','NO_QUALIFIED_INTERIOR_ANCHOR_TRACE')
    points,adj,records=retry_short_boundary(solid,points,adj,records,ownership)
    for record in records:
        if record['status']=='geometric_surface_estimate':continue
        record['endpoint_adjudication']=dict(
            status='UNRESOLVED_UNDER_FIXED_ESTIMATOR_POLICY',
            pure_source_trace=dict(reason=record.get('reason'),
                                   initial_section=record.get('initial_section_evidence')),
            interior_anchor=record.get('interior_retry_evidence'),
            short_boundary=record.get('short_boundary_evidence'),
            exhaustive_scope='the existing source trace, eligible interior anchors and five fixed short-boundary directions only',
            geometric_impossibility_proven=False,numeric_endpoint_adopted=False,
            full_length_acceptance='BLOCKED',radius_qualification_changed=False)
    after=graph_measurement_ledger(points,adj,coordinate_units='model_units')
    for key in ['source_components','source_cycle_rank','endpoint_node_ids']:
        assert after[key]==before[key],('TERMINAL_RETRY_TOPOLOGY_CHANGED',key)
    assert len(after['branches'])==len(before['branches'])
    return dict(points=points,adjacency=adj),records
