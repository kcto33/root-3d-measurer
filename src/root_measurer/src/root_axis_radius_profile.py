"""Branch-associated geometric area-equivalent radii on a supplied ID axis."""
from collections import Counter
import math
import numpy as np
from root_graph_measurement_ledger import graph_measurement_ledger
from root_section_centroid import _normal
from root_branch_section_measurement import BranchOwnership,owned_section_class
from root_cubical_disk_section import interior_boundary_section_class


def axis_radius_profile(solid,points,adjacency,*,compare_previous=False,progress=None):
    p=np.asarray(points,float);ledger=graph_measurement_ledger(p,adjacency,coordinate_units='caller_solid_units')
    owner=BranchOwnership(p,adjacency);old_engine=owned_section_class()(solid) if compare_previous else None
    engine=interior_boundary_section_class()(solid);normals={};ambiguous=set()
    for branch in ledger['branches']:
        ids=branch['source_node_ids'];path=p[ids]
        for at,node in enumerate(ids):
            if len(adjacency[node])>=3:continue
            n=_normal(path,at,solid.pitch,branch['closed'])
            if n is None:continue
            if node in normals and not np.isclose(abs(float(normals[node]@n)),1.,atol=1e-12):ambiguous.add(node)
            else:normals[node]=n
    weights=np.zeros(len(p))
    for a,b in ledger['source_edge_ids']:
        length=math.dist(p[a],p[b]);weights[a]+=length/2;weights[b]+=length/2
    rows=[];before=[];after=[]
    def radius(result,node):
        if result['radius'] is None:return None,result['status']
        for q in result['section_centroid_evidence']['_interior_polygon_witnesses']:
            if owner.nearest_owners(q)!=owner.node_branches[node]:return None,'MISSING_MIXED_BRANCH_SECTION'
        return result['radius'],result['status']
    for node,point in enumerate(p):
        n=normals.get(node)
        if n is None or node in ambiguous:
            row=dict(node_id=node,position=point.tolist(),normal=None,old_radius=None,new_radius=None,
                old_status='MISSING_UNIQUE_BRANCH_NORMAL',new_status='MISSING_UNIQUE_BRANCH_NORMAL')
        else:
            new=engine.section(point,n);old=old_engine.section(point,n) if old_engine is not None else None
            if old is not None and old['radius'] is not None:
                assert new['radius']==old['radius'] and new['area']==old['area'],'PREVIOUSLY_ACCEPTED_RADIUS_CHANGED'
            nr,ns=radius(new,node);orr,os=(radius(old,node) if old is not None else (None,'NOT_REQUESTED'))
            if orr is not None:assert nr==orr
            row=dict(node_id=node,position=point.tolist(),normal=n.tolist(),old_radius=orr,new_radius=nr,
                old_status=os,new_status=ns,old_kernel_status=None if old is None else old['status'],
                new_kernel_status=new['status'],query_boundary_resolution=new.get('query_boundary_resolution'))
        before.append(row['old_radius']);after.append(row['new_radius']);rows.append(row)
        if progress is not None and (node+1)%100==0:progress(dict(processed=node+1,total=len(p)))
    def coverage(values):
        available=np.array([v is not None for v in values]);total=math.fsum(map(float,weights))
        return dict(nodes_available=int(available.sum()),nodes=len(values),
            node_fraction=float(available.mean()) if len(available) else None,
            arc_fraction=math.fsum(map(float,weights[available]))/total if total else None)
    return dict(status='COMPLETE_BRANCH_ASSOCIATED_GEOMETRIC_AXIS_RADII',points=p.tolist(),rows=rows,
        radii=after,previous_radii=before if compare_previous else None,node_halfedge_arc_weights=weights.tolist(),
        new_coverage=coverage(after),old_coverage=coverage(before) if compare_previous else None,
        status_counts=dict(Counter(r['new_status'] for r in rows)),
        shared_previously_accepted_radii_exact=True if compare_previous else None,
        source_graph_ledger=ledger,source_geometry_and_IDs_unchanged=True,GT_parameters_accepted=False,
        radius_definition='area-equivalent radius of one branch-associated closed disk section of the cube union',
        physical_or_biological_radius_qualified=False)
