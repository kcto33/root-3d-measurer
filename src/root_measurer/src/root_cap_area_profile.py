"""Round-cap area profiles independent of closed-disk radius admission."""
import math
import types
import numpy as np
from root_cube_plane_integrals import cube_union_plane_integrals
from root_branch_section_measurement import BranchOwnership
import root_round_cap_measurement as original


def cap_area_profile_measurement(solid,points,adjacency,*,round_caps_assumed,body_points=None):
    p=np.asarray(points,float);ownership=BranchOwnership(p,adjacency);events=[]
    class AreaProfiles:
        def __init__(self,source):self.source=source;self.owner=None;self.leaf=None
        def section(self,point,normal):
            matches=np.flatnonzero(np.all(p==point,axis=1))
            leaves=[int(i) for i in matches if len(adjacency[i])==1]
            if leaves:
                self.leaf=leaves[0] if len(leaves)==1 else None
                self.owner=ownership.node_branches[self.leaf] if self.leaf is not None else None
            integral=cube_union_plane_integrals(self.source,point,normal)
            area=integral['area'];reason=None
            if area is None or area<=0:reason=integral['status']
            elif integral['plane_cell_components']!=1:reason='MISSING_MULTIPLE_CAP_PLANE_CELL_GROUPS'
            elif self.owner is None:reason='MISSING_UNIQUE_CAP_SOURCE_BRANCH'
            else:
                for q in integral['polygon_witnesses']:
                    if ownership.nearest_owners(q)!=self.owner:
                        reason='MISSING_FOREIGN_CAP_SECTION_WITNESS';break
            events.append(dict(node_id=self.leaf,status=reason or integral['status'],area=area,
                disk_radius_admission_used=False,groups=integral.get('plane_cell_components')))
            if reason:return dict(area=area,radius=None,status=reason)
            # Compatibility key is an internal sampling scale, not an admitted
            # cross-section radius. Public results retain this distinction.
            return dict(area=area,radius=math.sqrt(area/math.pi),status='COMPLETE_CAP_MODEL_AREA_SAMPLE',
                section_centroid_evidence=dict(centroid=integral['centroid']))
    namespace=dict(vars(original));namespace['centroid_section_class']=lambda:AreaProfiles
    method=types.FunctionType(original.round_cap_graph_measurement.__code__,namespace,
        name='round_cap_graph_measurement',argdefs=original.round_cap_graph_measurement.__defaults__)
    method.__kwdefaults__=original.round_cap_graph_measurement.__kwdefaults__
    assert method.__code__ is original.round_cap_graph_measurement.__code__
    result=method(solid,p,adjacency,round_caps_assumed=round_caps_assumed,body_points=body_points)
    result['area_profile_evidence']=dict(events=events,original_round_cap_codeobject_exact=True,
        complete_area_integral_not_closed_disk_radius=True,original_radius_admission_gates_unchanged=True,
        source_material_or_unknowns_modified=False)
    result['status']='COMPLETE_CONDITIONAL_ROUND_CAP_FROM_AREA_PROFILES'
    return result
