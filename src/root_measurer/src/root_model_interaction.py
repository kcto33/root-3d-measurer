"""User-specified continuous branch intervals and fresh original-solid sections.

No graph route search, stored-radius interpolation or reference data access.
"""
import math
import numpy as np
from root_model_delivery_exports import file_reference


def number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError('FINITE_NUMBER_REQUIRED')
    return float(value)


class ModelInteraction:
    def __init__(self, measurement):
        self.report = measurement
        graph = measurement.get('candidate_centerline')
        if not graph or not graph.get('points'):
            raise ValueError('MEASURED_AXIS_REQUIRED')
        from root_graph_measurement_ledger import graph_measurement_ledger
        self.points = np.asarray(graph['points'], float)
        self.adjacency = graph['adjacency']
        self.branches = graph_measurement_ledger(self.points, self.adjacency,
                                                 coordinate_units='model_units')['branches']
        self.paths = {}
        for i, branch in enumerate(self.branches):
            points = self.points[branch['source_node_ids']]
            arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
            self.paths[f'path{i:05d}'] = (branch, points, arc)
        self.solid = self.engine = self.owner = self.faces = None
        self.slabs = {}

    def path(self, name):
        if name not in self.paths:
            raise ValueError('EXPLICIT_EXISTING_PATH_REQUIRED')
        return self.paths[name]

    def at(self, name, arc):
        branch, points, cumulative = self.path(name)
        arc = number(arc)
        if not 0 <= arc <= cumulative[-1]:
            raise ValueError('ARC_OUTSIDE_EXPLICIT_PATH')
        p = np.array([np.interp(arc, cumulative, points[:, axis]) for axis in range(3)])
        return p, branch, points, cumulative

    def length(self, intervals):
        if not isinstance(intervals, list) or not 1 <= len(intervals) <= 256:
            raise ValueError('EXPLICIT_ORDERED_INTERVALS_REQUIRED')
        parts = []; result = []; used = {}; previous = None
        for item in intervals:
            name = item['path_id']; start = number(item['start']); end = number(item['end'])
            a, branch, points, arc = self.at(name, start)
            b, _, _, _ = self.at(name, end)
            if start == end:
                raise ValueError('POSITIVE_INTERVAL_REQUIRED')
            tolerance = 128*np.finfo(float).eps*max(1., float(np.abs(self.points).max()))
            if previous is not None and math.dist(previous, a) > tolerance:
                raise ValueError('USER_INTERVALS_ARE_DISCONNECTED_NO_AUTO_ROUTING')
            lo, hi = sorted([start, end])
            if any(min(hi, v) > max(lo, u)+tolerance for u, v in used.get(name, [])):
                raise ValueError('REPEATED_ROOT_INTERVAL_NOT_COUNTED_TWICE')
            used.setdefault(name, []).append((lo, hi))
            interior = points[(arc > lo) & (arc < hi)]
            if end < start: interior = interior[::-1]
            chain = [a.tolist(), *interior.tolist(), b.tolist()]
            result.extend(chain if not result else chain[1:]); previous = b
            parts.append(dict(path_id=name, start=start, end=end, length=abs(end-start)))
        scale = self.report.get('scale_mm_per_model_unit') or 1.
        return dict(kind='USER_SPECIFIED_ROOT_LENGTH', length=math.fsum(p['length'] for p in parts)*scale,
                    unit='mm' if self.report.get('scale_mm_per_model_unit') else 'model_units',
                    intervals=parts, polyline=result, graph_route_searched=False,
                    clicked_screen_trace_used_as_length=False)

    def initialize_solid(self):
        if self.solid is not None:return
        from root_native_entity_measurement import load_entity, NativeSectionIndex
        from root_cubical_disk_section import interior_boundary_section_class
        from root_branch_section_measurement import BranchOwnership
        binding = self.report.get('entity_binding')
        if not binding:raise ValueError('INTERACTIVE_NATIVE_ENTITY_REQUIRED')
        solid, actual, _ = load_entity(binding['metadata']['path'], stage=binding['stage'])
        if actual['cells_sha256'] != binding['cells_sha256']:
            raise ValueError('ORIGINAL_ENTITY_CHANGED')
        self.solid = solid
        self.engine = NativeSectionIndex(solid, interior_boundary_section_class())
        self.owner = BranchOwnership(self.points, self.adjacency)

    def section_outline(self, point, normal):
        # Display-only intersection with original exposed faces; radius comes
        # exclusively from the qualified measurement kernel above.
        from root_reconstruction_measurement import boundary_faces
        if self.faces is None:
            centers, axes = boundary_faces(self.solid)
            corners = []
            for axis in range(3):
                selected = centers[axes == axis].astype(float)
                others = [a for a in range(3) if a != axis]
                offsets = np.zeros((4, 3))
                offsets[:, others] = [[-.5,-.5],[.5,-.5],[.5,.5],[-.5,.5]]
                corners.append((selected[:, None, :]+offsets)*self.solid.pitch)
            self.faces = np.concatenate(corners)
        signed = (self.faces-point)@normal
        tolerance = self.solid.pitch*1e-8
        candidates = self.faces[(signed.min(axis=1)<=tolerance)&(signed.max(axis=1)>=-tolerance)]
        basis0 = np.cross(normal, np.eye(3)[np.argmin(np.abs(normal))]);basis0 /= np.linalg.norm(basis0)
        basis = np.array([basis0, np.cross(normal, basis0)])
        segments = set()
        for face in candidates:
            hits = []
            for a, b in zip(face, np.roll(face, -1, axis=0)):
                da, db = float((a-point)@normal), float((b-point)@normal)
                if abs(da)<=tolerance:hits.append(a)
                if da*db<0:hits.append(a+da/(da-db)*(b-a))
            unique = {tuple(np.round((v-point)@basis.T/tolerance).astype(np.int64)):v for v in hits}
            if len(unique)==2:
                keys=sorted(unique)
                segments.add(tuple(keys))
        adjacency={}
        for a,b in segments:
            adjacency.setdefault(a,set()).add(b);adjacency.setdefault(b,set()).add(a)
        visited=set();containing=[];all_polygons=[]
        for start in sorted(adjacency):
            if start in visited or len(adjacency[start])!=2:continue
            chain=[start];previous=None;current=start
            for _ in range(len(adjacency)+1):
                visited.add(current)
                nxt=next((v for v in sorted(adjacency[current]) if v!=previous),None)
                if nxt is None:break
                chain.append(nxt)
                if nxt==start:break
                previous,current=current,nxt
                if len(adjacency[current])!=2:break
            if chain[-1]!=start:continue
            polygon=np.asarray(chain,float)*tolerance
            all_polygons.append(polygon)
            inside=False
            for a,b in zip(polygon,polygon[1:]):
                if (a[1]>0)!=(b[1]>0) and a[0]+(b[0]-a[0])*(-a[1])/(b[1]-a[1])>0:
                    inside=not inside
            if inside:containing.append(polygon)
        self.last_section_outlines=[]
        if len(containing)!=1:return []
        def contains(poly,q):
            value=False
            for a,b in zip(poly,poly[1:]):
                if (a[1]>q[1])!=(b[1]>q[1]) and a[0]+(b[0]-a[0])*(q[1]-a[1])/(b[1]-a[1])>q[0]:value=not value
            return value
        outer=containing[0]
        inside=[p for p in all_polygons if p is not outer and contains(outer,p[0])]
        holes=[p for p in inside if not any(q is not p and contains(q,p[0]) for q in inside)]
        self.last_section_outlines=[(point+p@basis).tolist() for p in [outer,*holes]]
        return self.last_section_outlines[0]

    def radius(self, request):
        name = request['path_id']; branch, points, arc = self.path(name)
        picked = request.get('position')
        if 'arc' in request:
            location = number(request['arc'])
        else:
            q = np.array([number(v) for v in picked])
            if q.shape!=(3,):raise ValueError('FINITE_XYZ_REQUIRED')
            delta=np.diff(points,axis=0);den=np.einsum('ij,ij->i',delta,delta)
            t=np.clip(np.divide(np.einsum('ij,ij->i',q-points[:-1],delta),den,
                                out=np.zeros_like(den),where=den>0),0,1)
            projected=points[:-1]+t[:,None]*delta
            at=int(np.argmin(np.linalg.norm(projected-q,axis=1)))
            location=float(arc[at]+t[at]*(arc[at+1]-arc[at]))
        point, _, _, _ = self.at(name, location)
        self.initialize_solid()
        from root_native_entity_measurement import _tangent, axis_aligned_section
        normal = request.get('normal')
        if normal is None:normal=_tangent(points,arc,location,self.solid.pitch,branch['closed'])
        else:
            normal=np.array([number(v) for v in normal])
            if normal.shape!=(3,) or not np.linalg.norm(normal)>0:raise ValueError('NONZERO_XYZ_NORMAL_REQUIRED')
            normal=normal/np.linalg.norm(normal)
        if normal is None:raise ValueError('UNIQUE_SECTION_DIRECTION_REQUIRED')
        value=axis_aligned_section(self.solid,point,normal,self.slabs)
        if value is None:value=self.engine.section(point,normal)
        radius=value.get('radius')
        if radius is not None:
            witnesses=value['section_centroid_evidence']['_interior_polygon_witnesses']
            if any(self.owner.nearest_owners(p)!={branch['branch_id']} for p in witnesses):
                radius=None;value=dict(value,status='MISSING_MIXED_BRANCH_SECTION')
        branch_radius=radius
        from root_section_semantics import connected_section
        section=connected_section(dict(value,radius=radius))
        radius=section['radius']
        candidate_radius=radius;qualification=None
        if radius is not None:
            # Original-model geometry only: prevent an apparently disk-like
            # local component from hiding a point contact with another cut.
            from root_native_section_reference import NativePlaneReference
            qualification=NativePlaneReference(self.solid).section(point,normal,connected_area=True)
            if qualification.get('radius') is None:
                radius=None;section=dict(section,classification='unqualified_section')
                value=dict(value,status='SECTION_QUALIFICATION_'+qualification['status'])
        scale=self.report.get('scale_mm_per_model_unit') or 1.
        outline=self.section_outline(point,np.asarray(normal)) if radius is not None else []
        return dict(kind='FRESH_ORIGINAL_ENTITY_SECTION',path_id=name,arc=location,
                    requested_position=picked,position=point.tolist(),normal=np.asarray(normal).tolist(),
                    requested_to_axis_distance=None if picked is None else math.dist(picked,point),
                    radius=None if radius is None else radius*scale,
                    diameter=None if radius is None else radius*2*scale,
                    area=None if radius is None else value['area']*scale**2,
                    unit='mm' if self.report.get('scale_mm_per_model_unit') else 'model_units',
                    status=value['status'],outline=outline,outlines=getattr(self,'last_section_outlines',[]) if outline else [],
                    section_classification=section['classification'],holes_filled=False,
                    single_branch_radius=None if branch_radius is None else branch_radius*scale,
                    candidate_connected_radius=None if candidate_radius is None else candidate_radius*scale,
                    original_section_qualification=qualification,
                    display_outline_is_not_radius_estimator=True,stored_samples_interpolated=False,
                    source_cells_sha256=self.report['entity_binding']['cells_sha256'],
                    direction='user supplied' if request.get('normal') is not None else 'local branch tangent',
                    definition=section['definition'])
