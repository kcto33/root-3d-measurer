"""Scorer-only fixed reference bins and independent finite-entity section references."""
import argparse
import hashlib
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

PROJECT=Path(__file__).resolve().parents[1];WORKSPACE=PROJECT.parents[1]
sys.path.insert(0,str(PROJECT/'src'))
from root_graph_measurement_ledger import graph_measurement_ledger
from root_native_entity_measurement import load_entity
from root_model_delivery_exports import measurement_tables,write_json,file_reference


def polyline(points):
    a=np.asarray(points,float);d=np.diff(a,axis=0);lengths=np.linalg.norm(d,axis=1)
    return a,d,lengths,np.r_[0,np.cumsum(lengths)]


def nearest(point,curve):
    p,d,lengths,arc=curve
    squared=lengths**2;t=np.zeros(len(d));positive=squared>0
    t[positive]=np.clip(np.einsum('ij,ij->i',point-p[:-1],d)[positive]/squared[positive],0,1)
    q=p[:-1]+t[:,None]*d;dist=np.linalg.norm(q-point,axis=1);i=int(np.argmin(dist))
    return q[i],float(dist[i]),i,float(arc[i]+t[i]*lengths[i])


def at_arc(curve,s):
    p,d,lengths,arc=curve;i=min(len(d)-1,int(np.searchsorted(arc,s,side='right')-1));i=max(0,i)
    t=(s-arc[i])/lengths[i] if lengths[i]>0 else 0
    return p[i]+t*d[i],d[i]/lengths[i] if lengths[i]>0 else np.zeros(3)


def prepare_contract(pred,ref):
    if (pred.get('radius_curve') or {}).get('sampling')!='edge_midpoint':
        raise ValueError('SCORING_REQUIRES_FINAL_EDGE_MIDPOINT_SAMPLES')
    pg=pred['candidate_centerline'];pp=np.asarray(pg['points'],float)
    rl=graph_measurement_ledger(ref['points'],ref['adjacency'],coordinate_units=pred['coordinate_units'])
    pl=graph_measurement_ledger(pp,pg['adjacency'],coordinate_units=pred['coordinate_units'])
    rpaths=rl['branches'];ppaths=pl['branches']
    rc=[polyline(np.asarray(ref['points'])[p['source_node_ids']]) for p in rpaths]
    pc=[polyline(pp[p['source_node_ids']]) for p in ppaths]
    costs=np.full((len(rpaths),len(ppaths)),1e6);reverse={}
    for i,r in enumerate(rpaths):
        for j,p in enumerate(ppaths):
            if r['closed']!=p['closed']:continue
            if not r['closed']:
                rd=sorted(len(ref['adjacency'][n]) for n in [r['source_node_ids'][0],r['source_node_ids'][-1]])
                pd=sorted(len(pg['adjacency'][n]) for n in [p['source_node_ids'][0],p['source_node_ids'][-1]])
                if rd!=pd:continue
                ra,rb=rc[i][0][0],rc[i][0][-1];pa,pb=pc[j][0][0],pc[j][0][-1]
                direct=math.dist(ra,pa)+math.dist(rb,pb);flipped=math.dist(ra,pb)+math.dist(rb,pa)
                reverse[i,j]=flipped<direct;costs[i,j]=min(direct,flipped)
            else:
                costs[i,j]=sum(nearest(at_arc(rc[i],s)[0],pc[j])[1] for s in np.linspace(0,rc[i][-1][-1],21))/21
    mapping={int(i):int(j) for i,j in zip(*linear_sum_assignment(costs)) if costs[i,j]<1e5}
    endpoint_map={};back_map={};consistent=True
    for i,j in mapping.items():
        if rpaths[i]['closed']:continue
        rn=[rpaths[i]['source_node_ids'][0],rpaths[i]['source_node_ids'][-1]]
        pn=[ppaths[j]['source_node_ids'][0],ppaths[j]['source_node_ids'][-1]]
        if reverse.get((i,j)):pn=pn[::-1]
        for a,b in zip(rn,pn):
            if (a in endpoint_map and endpoint_map[a]!=b) or (b in back_map and back_map[b]!=a):consistent=False
            endpoint_map[a]=b;back_map[b]=a
    rows=pred['radius_curve']['rows'];by_edge={tuple(sorted(r['node_ids'])):i for i,r in enumerate(rows)}
    radius_by_edge={tuple(sorted(e)):radius for e,radius in zip(ref['edges'],ref['construction_segment_radii'])}
    bins=[];queries={};per_path=[];h=ref['pitch']
    for i,r in enumerate(rpaths):
        total=r['source_length'];count=max(1,math.ceil(total/(2*h)));weight=total/count
        j=mapping.get(i)
        per_path.append(dict(reference_path=i,predicted_path=j,reference_length=total,
             predicted_length=None if j is None else ppaths[j]['source_length'],
             relative_error=None if j is None else abs(ppaths[j]['source_length']-total)/total))
        for k in range(count):
            s=(k+.5)*weight;q,n=at_arc(rc[i],s)
            cell=min(len(r['source_edge_ids'])-1,int(np.searchsorted(rc[i][-1],s,side='right')-1))
            radius=radius_by_edge.get(tuple(sorted(r['source_edge_ids'][max(0,cell)])),ref['max_radius'])
            tolerance=max(math.sqrt(3)*h,.25*radius)
            item=dict(reference_path=i,reference_arc=s,weight=weight,axis_covered=False,
                      position=q.tolist(),normal=n.tolist(),distance_tolerance=tolerance)
            if j is not None:
                projected,distance,edge,position=nearest(q,pc[j])
                ids=ppaths[j]['source_node_ids'];row_id=by_edge.get(tuple(sorted((ids[edge],ids[edge+1]))))
                if row_id is not None:
                    row=rows[row_id];normal=np.asarray(row.get('normal') or [0,0,0]);norm=math.hypot(*normal)
                    angle=abs(float(normal@n))/norm if norm else 0
                    if distance<=tolerance and angle>=math.cos(math.pi/4)-1e-12:
                        item.update(axis_covered=True,prediction_row=row_id,distance=distance)
                        ref_point,ref_distance,at,_=nearest(np.asarray(row['position']),rc[i])
                        ref_normal=rc[i][1][at]/rc[i][2][at]
                        if ref_distance<=tolerance:
                            queries[row_id]=dict(prediction_row=row_id,position=ref_point.tolist(),
                                                 normal=ref_normal.tolist(),reference_path=i)
                            item['reference_query']=row_id
            bins.append(item)
    topology=dict(reference_components=rl['source_components'],predicted_components=pl['source_components'],
                  reference_cycles=rl['source_cycle_rank'],predicted_cycles=pl['source_cycle_rank'],
                  reference_paths=len(rpaths),predicted_paths=len(ppaths),
                  all_paths_matched=len(mapping)==len(rpaths)==len(ppaths),
                  endpoint_junction_correspondence_consistent=consistent)
    return dict(bins=bins,queries=list(queries.values()),
                same_plane_queries=[dict(prediction_row=i,position=rows[i]['position'],normal=rows[i]['normal'])
                                    for i in queries if rows[i].get('normal') is not None],
                per_path=per_path,topology=topology,
                reference_total_length=rl['measurement_graph_length'],
                predicted_total_length=pl['measurement_graph_length'],
                selection_uses_radius_success_or_error=False,
                coverage_bin_step='at most two native pitches; fixed reference arc denominator',
                radius_comparison='prediction station projected to reference axis; reference tangent and same finite entity',
                distance_tolerance='max(sqrt(3)*native_pitch, .25*construction_local_radius)',
                tangent_angle_max_degrees=45)


class SectionOracle:
    def __init__(self,solid,max_radius):
        self.solid=solid;self.occupied=set(map(tuple,solid.cells))
        self.centers=(solid.cells+.5)*solid.pitch;self.tree=cKDTree(self.centers)
        self.bound=max(3*max_radius,6*solid.pitch);self.cache={}

    def axis_reference(self,point,normal):
        active=np.flatnonzero(normal)
        if len(active)!=1:return None
        axis=int(active[0]);other=[i for i in range(3) if i!=axis]
        q=[Fraction(float(x))/Fraction(self.solid.pitch) for x in point]
        at=q[axis];lo=at.numerator//at.denominator
        levels=[lo-1,lo] if at.denominator==1 else [lo]
        chosen=self.solid.cells[np.isin(self.solid.cells[:,axis],levels)][:,other]
        if not len(chosen):return dict(radius=None,status='REFERENCE_NO_SECTION')
        minimum=chosen.min(axis=0)-1;shape=chosen.max(axis=0)-minimum+2
        image=np.zeros(tuple(shape),bool);image[tuple((chosen-minimum).T)]=True
        labels,_=ndimage.label(image,structure=ndimage.generate_binary_structure(2,1))
        incident=[]
        for v in [q[a] for a in other]:
            n=v.numerator//v.denominator;incident.append([n-1,n] if v.denominator==1 else [n])
        values=[]
        for a in incident[0]:
            for b in incident[1]:
                x,y=np.array([a,b])-minimum
                values.append(int(labels[x,y]) if 0<=x<shape[0] and 0<=y<shape[1] else 0)
        if 0 in values or len(set(values))!=1:return dict(radius=None,status='REFERENCE_QUERY_NOT_STRICTLY_INSIDE_ONE_COMPONENT')
        mask=labels==values[0]
        if np.any(ndimage.binary_fill_holes(mask)!=mask):return dict(radius=None,status='REFERENCE_NONDISK_HOLE')
        diagonal=(mask[:-1,:-1]&mask[1:,1:]&~mask[1:,:-1]&~mask[:-1,1:])|(
                  mask[1:,:-1]&mask[:-1,1:]&~mask[:-1,:-1]&~mask[1:,1:])
        if diagonal.any():return dict(radius=None,status='REFERENCE_POINT_CONTACT')
        return dict(radius=math.sqrt(int(mask.sum())*self.solid.pitch**2/math.pi),
                    status='REFERENCE_DENSE_2D_SQUARE_UNION',independent_of_predictor_disk_kernel=True)

    def section(self,point,normal):
        point=np.asarray(point,float);normal=np.asarray(normal,float)
        axis=self.axis_reference(point,normal)
        if axis is not None:return axis
        from root_triangle_measurement_input import TriangleInput
        from root_disconnected_cut_reference import DisconnectedCutSectionReference
        h=self.solid.pitch;bucket=tuple(np.floor(point/(self.bound/2)).astype(int))
        if bucket not in self.cache:
            center=(np.asarray(bucket)+.5)*(self.bound/2)
            ids=self.tree.query_ball_point(center,math.sqrt(3)*(self.bound+h))
            selected={tuple(map(int,c)) for c in self.solid.cells[ids]
                      if np.all(np.abs((c+.5)*h-center)<=self.bound)}
            if len(selected)>10000:return dict(radius=None,status='REFERENCE_LOCAL_COMPLEXITY_LIMIT')
            vertices=[];faces=[];lookup={};artificial=[]
            for cell in sorted(selected):
                for a in range(3):
                    other=[i for i in range(3) if i!=a]
                    for sign in [-1,1]:
                        neighbor=list(cell);neighbor[a]+=sign;neighbor=tuple(neighbor)
                        if neighbor in selected:continue
                        quad=[]
                        for u,v in [(0,0),(1,0),(1,1),(0,1)]:
                            key=list(cell);key[a]+=int(sign>0);key[other[0]]+=u;key[other[1]]+=v;key=tuple(key)
                            if key not in lookup:lookup[key]=len(vertices);vertices.append(key)
                            quad.append(lookup[key])
                        if sign!=(1 if a in [0,2] else -1):quad=quad[::-1]
                        faces.extend([[quad[0],quad[1],quad[2]],[quad[0],quad[2],quad[3]]])
                        if neighbor in self.occupied:artificial.append([vertices[i] for i in quad])
            v=np.asarray(vertices,float).reshape(-1,3)*h;f=np.asarray(faces,np.int64).reshape(-1,3)
            v.setflags(write=False);f.setflags(write=False)
            digest=hashlib.sha256(v.tobytes()+f.tobytes()).hexdigest()
            try:engine=DisconnectedCutSectionReference(TriangleInput(v,f,digest,dict(reference_subset=True)))
            except (ValueError,AssertionError):engine=None
            if len(self.cache)>=16:self.cache.pop(next(iter(self.cache)))
            self.cache[bucket]=(engine,np.asarray(artificial,float).reshape(-1,4,3)*h)
        engine,artificial=self.cache[bucket]
        if engine is None:return dict(radius=None,status='REFERENCE_BOUNDARY_NOT_QUALIFIED')
        if len(artificial):
            values=(artificial-point)@normal
            if np.any((values.min(axis=1)<=0)&(values.max(axis=1)>=0)):
                return dict(radius=None,status='REFERENCE_SECTION_TOUCHES_SUBSET_BOUNDARY')
        try:
            answer=engine.section(point,normal)
            return dict(radius=answer.get('radius'),status=answer['status'],
                        independent_of_predictor_disk_kernel=True)
        except (ValueError,AssertionError,OverflowError) as exc:
            return dict(radius=None,status='REFERENCE_UNQUALIFIED',reason=str(exc))


def score(measurement,reference,output,previous_score=None,reference_engine='native_plane_geos_owned_v003'):
    start=time.perf_counter();measurement=Path(measurement);reference=Path(reference);output=Path(output)
    pred=json.loads(measurement.read_text(encoding='utf-8'));ref=json.loads(reference.read_text(encoding='utf-8'))
    whole_component=reference_engine=='native_plane_geos_connected_v004'
    if whole_component:
        from root_section_semantics import connected_section
        for row in pred['radius_curve']['rows']:
            value=connected_section(row);row['radius']=value['radius'];row['available']=value['available']
    contract=prepare_contract(pred,ref)
    contract.update(measurement=file_reference(measurement),reference=file_reference(reference),
                    scorer=file_reference(Path(__file__).resolve(),WORKSPACE),reference_engine=reference_engine,
                    reference_implementation=file_reference(PROJECT/'src/root_native_section_reference.py',WORKSPACE),
                    reference_admission_implementation=file_reference(PROJECT/'src/root_native_reference_admission.py',WORKSPACE),
                    reference_ownership_implementation=file_reference(PROJECT/'src/root_branch_section_measurement.py',WORKSPACE),
                    radius_reference_eligibility=('entire original connected material section, holes deducted; no branch ownership exclusion' if whole_component else 'strict disk and no foreign reference branch at independent cube-cut centroids; full reference arc denominator'))
    write_json(output.with_name(output.stem+'_queries.json'),contract)
    solid,binding,_=load_entity(reference.parent)
    if binding['cells_sha256']!=pred['entity_binding']['cells_sha256']:raise ValueError('REFERENCE_AND_PREDICTION_DIFFERENT_ENTITIES')
    from root_native_reference_admission import admit_reference
    admission=admit_reference(solid,ref)
    write_json(output.with_name(output.stem+'_reference_admission.json'),admission)
    if reference_engine.startswith('native_plane_geos'):
        from root_native_section_reference import NativePlaneReference
        oracle=NativePlaneReference(solid,whole_plane=reference_engine=='native_plane_geos_v001')
    else:oracle=SectionOracle(solid,ref['max_radius'])
    from root_branch_section_measurement import BranchOwnership
    reference_ownership=BranchOwnership(ref['points'],ref['adjacency'])
    answers={};same_plane={};reused=0
    if previous_score is not None:
        old=json.loads(Path(previous_score).read_text(encoding='utf-8'))
        bound=old['contract'];actual=file_reference(bound['path'])
        if any(actual[k]!=bound[k] for k in ['sha256','bytes']):raise ValueError('REFERENCE_CACHE_CONTRACT_CHANGED')
        old_contract=json.loads(Path(bound['path']).read_text(encoding='utf-8'))
        if old_contract.get('reference_engine','triangle_subset_v001')!=reference_engine:
            raise ValueError('REFERENCE_CACHE_ENGINE_MISMATCH')
        for key in ['reference_implementation','reference_admission_implementation','reference_ownership_implementation']:
            if old_contract.get(key)!=contract[key]:raise ValueError('REFERENCE_CACHE_IMPLEMENTATION_MISMATCH')
        for key in ['measurement','reference','queries','bins','per_path','topology']:
            if old_contract[key]!=contract[key]:raise ValueError('REFERENCE_CACHE_QUERY_OR_INPUT_MISMATCH')
        answers={int(k):v for k,v in old['reference_answers'].items()};reused=len(answers)
    for i,query in enumerate(contract['queries']):
        if query['prediction_row'] not in answers:
            answers[query['prediction_row']]=(
                oracle.owned_section(query['position'],query['normal'],reference_ownership,query['reference_path'])
                if reference_engine=='native_plane_geos_owned_v003' else oracle.section(query['position'],query['normal'],**({'connected_area':True} if whole_component else {})))
        if i%50==0:print(json.dumps(dict(reference_queries=i+1,total=len(contract['queries']),seconds=time.perf_counter()-start)),flush=True)
    for i,query in enumerate(contract['same_plane_queries']):
        same_plane[query['prediction_row']]=oracle.section(query['position'],query['normal'],**({'connected_area':True} if whole_component else {}))
        if i%50==0:print(json.dumps(dict(same_plane_queries=i+1,total=len(contract['same_plane_queries']),seconds=time.perf_counter()-start)),flush=True)
    total=contract['reference_total_length'];axis=radius_arc=scorable=paired=errors=0.
    same_paired=same_errors=0.
    rows=pred['radius_curve']['rows']
    for item in contract['bins']:
        if not item['axis_covered']:continue
        w=item['weight'];axis+=w;row=rows[item['prediction_row']]
        same=same_plane.get(item['prediction_row'],{}).get('radius')
        if same is not None and same>0 and row['radius'] is not None:
            same_paired+=w;same_errors+=w*abs(row['radius']-same)/same
        answer=answers.get(item.get('reference_query'),{})
        r=answer.get('radius')
        if r is not None and r>0:
            scorable+=w
            if row['radius'] is not None:
                radius_arc+=w
                paired+=w;errors+=w*abs(row['radius']-r)/r
    coverage=measurement_tables(pred)[2]
    topology=contract['topology']
    topology_pass=topology['all_paths_matched'] and topology['endpoint_junction_correspondence_consistent'] and (
        topology['reference_components']==topology['predicted_components']) and (
        topology['reference_cycles']==topology['predicted_cycles'])
    full_length_error=abs(contract['predicted_total_length']-total)/total
    path_pass=all(p['relative_error'] is not None and p['relative_error']<=.05 for p in contract['per_path'])
    endpoint_reference=all(r['status']=='RESOLVED' for r in ref['endpoints'])
    endpoint_prediction=all(r['status']=='geometric_surface_estimate' for r in pred.get('geometric_endpoints',[]))
    def fraction(value):return min(1.,max(0.,value/total))
    measured_error=errors/paired if paired else None
    metrics=dict(total_length_relative_error=full_length_error,C_axis_ref=fraction(axis),
                 C_radius_ref=fraction(radius_arc),V_radius_out=coverage['V_radius_out'],
                 reference_radius_scorable_fraction=fraction(scorable),
                 radius_error_paired_reference_fraction=fraction(paired),
                 radius_weighted_relative_error=measured_error,
                 same_plane_weighted_relative_error=same_errors/same_paired if same_paired else None,
                 same_plane_paired_reference_fraction=fraction(same_paired),
                 output_total_length=contract['predicted_total_length'],reference_total_length=total)
    correspondence_qualified=topology_pass and admission['qualified']
    if not correspondence_qualified:
        metrics['candidate_C_axis_ref']=metrics['C_axis_ref'];metrics['candidate_C_radius_ref']=metrics['C_radius_ref']
        metrics['C_axis_ref']=None;metrics['C_radius_ref']=None
    checks=dict(reference_comparable=admission['qualified'],topology=topology_pass,reference_endpoints=endpoint_reference,
                predicted_geometric_endpoints=endpoint_prediction,total_length=full_length_error<=.03,
                per_path_length=path_pass,reference_axis_coverage=correspondence_qualified and metrics['C_axis_ref']>=.95,
                reference_radius_coverage=correspondence_qualified and metrics['C_radius_ref']>=.90,
                output_radius_validity=coverage['V_radius_out'] is not None and coverage['V_radius_out']>=.90,
                reference_radius_scorable=metrics['reference_radius_scorable_fraction']>=.90,
                radius_error=measured_error is not None and measured_error<=.05)
    result=dict(schema_version='native-model-reference-score-v003',metrics=metrics,checks=checks,
                precision_pass=all(checks.values()),topology=topology,per_path=contract['per_path'],
                reference_admission=admission,reference_engine=reference_engine,
                coverage_correspondence_qualified=correspondence_qualified,
                reference_answers=answers,same_plane_reference_answers=same_plane,
                contract=file_reference(output.with_name(output.stem+'_queries.json')),
                elapsed_seconds=time.perf_counter()-start,
                missing_branches_or_reference_intervals_removed=False,
                reference_is_actual_finite_entity=True,
                continuous_generator_radius_not_used_as_section_truth=True)
    result['identical_prior_reference_queries_reused']=reused
    result['radius_quantity']='whole_original_connected_section_area_equivalent' if whole_component else 'legacy_single_branch_radius'
    if whole_component:
        result['single_branch_acceptance_replaced']=False
    write_json(output,result)
    print(json.dumps(dict(metrics=metrics,checks=checks,precision_pass=result['precision_pass']),ensure_ascii=False),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--measurement',required=True,type=Path)
    parser.add_argument('--reference',required=True,type=Path);parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--previous-score',type=Path)
    parser.add_argument('--reference-engine',choices=['native_plane_geos_connected_v004','native_plane_geos_owned_v003','native_plane_geos_component_v002','native_plane_geos_v001','triangle_subset_v001'],default='native_plane_geos_owned_v003')
    args=parser.parse_args();score(args.measurement,args.reference,args.output,args.previous_score,args.reference_engine)


if __name__=='__main__':main()
