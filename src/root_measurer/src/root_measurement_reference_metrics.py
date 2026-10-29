"""Reference-only metrology diagnostics; never imported by prediction cores."""
import math
import xml.etree.ElementTree as ET
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching


def native_segments_from_vtp(text, normalization):
    arrays={e.get('Name'):np.fromstring(e.text or '',sep=' ')
            for e in ET.fromstring(text).iter('DataArray') if e.get('Name')}
    raw=arrays['Coordinates'].reshape(-1,3)
    bounds=np.asarray(normalization['original_bounds'],float)
    scale=float(normalization['normalized_scale'])
    if bounds.shape!=(6,) or not np.isfinite(bounds).all() or not math.isfinite(scale) or scale<=0:
        raise ValueError('REFERENCE_NORMALIZATION')
    center=(bounds[::2]+bounds[1::2])/2
    points=(raw-center)*scale
    offsets=arrays['offsets'].astype(np.int64);connectivity=arrays['connectivity'].astype(np.int64)
    if (not len(offsets) or offsets[-1]!=len(connectivity) or
        not np.all(np.diff(np.r_[0,offsets])==2)):
        raise ValueError('NATIVE_TWO_POINT_SEGMENTS_REQUIRED')
    edges=connectivity.reshape(-1,2)
    if np.any(edges<0) or np.any(edges>=len(points)) or not np.isfinite(points).all():
        raise ValueError('REFERENCE_NATIVE_INDICES')
    radii=arrays['radius']*scale
    if radii.shape!=(len(edges),) or not np.isfinite(radii).all() or np.any(radii<=0):
        raise ValueError('REFERENCE_RADIUS_ATTRIBUTE')
    lengths=np.linalg.norm(points[edges[:,1]]-points[edges[:,0]],axis=1)
    if np.any(lengths<=0):raise ValueError('ZERO_LENGTH_REFERENCE_SEGMENT')
    degree=np.bincount(edges.reshape(-1),minlength=len(points))
    return dict(points=points,edges=edges,radii=radii,lengths=lengths,
                endpoint_ids=np.flatnonzero(degree==1),
                branchpoint_ids=np.flatnonzero(degree>=3),
                total_length=float(lengths.sum()),
                arc_weighted_native_radius=float(np.dot(lengths,radii)/lengths.sum()))


def nearest_segment_attributes(queries, starts, ends, radii, batch=256):
    q=np.asarray(queries,float).reshape(-1,3);a=np.asarray(starts,float).reshape(-1,3)
    b=np.asarray(ends,float).reshape(-1,3);r=np.asarray(radii,float)
    if a.shape!=b.shape or not len(a) or r.shape!=(len(a),) or batch<=0:
        raise ValueError('REFERENCE_SEGMENT_SHAPE')
    if not all(np.isfinite(x).all() for x in [q,a,b,r]) or np.any(r<=0):
        raise ValueError('FINITE_REFERENCE_SEGMENTS_REQUIRED')
    direction=b-a;squared=np.einsum('ij,ij->i',direction,direction)
    if np.any(squared<=0):raise ValueError('ZERO_REFERENCE_SEGMENT')
    distances=[];indices=[];fractions=[]
    for start in range(0,len(q),batch):
        delta=q[start:start+batch,None,:]-a[None,:,:]
        t=np.clip(np.einsum('qij,ij->qi',delta,direction)/squared,0,1)
        residual=delta-t[:,:,None]*direction[None,:,:]
        d2=np.einsum('qij,qij->qi',residual,residual)
        ix=np.argmin(d2,axis=1);row=np.arange(len(ix))
        distances.extend(np.sqrt(d2[row,ix]));indices.extend(ix);fractions.extend(t[row,ix])
    ix=np.asarray(indices,np.int64)
    return dict(distances=np.asarray(distances),segment_indices=ix,
                segment_fractions=np.asarray(fractions),native_radius=r[ix])


def weighted_error_summary(prediction,reference,weights):
    p=np.asarray(prediction,float);r=np.asarray(reference,float);w=np.asarray(weights,float)
    if p.shape!=r.shape or p.shape!=w.shape or not all(np.isfinite(x).all() for x in [p,r,w]) or np.any(r<=0) or np.any(w<0):
        raise ValueError('FINITE_WEIGHTED_MEASUREMENT_VALUES')
    if not len(p) or w.sum()<=0:return dict(status='NO_VALID_WEIGHT',samples=0)
    error=p-r;relative=error/r
    return dict(status='COMPLETE_CONDITIONAL_RADIUS_DISCREPANCY',samples=len(p),weight=float(w.sum()),
                weighted_prediction_mean=float(np.dot(w,p)/w.sum()),
                weighted_reference_mean=float(np.dot(w,r)/w.sum()),
                signed_bias=float(np.dot(w,error)/w.sum()),
                MAE=float(np.dot(w,np.abs(error))/w.sum()),
                RMSE=float(np.sqrt(np.dot(w,error*error)/w.sum())),
                signed_relative_bias=float(np.dot(w,relative)/w.sum()),
                MARE=float(np.dot(w,np.abs(relative))/w.sum()),
                no_reference_distance_filter=True,biological_accuracy_admitted=False)


def endpoint_diagnostics(predicted,reference,tolerances):
    p=np.asarray(predicted,float).reshape(-1,3);r=np.asarray(reference,float).reshape(-1,3)
    if not np.isfinite(p).all() or not np.isfinite(r).all() or any(not math.isfinite(t) or t<=0 for t in tolerances):
        raise ValueError('ENDPOINT_DIAGNOSTIC_INPUT')
    def summary(d):
        return None if d is None or not len(d) else dict(mean=float(np.mean(d)),p50=float(np.quantile(d,.5)),
            p95=float(np.quantile(d,.95)),maximum=float(np.max(d)))
    forward=cKDTree(r).query(p)[0] if len(p) and len(r) else None
    reverse=cKDTree(p).query(r)[0] if len(p) and len(r) else None
    rows=[]
    for tol in tolerances:
        distances=np.linalg.norm(p[:,None,:]-r[None,:,:],axis=2)
        matching=maximum_bipartite_matching(csr_matrix(distances<=tol),perm_type='column')
        count=int(np.count_nonzero(matching>=0))
        rows.append(dict(tolerance=float(tol),one_to_one_matches=count,
            unmatched_predicted=len(p)-count,unmatched_reference=len(r)-count,
            precision=count/len(p) if len(p) else None,recall=count/len(r) if len(r) else None,
            diagnostic_not_acceptance_gate=True))
    return dict(predicted_endpoints=len(p),all_native_degree1_reference_including_collar=len(r),
                prediction_to_reference_distance=summary(forward),reference_to_prediction_distance=summary(reverse),
                one_to_one_diagnostics=rows,biological_tip_identity_qualified=False)
