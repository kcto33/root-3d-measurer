"""Lossless cubical entity boundary; singular contacts are exposed, never erased."""
import itertools
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def grid(indices,pitch):
    p=np.asarray(indices)
    if p.ndim!=2 or p.shape[1]!=3 or not len(p) or not np.issubdtype(p.dtype,np.integer) or pitch<=0 or not np.isfinite(pitch):
        raise ValueError('NONEMPTY_INTEGER_ENTITY_SCOPE')
    if len(np.unique(p,axis=0))!=len(p):raise ValueError('DUPLICATE_OCCUPIED_CELL')
    low=p.min(axis=0)-1;shape=p.max(axis=0)-low+2
    keys=np.ravel_multi_index((p-low).T,shape);order=np.argsort(keys)
    return p,shape,keys,order,np.array([shape[1]*shape[2],shape[2],1],np.int64)


def component_labels(indices):
    p,shape,keys,order,strides=grid(indices,1.)
    sorted_keys=keys[order];rows=[];cols=[]
    for delta in itertools.product([-1,0,1],repeat=3):
        if delta==(0,0,0) or delta<(0,0,0):continue
        target=keys+np.dot(delta,strides);where=np.searchsorted(sorted_keys,target)
        valid=where<len(p);valid[valid]=sorted_keys[where[valid]]==target[valid]
        rows.append(np.flatnonzero(valid));cols.append(order[where[valid]])
    row=np.concatenate(rows);col=np.concatenate(cols)
    graph=coo_matrix((np.ones(len(row),np.uint8),(row,col)),shape=(len(p),len(p))).tocsr()
    n,labels=connected_components(graph,directed=False)
    return int(n),labels


def boundary(indices,origin,pitch):
    p,shape,keys,order,strides=grid(indices,pitch)
    origin=np.asarray(origin,float)
    if origin.shape!=(3,) or not np.isfinite(origin).all():raise ValueError('FINITE_ENTITY_ORIGIN')
    corners=[];owners=[]
    for axis in range(3):
        b=(axis+1)%3;c=(axis+2)%3
        for sign in [-1,1]:
            exposed=~np.isin(keys+sign*strides[axis],keys,assume_unique=True)
            use=np.flatnonzero(exposed)
            quad=np.repeat((2*p[use])[:,None],4,axis=1)
            quad[:,:,axis]+=sign
            quad[:,:,b]+=np.array([-1,1,1,-1])
            quad[:,:,c]+=np.array([-1,-1,1,1])
            if sign<0:quad=quad[:,::-1]
            corners.append(quad);owners.append(use)
    integer_corners,inverse=np.unique(np.concatenate(corners).reshape(-1,3),axis=0,return_inverse=True)
    vertices=origin+integer_corners*pitch/2
    faces=inverse.reshape(-1,4);owner=np.concatenate(owners)
    return vertices,faces,owner


def audit(vertices,faces,cell_count,pitch):
    v=np.asarray(vertices,float);f=np.asarray(faces,np.int64)
    triangles=np.concatenate([f[:,[0,1,2]],f[:,[0,2,3]]])
    terms=np.einsum('ij,ij->i',v[triangles[:,0]],np.cross(v[triangles[:,1]],v[triangles[:,2]]))/6
    volume=float(terms.sum());expected=cell_count*pitch**3
    bound=128*np.finfo(float).eps*max(float(abs(terms).sum()),expected)
    if abs(volume-expected)>bound:raise AssertionError('CUBICAL_VOLUME_OR_ORIENTATION_CHANGED')
    edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,3]],f[:,[3,0]]])
    _,inverse,counts=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    orientation=np.bincount(inverse,weights=np.where(edges[:,0]<edges[:,1],1,-1))
    if np.any(orientation!=0):raise AssertionError('OUTWARD_BOUNDARY_ORIENTATION')
    links={}
    for quad in f:
        for i,node in enumerate(quad):links.setdefault(int(node),[]).append((int(quad[(i-1)%4]),int(quad[(i+1)%4])))
    singular=0
    for pairs in links.values():
        adjacency={}
        for a,b in pairs:adjacency.setdefault(a,[]).append(b);adjacency.setdefault(b,[]).append(a)
        if any(len(n)!=2 for n in adjacency.values()):singular+=1;continue
        visited=set();todo=[next(iter(adjacency))]
        while todo:
            node=todo.pop()
            if node in visited:continue
            visited.add(node);todo.extend(adjacency[node])
        if len(visited)!=len(adjacency):singular+=1
    return {'status':'PASS_EXACT_CUBICAL_BOUNDARY','vertices':len(v),'quads':len(f),'triangles':len(triangles),
        'expected_cell_volume':float(expected),'oriented_boundary_volume':volume,'volume_numerical_bound':bound,
        'open_edges':int(np.count_nonzero(counts==1)),'nonmanifold_edges':int(np.count_nonzero(counts>2)),
        'nonmanifold_vertices':singular,'regular_manifold':bool(not singular and np.all(counts==2)),
        'coordinate_units':'normalized_model_coordinates','smooth_physical_surface_claim':False}


def qualify():
    fixtures=[]
    fixtures.append(('singleton',np.array([[0,0,0]]),False))
    fixtures.append(('solid',np.argwhere(np.ones((3,3,3),bool)),False))
    fixtures.append(('separate',np.array([[0,0,0],[4,0,0]]),False))
    a=np.ones((3,3,1),bool);a[1,1]=0;fixtures.append(('ring',np.argwhere(a),False))
    a=np.ones((3,3,3),bool);a[1,1,1]=0;fixtures.append(('cavity',np.argwhere(a),False))
    fixtures.append(('corner_contact',np.array([[0,0,0],[1,1,1]]),True))
    fixtures.append(('edge_contact',np.array([[0,0,0],[1,1,0]]),True))
    fixtures.append(('Y',np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]),False))
    controls=[]
    for name,cells,singular in fixtures:
        for perm in itertools.permutations(range(3)):
            for signs in itertools.product([-1,1],repeat=3):
                p=cells[:,perm]*signs;origin=np.array([.2,-.8,1.4]);pitch=.07
                v,f,owner=boundary(p,origin,pitch);q=audit(v,f,len(p),pitch)
                assert np.array_equal(np.bincount(owner,minlength=len(p))>0,np.ones(len(p),bool)) or name=='solid'
                assert bool(q['nonmanifold_vertices'])==singular and q['open_edges']==0
                assert np.allclose(v.min(0),origin+(p.min(0)-.5)*pitch)
                assert np.allclose(v.max(0),origin+(p.max(0)+.5)*pitch)
                n,l=component_labels(p)
                assert n==(2 if name=='separate' else 1)
                controls.append({'name':name,'perm':list(perm),'signs':list(signs),'status':'PASS'})
    invalid=0
    for p,pitch in [(np.array([[0,0,0],[0,0,0]]),1.),(np.array([[.1,0,0]]),1.),(np.empty((0,3),int),1.),(np.array([[0,0,0]]),0.)]:
        try:boundary(p,[0,0,0],pitch)
        except ValueError:invalid+=1
        else:raise AssertionError('INVALID_SURFACE_INPUT_ACCEPTED')
    return {'status':'PASS_CUBICAL_ENTITY_CONTROLS','controls':controls,'invalid_rejections':invalid,'cells_changed':0}
