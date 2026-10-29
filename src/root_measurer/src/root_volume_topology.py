"""Exact topology of a finite union of closed observed-occupancy voxels.

These descriptors treat unoccupied/unknown cells as empty in the *derived*
field. They do not establish the true topology of an unobserved root.
"""
import itertools
import numpy as np
from scipy import ndimage


def cubical_euler(indices):
    p=np.asarray(indices)
    if p.ndim!=2 or p.shape[1]!=3 or not np.issubdtype(p.dtype,np.integer):
        raise ValueError('INTEGER_CUBICAL_INDICES_REQUIRED')
    if len(np.unique(p,axis=0))!=len(p):raise ValueError('DUPLICATE_CUBICAL_CELL')
    if not len(p):return {'vertices':0,'edges':0,'faces':0,'cubes':0,'euler':0}
    def count(deltas):
        return len(np.unique(np.concatenate([p+np.array(d) for d in deltas]),axis=0))
    vertices=count(itertools.product([0,1],repeat=3));edges=0;faces=0
    for axis in range(3):
        # Positive oriented edge starts; orientation separates cell types.
        edges+=count(d for d in itertools.product([0,1],repeat=3) if d[axis]==0)
        delta=np.zeros(3,dtype=int);delta[axis]=1
        faces+=count([np.zeros(3,dtype=int),delta])
    cubes=len(p)
    return {'vertices':vertices,'edges':edges,'faces':faces,'cubes':cubes,
        'euler':vertices-edges+faces-cubes}


def topology(occupied, observed=None):
    a=np.asarray(occupied)
    if a.ndim!=3 or a.dtype!=np.bool_:raise ValueError('BOOLEAN_VOLUME_REQUIRED')
    if observed is not None:
        o=np.asarray(observed)
        if o.shape!=a.shape or o.dtype!=np.bool_ or np.any(a&~o):
            raise ValueError('VOLUME_EVIDENCE_BOUNDARY')
    cells=np.argwhere(a)
    counts=cubical_euler(cells)
    foreground,b0=ndimage.label(a,structure=np.ones((3,3,3),bool))
    labels=foreground[tuple(cells.T)] if len(cells) else np.empty(0,np.int32)
    ids,sizes=np.unique(labels,return_counts=True)
    del foreground
    # Padding is mathematically exterior, not added occupied evidence.
    background,bg_count=ndimage.label(np.pad(~a,1,constant_values=True),
        structure=ndimage.generate_binary_structure(3,1))
    exterior=int(background[0,0,0]);assert exterior>0
    b2=int(bg_count-1)
    cavity_ids=[i for i in range(1,bg_count+1) if i!=exterior]
    cavity_sizes=[]
    for i in cavity_ids:cavity_sizes.append(int(np.count_nonzero(background==i)))
    del background
    b1=int(b0+b2-counts['euler'])
    if b1<0:raise AssertionError('NEGATIVE_VOLUME_BETTI1')
    return dict(counts,foreground_components=int(b0),enclosed_empty_components=b2,
        tunnel_dimension=b1,enclosed_empty_sizes=cavity_sizes,
        component_occupied_counts={str(int(i)):int(n) for i,n in zip(ids,sizes)},
        unknown_as_empty_for_derived_topology=True,physical_topology_qualified=False)


def qualify():
    fixtures=[]
    def add(name,a,expected):fixtures.append((name,a,expected))
    a=np.zeros((11,11,11),bool);a[3:8,3:8,3:8]=1;add('solid_cube',a,(1,0,0))
    a=np.zeros((11,11,11),bool);a[3:8,3:8,5]=1;a[4:7,4:7,5]=0;add('true_ring',a,(1,1,0))
    a=np.zeros((11,11,11),bool);a[3:8,3:8,3:8]=1;a[4:7,4:7,4:7]=0;add('closed_cavity',a,(1,0,1))
    a=np.zeros((11,11,11),bool);a[5,5,2:9]=1;a[2:9,5,5]=1;add('solid_junction',a,(1,0,0))
    a=np.zeros((11,11,11),bool);a[2:4,2:4,2:4]=1;a[7:9,7:9,7:9]=1;add('separate',a,(2,0,0))
    a=np.zeros((11,11,11),bool);a[4,4,4]=1;a[5,5,5]=1;add('corner_contact',a,(1,0,0))
    a=np.zeros((11,11,11),bool);a[5,5,5]=1;add('singleton',a,(1,0,0))
    a=np.zeros((11,11,11),bool);add('empty',a,(0,0,0))
    rows=[]
    for name,volume,expected in fixtures:
        for axes in itertools.permutations(range(3)):
            for signs in itertools.product([-1,1],repeat=3):
                v=np.transpose(volume,axes)[tuple(slice(None,None,s) for s in signs)].copy()
                before=v.copy();q=topology(v,np.ones_like(v))
                got=(q['foreground_components'],q['tunnel_dimension'],q['enclosed_empty_components'])
                assert got==expected and np.array_equal(v,before),(name,got,expected)
                rows.append({'fixture':name,'axes':axes,'reflection':signs,'Betti':got,'arrays_unchanged':True})
    invalid=0
    for bad,observed in [(np.zeros((2,2),bool),None),
        (np.zeros((2,2,2),np.uint8),None),
        (np.ones((2,2,2),bool),np.zeros((2,2,2),bool))]:
        try:topology(bad,observed)
        except ValueError:invalid+=1
        else:raise AssertionError('INVALID_TOPOLOGY_INPUT_ACCEPTED')
    return {'status':'PASS_CUBICAL_VOLUME_TOPOLOGY_CONTROLS','controls':len(rows),
        'invalid_input_gates':invalid,'rows':rows,'GT_opened':False,'physical_root_accuracy_qualified':False}
