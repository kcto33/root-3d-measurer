"""Measure an occupied solid without access to its generating axes or radii."""
import math
import numpy as np
from root_guarded_thinning import guarded_skeletonize
from root_reconstruction_measurement import ReconstructedSolid,neighbor_graph,cubical_radii
from root_lattice_graph import reduce_lattice_graph
from root_centerline_measurement import canonicalize_corner_junctions
from root_graph_measurement_ledger import graph_measurement_ledger


def measure_tube_solid(solid,*,max_dense_cells=64000000):
    if not isinstance(solid,ReconstructedSolid):raise ValueError('RECONSTRUCTED_SOLID_REQUIRED')
    cells=solid.cells;h=solid.pitch
    if len(cells):
        lo=cells.min(axis=0)-1;shape=cells.max(axis=0)-lo+2
        if math.prod(map(int,shape+2))>max_dense_cells:raise ValueError('TUBE_MEASUREMENT_DENSE_GUARD')
        volume=np.zeros(tuple(shape),bool);volume[tuple((cells-lo).T)]=True
        thinning={};out=guarded_skeletonize(volume,receipt=thinning,max_dense_cells=max_dense_cells)
        nodes=np.argwhere(out)+lo
        del volume,out
        pairs,_=neighbor_graph(nodes);adj=[set() for unused in nodes]
        for i,j in pairs:adj[i].add(int(j));adj[j].add(int(i))
        reduced,reduction=reduce_lattice_graph((nodes+.5)*h,adj,h)
        points,graph,representation=canonicalize_corner_junctions(nodes,reduced,h)
        radii=list(map(float,cubical_radii(solid,nodes)))+[None]*(len(points)-len(nodes))
    else:
        nodes=np.empty((0,3),int);points=np.empty((0,3));graph=[];radii=[]
        thinning={};reduction={};representation=dict(status='EMPTY_SOLID')
    ledger=graph_measurement_ledger(points,graph,coordinate_units='caller_solid_coordinate_units',radii=radii)
    weights=np.zeros(len(points))
    for i,j in ledger['source_edge_ids']:
        length=math.dist(points[i],points[j]);weights[i]+=length/2;weights[j]+=length/2
    return dict(status='COMPLETE_UNCALIBRATED_TUBE_GEOMETRIC_MEASUREMENT',
        occupied_cells=len(cells),pitch=h,thinning=thinning,reduction=reduction,representation=representation,
        ledger=ledger,points=points.tolist(),adjacency=[sorted(ns) for ns in graph],
        source_skeleton_cells=nodes.tolist(),cubical_boundary_radii=radii,node_arc_weights=weights.tolist(),
        reference_axes_or_radii_accepted_by_inference=False,material_geometry_modified=False,
        radius_definition='distance from source skeleton cell center to occupied cube-union boundary',
        centerline_candidate_available=representation['status']=='PASS_SUPPORTED_GEOMETRIC_GRAPH_TOPOLOGY',
        physical_or_biological_measurements_admitted=False)
