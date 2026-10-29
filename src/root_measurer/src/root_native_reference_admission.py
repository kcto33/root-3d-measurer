"""Scorer-only topology compatibility of a reference graph and finite input."""
import math
import numpy as np
from root_volume_topology import topology
from root_graph_measurement_ledger import graph_measurement_ledger


def admit_reference(solid,reference,max_dense_cells=350000000):
    lo=solid.cells.min(axis=0)-1;shape=solid.cells.max(axis=0)-lo+2
    graph=graph_measurement_ledger(reference['points'],reference['adjacency'],coordinate_units='reference_units')
    expected=[graph['source_components'],graph['source_cycle_rank'],0]
    if math.prod(map(int,shape+2))>max_dense_cells:
        return dict(qualified=False,status='REFERENCE_TOPOLOGY_RESOURCE_UNRESOLVED',expected_Betti=expected)
    volume=np.zeros(tuple(shape),bool);volume[tuple((solid.cells-lo).T)]=True
    material=topology(volume)
    actual=[material['foreground_components'],material['tunnel_dimension'],material['enclosed_empty_components']]
    return dict(qualified=actual==expected,
                status='REFERENCE_TOPOLOGY_COMPATIBLE' if actual==expected else 'REFERENCE_NOT_TOPOLOGICALLY_COMPARABLE',
                actual_Betti=actual,expected_Betti=expected,material_topology=material,
                input_geometry_changed=False,reference_graph_changed=False,
                topology_compatibility_is_necessary_not_sufficient_for_axis_accuracy=True)
