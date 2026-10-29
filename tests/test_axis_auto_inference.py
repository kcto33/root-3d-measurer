"""Multi-cue automatic axis inference: synthetic truth and conservative fallbacks."""
import copy
import math

import pytest

from root_axis_architecture import build_axis_architecture


def tree(points, edges):
    adjacency = [[] for _ in points]
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    return dict(points=points, adjacency=adjacency)


def samples(edges, radii):
    return [dict(edge_nodes=[a,b], sampling="edge_midpoint", node_id=None,
                 radius=r, available=r is not None, status="CLOSED_DISK_CONTROL")
            for (a,b),r in zip(edges, radii)]


def run(points, edges, radii, annotations=None):
    graph=tree(points,edges)
    report=build_axis_architecture(graph, annotations=annotations or
        dict(root_bases=[dict(node_id=0,root_type="primary")]),
        radius_samples=samples(edges,radii))
    actual=[tuple(e) for axis in report["axes"] for e in axis["edge_ids"]]
    assert sorted(actual)==sorted(tuple(sorted(e)) for e in edges)
    assert math.isclose(report["edge_partition"]["total_length_model_units"],
                        sum(math.dist(points[a],points[b]) for a,b in edges))
    return report


def main(result):
    return next(a for a in result["axes"] if a["parent_axis_id"] is None)


def test_automatic_continuation_and_lateral_order():
    edges=[(0,1),(1,2),(2,3),(1,4),(4,5),(2,6)]
    points=[[0,0,0],[0,0,1],[0,0,2],[0,0,3],[1,0,1.2],[2,0,1.3],[1,0,2.2]]
    before=copy.deepcopy(points)
    result=run(points,edges,[1,1,.95,.5,.45,.4])
    assert main(result)["node_ids"]==[0,1,2,3]
    assert result["status"]=="AUTO_INFERRED_ARCHITECTURE"
    assert all(a["confirmation_status"]=="AUTO_INFERRED" for a in result["axes"])
    assert [a["developmental_order"] for a in result["axes"]]==[0,1,1]
    assert all(d["confidence_tier"]=="HIGH_GEOMETRIC" for d in result["junction_decisions"])
    assert all(d["biological_identity_validated"] is False for d in result["junction_decisions"])
    assert points==before


def test_radius_overrules_longer_branch_when_direction_agrees():
    edges=[(0,1),(1,2),(2,3),(1,4),(4,5),(5,6)]
    points=[[0,0,0],[0,0,1],[0,0,2],[0,0,3],[1,0,1.3],[2,0,1.5],[3,0,1.7]]
    result=run(points,edges,[1,1,1,.25,.2,.2])
    assert main(result)["node_ids"]==[0,1,2,3]
    decision=result["junction_decisions"][0]
    assert decision["outgoing_node_id"]==2
    assert decision["alternatives"][0]["downstream_max_path_length"] < decision["alternatives"][1]["downstream_max_path_length"]


def test_missing_radius_still_can_suggest_but_never_claims_validation():
    edges=[(0,1),(1,2),(1,3)]
    points=[[0,0,0],[0,0,1],[0,0,3],[1,0,1.4]]
    result=run(points,edges,[None,None,None])
    decision=result["junction_decisions"][0]
    assert all(x["radius_ratio_to_mother"] is None for x in decision["alternatives"])
    assert decision["confidence_tier"] in ("REVIEW_RECOMMENDED","HIGH_GEOMETRIC")
    assert decision["biological_identity_validated"] is False


def test_symmetric_evidence_tie_remains_unresolved():
    edges=[(0,1),(1,2),(1,3)]
    points=[[0,0,0],[0,0,1],[1,0,2],[-1,0,2]]
    result=run(points,edges,[1,1,1])
    assert result["status"]=="UNRESOLVED_ARCHITECTURE"
    assert result["components"][0]["unresolved_junction_node_ids"]==[1]


def test_manual_override_wins_and_confirmation_remains_manual():
    edges=[(0,1),(1,2),(1,3)]
    points=[[0,0,0],[0,0,1],[0,0,2],[1,0,1.2]]
    result=run(points,edges,[1,1,.4],dict(
        root_bases=[dict(node_id=0)],confirmed_components=[0],
        junction_continuations=[dict(node_id=1,incoming_node_id=0,outgoing_node_id=3)]))
    assert result["status"]=="CONFIRMED_ARCHITECTURE"
    assert main(result)["node_ids"]==[0,1,3]
    assert result["junction_decisions"][0]["status"]=="EXPLICIT_CONTINUATION"


def test_invalid_radius_is_excluded_without_fake_measurement():
    edges=[(0,1),(1,2),(1,3)]
    points=[[0,0,0],[0,0,1],[0,0,2],[1,0,1.2]]
    result=run(points,edges,[1,None,.4])
    assert result["junction_decisions"][0]["alternatives"][0]["radius_ratio_to_mother"] is None
