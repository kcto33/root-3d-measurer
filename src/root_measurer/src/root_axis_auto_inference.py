"""Conservative multi-cue geometric continuation inference on an unchanged rooted tree.

This estimates *axis continuity*, not anatomical origin, and never selects a
root base. Scores are uncalibrated ranking utilities, not probabilities.
Missing/invalid radii are excluded rather than interpolated.
"""
import math
from collections import defaultdict

from root_trait_distributions import is_qualified_radius_sample


METHOD = "root-axis-multicue-v001"
WEIGHTS = {"direction": 0.45, "radius": 0.35, "length": 0.15, "hierarchy": 0.05}


def radius_support_by_edge(samples):
    """Index qualified original radius support, preserving sampling positions."""
    result = defaultdict(list)
    for sample in samples:
        edge = sample.get("edge_nodes")
        if not isinstance(edge, (list, tuple)) or len(edge) != 2:
            raise ValueError("RADIUS_SAMPLE_EDGE_REQUIRED")
        a, b = edge
        if (isinstance(a, bool) or isinstance(b, bool) or not isinstance(a, int)
                or not isinstance(b, int) or a == b):
            raise ValueError("INVALID_RADIUS_SAMPLE_EDGE")
        key = tuple(sorted((a, b)))
        sampling = sample.get("sampling", "edge_midpoint")
        node = sample.get("node_id")
        if sampling == "edge_midpoint" and node is None:
            lo, hi = 0.0, 1.0
        elif sampling == "node_halfedge" and node in key:
            lo, hi = (0.0, 0.5) if node == key[0] else (0.5, 1.0)
        else:
            raise ValueError("UNSUPPORTED_AXIS_RADIUS_SAMPLING")
        if is_qualified_radius_sample(sample) and float(sample["radius"]) > 0:
            result[key].append((lo, hi, float(sample["radius"])))
    return result


def _arm(node, first, parent, children, lengths, *, upstream):
    """Walk from a junction to the next junction or endpoint without branching."""
    segments = []
    previous, current = node, first
    traveled = 0.0
    while True:
        key = tuple(sorted((previous, current)))
        step = lengths[key]
        segments.append((previous, current, traveled, traveled + step))
        traveled += step
        if upstream:
            successor = parent[current]
        else:
            successor = children[current][0] if len(children[current]) == 1 else None
        if successor is None:
            return segments, traveled
        previous, current = current, successor


def _tangent(node, segments, points, window, minimum_fraction):
    """Physical-length local chord; unchanged by subdividing straight edges."""
    if not segments:
        return None, 0.0
    total = segments[-1][3]
    used = min(window, total)
    if used <= 0 or used < minimum_fraction * window - 1e-12 * window:
        return None, used
    for a, b, lo, hi in segments:
        if hi <= lo or used > hi + 1e-12 * max(window, hi):
            continue
        t = min(1.0, max(0.0, (used - lo) / (hi - lo)))
        end = points[a] + (points[b] - points[a]) * t
        vector = end - points[node]
        norm = math.sqrt(float(vector @ vector))
        return (vector / norm if norm > 0 else None), used
    return None, used


def _radius_window(segments, support, offset, window, min_coverage, factor):
    if not segments:
        return dict(value=None, coverage=0.0, effective_window=0.0)
    total = segments[-1][3] * factor
    # A branch with a short internode must not be rejected solely because the
    # requested physical exclusion offset is longer than the observed arm.
    start = min(offset, 0.2 * total)
    width = min(window, total - start)
    if width <= 0:
        return dict(value=None, coverage=0.0, effective_window=width)
    numerator, covered = 0.0, 0.0
    for a, b, lo, hi in segments:
        length = (hi - lo) * factor
        if length <= 0:
            continue
        key = tuple(sorted((a, b)))
        forward = a == key[0]
        for frac_lo, frac_hi, radius in support.get(key, ()):
            left = lo * factor + length * (frac_lo if forward else 1 - frac_hi)
            right = lo * factor + length * (frac_hi if forward else 1 - frac_lo)
            overlap = max(0.0, min(right, start + width) - max(left, start))
            numerator += overlap * radius * factor
            covered += overlap
    coverage = min(1.0, covered / width)
    return dict(value=numerator / covered if covered > 0 and coverage >= min_coverage - 1e-12 else None,
                coverage=coverage, effective_window=width,
                offset_used=start)


def _normalized_choice(rows, key, epsilon):
    values = [(row[key], row["outgoing_node_id"]) for row in rows if row[key] is not None]
    if len(values) != len(rows) or len(values) < 2:
        return None
    values.sort(reverse=True)
    return values[0][1] if values[0][0] - values[1][0] >= epsilon else None


def infer_continuation(node, incoming, child_nodes, *, points, parent, children,
                       lengths, node_orders, longest, radius_support, config, factor):
    """Return selected outgoing node, status, ranking evidence and review tier.

    HIGH_GEOMETRIC requires at least two discriminative cues and no strong
    direction-versus-diameter disagreement. A selection is never a biological
    validation, and uncertain cases always remain review candidates.
    """
    w = float(config["tangent_window"])
    fraction = float(config["minimum_window_fraction"])
    upstream, _ = _arm(node, incoming, parent, children, lengths, upstream=True)
    incoming_vector, incoming_used = _tangent(node, upstream, points, w / factor, fraction)
    mother = _radius_window(upstream, radius_support, config["diameter_offset"],
                            config["diameter_window"], config["minimum_radius_coverage"], factor)
    records = []
    for child in child_nodes:
        outgoing, length_to_junction = _arm(node, child, parent, children, lengths, upstream=False)
        direction, direction_used = _tangent(node, outgoing, points, w / factor, fraction)
        cosine = None if incoming_vector is None or direction is None else max(-1.0, min(
            1.0, float((-incoming_vector) @ direction)))
        # Incoming arm coordinates point toward the root, hence the minus sign.
        caliber = _radius_window(outgoing, radius_support, config["diameter_offset"],
                                 config["diameter_window"], config["minimum_radius_coverage"], factor)
        ratio = None if mother["value"] is None or caliber["value"] is None else caliber["value"] / mother["value"]
        similarity = None if ratio is None else math.exp(-abs(math.log(ratio)) / 0.40)
        incident_in = points[node] - points[incoming]
        incident_out = points[child] - points[node]
        norm_in = math.sqrt(float(incident_in @ incident_in))
        norm_out = math.sqrt(float(incident_out @ incident_out))
        incident_cos = None if not norm_in or not norm_out else float(incident_in @ incident_out / norm_in / norm_out)
        records.append(dict(outgoing_node_id=child, alignment_cosine=cosine,
            incident_alignment_cosine=incident_cos,
            incoming_tangent_support=incoming_used * factor,
            outgoing_tangent_support=direction_used * factor,
            mother_radius=mother["value"], child_radius=caliber["value"],
            radius_ratio_to_mother=ratio, radius_similarity=similarity,
            mother_radius_coverage=mother["coverage"], child_radius_coverage=caliber["coverage"],
            downstream_max_path_length=(lengths[tuple(sorted((node, child)))] + longest[child]) * factor,
            downstream_strahler_order=node_orders[child],
            length_to_next_junction=length_to_junction * factor))
    if not records:
        return None, "UNRESOLVED_CONTINUATION", records, "UNRESOLVED", ["NO_OUTGOING_BRANCH"]

    lengths_out = [row["downstream_max_path_length"] for row in records]
    max_length = max(lengths_out)
    orders = [row["downstream_strahler_order"] for row in records]
    lo_order, hi_order = min(orders), max(orders)
    for row in records:
        row["direction_score"] = None if row["alignment_cosine"] is None else (row["alignment_cosine"] + 1.0) / 2.0
        row["length_score"] = row["downstream_max_path_length"] / max_length if max_length > 0 else None
        row["hierarchy_score"] = ((row["downstream_strahler_order"] - lo_order) / (hi_order - lo_order)
                                   if hi_order > lo_order else None)
    common = {"direction": "direction_score", "radius": "radius_similarity",
              "length": "length_score", "hierarchy": "hierarchy_score"}
    enabled = {cue: field for cue, field in common.items()
               if all(row[field] is not None for row in records)}
    divisor = sum(WEIGHTS[cue] for cue in enabled)
    for row in records:
        row["ranking_score"] = (sum(WEIGHTS[cue] * row[field] for cue, field in enabled.items()) / divisor
                                if divisor > 0 else None)
    ranked = sorted((row for row in records if row["ranking_score"] is not None),
                    key=lambda row: (-row["ranking_score"], row["outgoing_node_id"]))
    if len(ranked) != len(records) or len(ranked) < 2 or ranked[0]["ranking_score"] - ranked[1]["ranking_score"] <= 1e-8:
        return None, "UNRESOLVED_CONTINUATION", records, "UNRESOLVED", ["TIED_OR_MISSING_EVIDENCE"]

    best, second = ranked[:2]
    winner = best["outgoing_node_id"]
    votes = {}
    for cue, field in enabled.items():
        thresholds = {"direction": 0.12, "radius": 0.18, "length": 0.15, "hierarchy": 0.5}
        votes[cue] = _normalized_choice(records, field, thresholds[cue])
    supporting = [cue for cue, node_id in votes.items() if node_id == winner]
    conflicting = [cue for cue in ("direction", "radius") if votes.get(cue) not in (None, winner)]
    margin = best["ranking_score"] - second["ranking_score"]
    direction_ok = best["alignment_cosine"] is not None and best["alignment_cosine"] >= -0.1
    evidence_paired = (("direction" in supporting and any(k in supporting for k in ("radius", "length", "hierarchy")))
                       or ("radius" in supporting and "length" in supporting))
    auto = margin >= 0.12 and evidence_paired and not conflicting and direction_ok
    for row in records:
        row["enabled_cues"] = sorted(enabled)
        row["vote_for_best"] = sorted(supporting)
        row["ranking_score_is_probability"] = False
    if auto:
        return winner, "AUTO_INFERRED_CONTINUATION", records, "HIGH_GEOMETRIC", supporting
    return winner, "INFERRED_MULTI_CUE_CANDIDATE", records, "REVIEW_RECOMMENDED", (
        ["DIRECTION_RADIUS_CONFLICT"] if conflicting else
        ["INSUFFICIENT_INDEPENDENT_EVIDENCE_OR_MARGIN"])
