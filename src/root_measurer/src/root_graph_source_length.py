"""SOURCE-restricted H1 measure of a qualified embedded finite graph.

Consumes immutable proof receipts, not geometry. File identities must be checked
by the caller. An unsaturated COMPLETE ledger is essential: absence in a capped
or partial ledger is not a zero-boundary-length certificate.
"""
import hashlib
import json
from fractions import Fraction as F


def q(p):
    if not isinstance(p, (list, tuple)) or len(p) != 2 or any(not isinstance(s, str) or len(s) > 32772 for s in p):
        raise ValueError('BOUNDED_CANONICAL_RATIONAL_REQUIRED')
    n, d = (int(s, 16) for s in p)
    if d <= 0 or max(abs(n).bit_length(), d.bit_length()) > 131072:
        raise ValueError('BOUNDED_CANONICAL_RATIONAL_REQUIRED')
    x = F(n, d)
    if [hex(x.numerator), hex(x.denominator)] != list(p):
        raise ValueError('CANONICAL_RATIONAL_REQUIRED')
    return x


def pack(x):
    return [hex(x.numerator), hex(x.denominator)]


def bounds(b):
    lo, hi, width = (q(b[k]) for k in ('lower', 'upper', 'width'))
    if not 0 <= lo <= hi or hi - lo != width:
        raise ValueError('ORDERED_EXACT_LENGTH_ENCLOSURE_REQUIRED')
    return lo, hi


def restrict_graph_length(partition, support, complete, embedding, *, support_reference):
    """Prove zero boundary H1, then bound H1(G intersect int R)=H1(G intersect R).

    G must already be embedded, so summing its distinct original edge lengths
    agrees with geometric H1. Finite singleton events and isolated vertices have
    H1 zero. This never implies pointwise containment or biological identity.
    """
    a, ar = partition['case'], partition['result']
    b, br = support['case'], support['result']
    c = complete['result']
    ec, er = embedding['case'], embedding['result']
    for r in (ar, br, c, er):
        if r.get('available') is not True:
            raise ValueError('ALL_QUALIFIED_RECEIPTS_REQUIRED')
    for key in ('root', 'mesh', 'source_frame', 'coordinate_units', 'graph_metadata'):
        if not a.get(key) or a[key] != b.get(key) or a[key] != ec.get(key):
            raise ValueError('EXACT_SOURCE_FRAME_GRAPH_NAMESPACE_REQUIRED')
    if a.get('index') != b.get('index') or a.get('support_reference') != support_reference or complete['input_pair'].get('support') != support_reference:
        raise ValueError('EXACT_INDEX_AND_IMMUTABLE_SUPPORT_REFERENCE_REQUIRED')
    m, x, v = a['graph_metadata'], br['classification'], ar['result']
    bind = a['length_binding']
    frame_sha = hashlib.sha256(json.dumps(a['source_frame'], sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    if (a['source_frame'].get('input_graph') != m['original_graph'] or
            a['index'].get('original_SOURCE_qualified') is not True or
            a['index']['source_binding'] != dict(original_source=a['mesh'], source_frame=a['source_frame'], coordinate_units=a['coordinate_units']) or
            v.get('available') is not True or v.get('binding') != bind or
            bind != dict(graph_sha256=m['original_graph']['sha256'], source_sha256=a['mesh']['sha256'],
                index_sha256=a['index']['sha256'], embedding_receipt_sha256=a['embedding_reference']['sha256'],
                support_receipt_sha256=support_reference['sha256'], source_frame_sha256=frame_sha, coordinate_units=a['coordinate_units']) or
            v.get('canonical_edge_order_sha256') != m['canonical_edge_labels_SHA256'] or
            er.get('declared_graph_embedding') is not True or er.get('original_graph') != m['original_graph'] or
            er['classification'].get('canonical_edge_labels_SHA256') != m['canonical_edge_labels_SHA256'] or
            er['classification'].get('nodes') != m['nodes'] or er['classification'].get('edges') != m['edges'] or
            c.get('root') != a['root'] or c.get('source') != a['mesh'] or c.get('original_graph') != m['original_graph'] or c.get('coordinate_units') != a['coordinate_units']):
        raise ValueError('ADMITTED_SOURCE_GRAPH_EMBEDDING_AND_LENGTH_BINDING_REQUIRED')
    if (c.get('complete_positive_exterior_length_coverage') is not True or
            c.get('singleton_events_have_zero_H1_measure') is not True or c.get('registered_capacity') != 16 or
            not isinstance(c.get('qualified_mask_kernel_sha256'), str) or len(c['qualified_mask_kernel_sha256']) != 64 or
            any(x.get(k) is not True for k in ('details_complete', 'complete_edge_status_bitset', 'pair_partition_complete')) or
            x.get('isolated_nodes_not_interior') != 0):
        raise ValueError('QUALIFIED_COMPLETE_UNSATURATED_EVENT_COVER_REQUIRED')
    n = m['edges']
    if not isinstance(n, int) or isinstance(n, bool) or not 0 <= n <= 2000000 or len(x['details']) > 1000:
        raise ValueError('REGISTERED_DENOMINATOR_GUARD')
    bits = bytes.fromhex(x['canonical_edge_interior_bitset_hex'])
    if len(bits) != (n+7)//8 or hashlib.sha256(bits).hexdigest() != x['canonical_edge_interior_bitset_SHA256'] or (n % 8 and bits[-1] >> (n % 8)):
        raise ValueError('EXACT_COMPLETE_CANONICAL_BITSET_REQUIRED')
    residual = {i for i in range(n) if not bits[i//8] & (1 << (i % 8))}
    if (n != x['edges'] or n != v['edges'] or n != c['graph_edges'] or v['nodes'] != m['nodes'] or
            v['interior_bitset_sha256'] != x['canonical_edge_interior_bitset_SHA256'] or
            v['residual_edges'] != len(residual) or c['residual_edges'] != len(residual) or
            x['new_edges_not_entirely_interior'] != len(residual) or len(x['details']) != len(residual) or
            v['interior_edges'] + len(residual) != n or v.get('sum_partition_exact') is not True):
        raise ValueError('ALL_ORIGINAL_EDGES_EXACTLY_COVERED_REQUIRED')
    seen = set(); boundary = 0; outside_singletons = 0; outside_atoms = 0; maxstored = 0
    for d in x['details']:
        if d.get('kind') != 'edge' or d['id'] not in residual or d['id'] in seen:
            raise ValueError('UNIQUE_RESIDUAL_EDGE_LEDGER_REQUIRED')
        seen.add(d['id']); ws = d['missing_witnesses']; maxstored = max(maxstored, len(ws))
        if not 0 < len(ws) < 16:
            raise ValueError('UNSATURATED_LEDGER_REQUIRED')
        for w in ws:
            mask = w['octant_mask']
            if not isinstance(mask, int) or isinstance(mask, bool) or not 0 <= mask < 255:
                raise ValueError('EXACT_MISSING_EVENT_MASK_REQUIRED')
            if w['kind'] == 'singleton':
                if not 0 <= q(w['parameter']) <= 1:
                    raise ValueError('CLOSED_PARAMETER_DOMAIN_REQUIRED')
                boundary += int(mask != 0); outside_singletons += int(mask == 0)
            elif w['kind'] == 'open_atom':
                if not 0 <= q(w['lower']) < q(w['upper']) <= 1:
                    raise ValueError('POSITIVE_OPEN_PARAMETER_DOMAIN_REQUIRED')
                if mask != 0:
                    raise ValueError('POSITIVE_BOUNDARY_ATOM_PREVENTS_ZERO_BOUNDARY_H1')
                outside_atoms += 1
            else:
                raise ValueError('EXACT_REGISTERED_EVENT_KIND_REQUIRED')
    if seen != residual or maxstored != c['max_stored_missing_witnesses'] or outside_atoms != c['complete_exterior_open_atoms']:
        raise ValueError('SAME_COMPLETE_LEDGER_CERTIFICATE_REQUIRED')
    total_lo, total_hi = bounds(v['declared_graph_length'])
    inside_lo, inside_hi = bounds(v['interior_graph_length'])
    residual_lo, residual_hi = bounds(v['residual_entire_edge_length'])
    out_lo, out_hi = bounds(c['whole_graph_SOURCE_exterior_length'])
    if (total_lo != inside_lo + residual_lo or total_hi != inside_hi + residual_hi or
            out_lo > total_hi or out_lo > residual_hi):
        raise ValueError('CONSISTENT_ORIGINAL_LENGTH_PARTITION_REQUIRED')
    lo, hi = max(F(0), total_lo-out_hi), total_hi-out_lo
    if lo > hi:
        raise ValueError('CONSISTENT_RESTRICTED_LENGTH_REQUIRED')
    interval = dict(lower=pack(lo), upper=pack(hi), width=pack(hi-lo))
    return dict(available=True, status='COMPLETE_SOURCE_RESTRICTED_ORIGINAL_GRAPH_H1_LENGTH', root=a['root'],
        original_graph=m['original_graph'], source=a['mesh'], coordinate_units=a['coordinate_units'],
        canonical_edge_order_sha256=m['canonical_edge_labels_SHA256'], graph_edges=n, residual_edges=len(residual),
        boundary_singleton_event_incidences=boundary, boundary_distinct_point_count=None,
        exterior_singleton_event_incidences=outside_singletons, positive_boundary_open_atoms=0,
        complete_boundary_H1_zero=True, boundary_H1_length=dict(lower=pack(F(0)), upper=pack(F(0)), width=pack(F(0))),
        Source_interior_graph_H1_length=interval, closed_Source_graph_H1_length=dict(interval),
        declared_graph_length=v['declared_graph_length'], complete_exterior_H1_length=c['whole_graph_SOURCE_exterior_length'],
        whole_graph_SOURCE_interior=br.get('wholegraph_SOURCE_interior'),
        whole_graph_closed_SOURCE_containment=False if out_lo > 0 else None,
        old_strict_refusals_unchanged=True, old_whole_inside_edge_provider_unchanged=True,
        new_Source_graph_index_queries=0, new_sqrt_evaluations=0, native_allocation_seconds=0, graph_modified=False,
        physical_length=None, biological_length_lower_bound=None, whole_Root_axis=None, branch_ownership=None, actual_GT_accuracy=None)
