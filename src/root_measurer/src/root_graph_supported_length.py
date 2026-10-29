"""Enclosed lengths of an immutable graph edge partition, never a repaired axis.

The caller must bind qualified full-graph embedding and SOURCE support receipts.
Residual length includes entire refused edges, not merely their exterior parts.
"""
import hashlib
import math
import os
import re
from fractions import Fraction as F


def _q(x):
    if isinstance(x, bool):
        raise ValueError('EXACT_COORDINATE_REQUIRED')
    q = F(x)
    if max(abs(q.numerator).bit_length(), q.denominator.bit_length()) > 4096:
        raise ValueError('COORDINATE_4096_BIT_GUARD')
    return q


def _packed(x):
    return [hex(x.numerator), hex(x.denominator)]


def _digest(s):
    return isinstance(s, str) and re.fullmatch('[0-9a-f]{64}', s)


def measure_graph_length_partition(points, edges, interior_bitset_hex, *,
                                   expected_nodes, expected_edges,
                                   expected_edge_order_sha256, expected_bitset_sha256,
                                   binding, expected_binding,
                                   execution_domain='local', max_work=10000000):
    """Every canonical edge once; integer square-root bounds at fixed 96 bits."""
    if (type(expected_nodes) is not int or not 0 < expected_nodes <= 1200000 or
            type(expected_edges) is not int or not 0 <= expected_edges <= 2000000 or
            len(points) != expected_nodes or len(edges) != expected_edges):
        raise ValueError('COMPLETE_GRAPH_DENOMINATORS')
    if execution_domain == 'local':
        if expected_nodes > 256 or expected_edges > 1000:
            raise ValueError('UNCHANGED_LOCAL_GRAPH_GUARDS')
    elif execution_domain != 'registered_slurm_complete_graph' or not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise ValueError('REGISTERED_COMPLETE_GRAPH_SLURM_REQUIRED')
    required = {'graph_sha256', 'source_sha256', 'index_sha256', 'embedding_receipt_sha256',
                'support_receipt_sha256', 'source_frame_sha256', 'coordinate_units'}
    if (not isinstance(binding, dict) or set(binding) != required or binding != expected_binding or
            any(not _digest(binding[k]) for k in required - {'coordinate_units'}) or
            not isinstance(binding['coordinate_units'], str) or not binding['coordinate_units']):
        raise ValueError('IMMUTABLE_SAME_SOURCE_GRAPH_PROOFS_BINDING')
    if not _digest(expected_edge_order_sha256) or not _digest(expected_bitset_sha256):
        raise ValueError('FROZEN_ORDER_AND_BITSET_SHA_REQUIRED')
    if type(max_work) is not int or not 0 <= max_work <= 10000000:
        raise ValueError('WORK_GUARD_ONLY_TIGHTENED')
    work = 0

    def use(n=1):
        nonlocal work
        work += n
        if work > max_work:
            raise ValueError('NO_PARTIAL_LENGTH_AT_WORK_GUARD')

    xyz = []
    for p in points:
        use()
        if len(p) != 3:
            raise ValueError('THREE_COORDINATES_REQUIRED')
        xyz.append(tuple(map(_q, p)))
    try:
        bits = bytes.fromhex(interior_bitset_hex)
    except (ValueError, TypeError):
        raise ValueError('CANONICAL_COMPLETE_BITSET_REQUIRED')
    if (interior_bitset_hex != bits.hex() or len(bits) != (expected_edges + 7) // 8 or
            hashlib.sha256(bits).hexdigest() != expected_bitset_sha256 or
            (expected_edges % 8 and bits[-1] >> (expected_edges % 8))):
        raise ValueError('CANONICAL_COMPLETE_BITSET_REQUIRED')
    h = hashlib.sha256()
    previous = None
    sums = [[0, 0], [0, 0]]  # residual, interior; scaled by 2**96
    counts = [0, 0]
    exact_counts = [0, 0]
    used = set()
    for i, edge in enumerate(edges):
        use(8)
        if (len(edge) != 2 or any(type(j) is not int or not 0 <= j < expected_nodes for j in edge) or
                not edge[0] < edge[1] or previous is not None and tuple(edge) <= previous):
            raise ValueError('COMPLETE_CANONICAL_UNIQUE_EDGE_ORDER')
        a, b = edge
        previous = tuple(edge)
        used.update(edge)
        h.update(('%d:%d\n' % (a, b)).encode())
        d2 = sum((xyz[b][k] - xyz[a][k]) ** 2 for k in range(3))
        if max(abs(d2.numerator).bit_length(), d2.denominator.bit_length()) > 131072:
            raise ValueError('UNCHANGED_INTERNAL_131072_BIT_GUARD')
        if not d2:
            raise ValueError('QUALIFIED_EMBEDDING_CANNOT_HAVE_ZERO_EDGE')
        numerator = d2.numerator << 192
        denominator = d2.denominator
        lo = math.isqrt(numerator // denominator)
        exact = lo * lo * denominator == numerator
        hi = lo if exact else lo + 1
        assert lo * lo * denominator <= numerator <= hi * hi * denominator
        group = (bits[i // 8] >> (i % 8)) & 1
        sums[group][0] += lo
        sums[group][1] += hi
        counts[group] += 1
        exact_counts[group] += int(exact)
    if h.hexdigest() != expected_edge_order_sha256:
        raise ValueError('ORIGINAL_CANONICAL_EDGE_SHA_MISMATCH')
    total = [sums[0][i] + sums[1][i] for i in range(2)]

    def enclosed(v):
        return dict(lower=_packed(F(v[0], 1 << 96)), upper=_packed(F(v[1], 1 << 96)),
                    width=_packed(F(v[1] - v[0], 1 << 96)))

    return dict(available=True, status='COMPLETE_IMMUTABLE_EDGE_LENGTH_PARTITION',
                binding=binding, nodes=expected_nodes, edges=expected_edges,
                canonical_edge_order_sha256=h.hexdigest(), interior_bitset_sha256=expected_bitset_sha256,
                interior_edges=counts[1], residual_edges=counts[0], isolated_nodes=expected_nodes - len(used),
                interior_graph_length=enclosed(sums[1]), residual_entire_edge_length=enclosed(sums[0]),
                declared_graph_length=enclosed(total), exact_sqrt_edges=sum(exact_counts),
                fixed_fraction_bits=96, sum_partition_exact=True, work=work,
                original_graph_modified=False, Source_queries=0, old_embedding_support_kernels_reexecuted=False,
                residual_is_exterior_length=False, physical_length=None, whole_Root_axis=None,
                biological_length_lower_bound=None, branch_ownership=None, actual_GT_accuracy=None)
