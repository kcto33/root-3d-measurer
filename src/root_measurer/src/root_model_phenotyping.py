"""Additive, measurement-bound root distributions and reviewed axis phenotypes.

The measured graph and original material are never edited. Root identities are
annotations: even a human-confirmed identity is not an accuracy certificate.
"""
import csv
import hashlib
import io
import json
from pathlib import Path

METHOD = 'root-phenotyping-v001'
ANNOTATION_SCHEMA = 'root-phenotype-annotations-v001'
DEFAULTS = {
    'vertical_axis': 'z', 'down_sign': -1,
    'diameter_edges': [0.5, 1.0, 2.0], 'depth_bins': 10,
    'inclination_edges_deg': list(range(0, 91, 15)),
    'azimuth_edges_deg': list(range(0, 181, 15)),
    'bearing_edges_deg': list(range(0, 361, 45)),
    'radial_edges': [1.0, 2.0, 5.0, 10.0], 'origin': None,
    'tangent_window': 1.0, 'diameter_offset': 1.0, 'diameter_window': 1.0,
    'minimum_window_fraction': 0.5, 'minimum_radius_coverage': 0.5,
}
DISTRIBUTION_OPTIONS = {'vertical_axis', 'down_sign', 'diameter_edges', 'depth_bins',
    'depth_edges', 'inclination_edges_deg', 'azimuth_edges_deg', 'bearing_edges_deg', 'radial_edges', 'origin'}
AXIS_OPTIONS = {'vertical_axis', 'down_sign', 'tangent_window', 'diameter_offset',
    'diameter_window', 'minimum_window_fraction', 'minimum_radius_coverage'}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                    allow_nan=False, separators=(',', ':')).encode()).hexdigest()


def _profile(profile, measurement_sha256):
    if profile is None:
        return {}, {}
    if isinstance(profile, (str, Path)):
        profile = json.loads(Path(profile).read_text(encoding='utf-8-sig'))
    if not isinstance(profile, dict) or set(profile) - {
            'schema_version', 'measurement_sha256', 'settings', 'annotations'}:
        raise ValueError('INVALID_PHENOTYPE_ANNOTATION_PROFILE')
    if profile.get('schema_version') != ANNOTATION_SCHEMA:
        raise ValueError('PHENOTYPE_ANNOTATION_SCHEMA_REQUIRED')
    if profile.get('measurement_sha256') != measurement_sha256:
        raise ValueError('ANNOTATION_MEASUREMENT_SHA256_MISMATCH')
    if not isinstance(profile.get('annotations'), dict) or not isinstance(profile.get('settings', {}), dict):
        raise ValueError('ANNOTATION_SETTINGS_AND_AXES_OBJECTS_REQUIRED')
    return profile.get('settings', {}), profile['annotations']


def _check_source(report):
    if report.get('entity_binding'):
        from root_native_entity_measurement import check_entity_unchanged
        check_entity_unchanged(report['entity_binding'])
        return 'ORIGINAL_ENTITY_HASH_VERIFIED'
    if report.get('input_file_binding'):
        from root_model_delivery_exports import file_reference
        binding = report['input_file_binding']
        actual = file_reference(binding['path'])
        if any(actual[k] != binding[k] for k in ('sha256', 'bytes')):
            raise ValueError('ORIGINAL_MESH_CHANGED')
        return 'ORIGINAL_MESH_HASH_VERIFIED'
    return 'NO_ORIGINAL_INPUT_BINDING'


def analyze_measurement(measurement, *, options=None, annotation_profile=None):
    from root_model_delivery_exports import file_reference, measurement_tables
    from root_trait_distributions import build_trait_distributions
    from root_axis_architecture import build_axis_architecture
    from root_axis_traits import build_axis_traits

    measured = Path(measurement).resolve()
    before = file_reference(measured)
    report = json.loads(measured.read_text(encoding='utf-8-sig'))
    inherited, annotations = _profile(annotation_profile, before['sha256'])
    if options is not None and not isinstance(options, dict):
        raise ValueError('PHENOTYPE_OPTIONS_OBJECT_REQUIRED')
    config = dict(DEFAULTS, **inherited)
    config.update(options or {})
    if set(config) - (DISTRIBUTION_OPTIONS | AXIS_OPTIONS):
        raise ValueError('UNKNOWN_PHENOTYPE_OPTION')
    source_status = _check_source(report)
    raw_graph = report.get('candidate_centerline')
    # The legacy table builder expects a nonempty N x 3 array. Empty results
    # remain explicit missing observations without changing that old contract.
    table_report = report if raw_graph and raw_graph.get('points') else dict(report, candidate_centerline=None)
    _, samples, coverage = measurement_tables(table_report)
    coverage = dict(coverage, coordinate_units=report.get('coordinate_units', 'model_units'),
                    length_unit=report.get('coordinate_units', 'model_units'),
                    scale_mm_per_model_unit=report.get('scale_mm_per_model_unit'),
                    scope='original measurement table coverage; length fields remain in model coordinates')
    distributions = build_trait_distributions(report, samples=samples,
        options={k: v for k, v in config.items() if k in DISTRIBUTION_OPTIONS},
        bounds=report.get('model_bounds'))
    graph = dict(report.get('candidate_centerline') or {'points': [], 'adjacency': []},
                 coordinate_units=report.get('coordinate_units', 'model_units'))
    scale = report.get('scale_mm_per_model_unit')
    architecture = build_axis_architecture(graph, annotations=annotations,
                                           scale_mm_per_model_unit=scale)
    traits = build_axis_traits(graph, architecture, samples,
        options={k: v for k, v in config.items() if k in AXIS_OPTIONS},
        scale_mm_per_model_unit=scale)
    profile = dict(schema_version=ANNOTATION_SCHEMA, measurement_sha256=before['sha256'],
                   settings=config, annotations=architecture.get('annotations', annotations))
    _check_source(report)
    if file_reference(measured) != before:
        raise ValueError('MEASUREMENT_CHANGED')
    return dict(schema_version=METHOD, measurement=before,
        source_binding=report.get('entity_binding') or report.get('input_file_binding'),
        source_status=source_status, settings=config, settings_sha256=_digest(config),
        annotation_profile=profile, annotations_sha256=_digest(profile),
        unit=distributions['unit'], partial_axis=distributions['partial_axis'],
        distributions=distributions, architecture=architecture, axis_traits=traits,
        measurement_coverage=coverage, source_geometry_unchanged=True,
        confirmation_is_accuracy_validation=False,
        scope='observed geometry and user-reviewed root axes; no mass, soil, time or crop-specific inference')


def _flatten(value, prefix=''):
    result = {}
    for key, item in value.items():
        name = prefix + str(key)
        if isinstance(item, dict):
            result.update(_flatten(item, name + '_'))
        elif isinstance(item, (list, tuple)):
            result[name] = json.dumps(item, ensure_ascii=False, allow_nan=False, separators=(',', ':'))
        else:
            result[name] = item
    return result


def phenotyping_csv(value, kind='distributions'):
    from root_trait_distributions import flatten_distribution_rows
    traits = value['axis_traits']
    if kind == 'distributions':
        records = flatten_distribution_rows(value['distributions'])
    elif kind == 'axes':
        records = [_flatten({k: v for k, v in row.items() if k not in ('diameter_profile', 'points')})
                   for row in traits['axes']]
    elif kind == 'branches':
        records = [_flatten(row) for row in traits['branch_relations']]
    elif kind == 'profiles':
        records = [dict(axis_id=axis['axis_id'], confirmation_status=axis.get('confirmation_status'),
                        **_flatten(row)) for axis in traits['axes'] for row in axis['diameter_profile']]
    elif kind == 'summaries':
        records = []
        for group, rows in traits['summaries'].items():
            if isinstance(rows, list):
                records.extend(dict(group=group, **_flatten(row)) for row in rows)
    else:
        raise ValueError('UNKNOWN_PHENOTYPE_EXPORT')
    metadata = dict(method=value['schema_version'], measurement_sha256=value['measurement']['sha256'],
                    settings_sha256=value['settings_sha256'], annotations_sha256=value['annotations_sha256'],
                    output_unit=value['unit'], partial_axis=value['partial_axis'],
                    vertical_axis=value['settings']['vertical_axis'], down_sign=value['settings']['down_sign'])
    rows = [dict(metadata, **r) for r in records]
    fields = list(dict.fromkeys([*metadata, *(key for row in rows for key in row)]))
    out = io.StringIO(newline='')
    writer = csv.DictWriter(out, fieldnames=fields, lineterminator='\r\n')
    writer.writeheader()
    for row in rows:
        # Root-type labels are user text, not spreadsheet formulas.
        safe = {k: "'" + v if isinstance(v, str) and v.lstrip().startswith(('=', '+', '-', '@')) else v
                for k, v in row.items()}
        writer.writerow(safe)
    return out.getvalue().encode('utf-8-sig')


def export_phenotyping(value, prefix):
    """Create a new snapshot and portable annotation profile; refuse overwrite."""
    prefix = Path(prefix).resolve()
    kinds = ('distributions', 'axes', 'branches', 'profiles', 'summaries')
    targets = {'phenotyping': prefix.with_name(prefix.name + '_phenotyping.json'),
               'annotations': prefix.with_name(prefix.name + '_annotations.json')}
    targets.update({kind: prefix.with_name(prefix.name + '_' + kind + '.csv') for kind in kinds})
    if any(path.exists() for path in targets.values()):
        raise ValueError('PHENOTYPE_OUTPUT_EXISTS')
    # Serialize all content first; a bad value cannot leave half an analysis.
    raw = {kind: phenotyping_csv(value, kind) for kind in kinds}
    raw['phenotyping'] = (json.dumps(value, ensure_ascii=False, allow_nan=False, indent=2) + '\n').encode()
    raw['annotations'] = (json.dumps(value['annotation_profile'], ensure_ascii=False, allow_nan=False, indent=2) + '\n').encode()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    created = []
    try:
        for kind, path in targets.items():
            with path.open('xb') as stream:
                created.append(path)
                stream.write(raw[kind])
    except (OSError, ValueError):
        for path in created:
            path.unlink(missing_ok=True)
        raise
    return {kind: str(path) for kind, path in targets.items()}


def latest_annotation_profile(measurement):
    """Restore only a profile explicitly saved for these exact measurement bytes."""
    from root_model_delivery_exports import file_reference
    measured = Path(measurement).resolve()
    digest = file_reference(measured)['sha256']
    folder = measured.parent / 'phenotyping'
    for path in sorted(folder.glob('run*/analysis_annotations.json'), reverse=True):
        try:
            profile = json.loads(path.read_text(encoding='utf-8-sig'))
            if isinstance(profile, dict) and profile.get('measurement_sha256') == digest:
                _profile(profile, digest)
                return profile
        except (OSError, ValueError):
            continue
    return None
