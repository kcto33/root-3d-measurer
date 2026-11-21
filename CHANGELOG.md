# Changelog

## 0.3.2 — 2026-11-14

- Added occupied-neighbor storage for large connected native solids, preserving Lee0.26 directional/candidate/recheck order and the existing whole-missing-component recovery.
- Retained dense storage for admitted small components, unchanged native pitch/material, graph, endpoint and section measurement definitions and parameters.
- Exact global sparse cavity accounting now uses actual occupied/interval storage limits instead of rejecting an unallocated empty component box.
- Added raw/guarded estimator, axis-order/full-kernel equality and large-empty-domain checks;348 source tests pass. Five-root execution/accuracy is recorded separately in experiments/sparse_thinning_20261114.

## 0.3.1 — 2026-10-29

- Recovered root-relative bridge-fragment analysis in cyclic components with an explicit root, preserving every original cyclic edge and the exact unique-edge length ledger.
- Retained unresolved cyclic cores as geometry; scoped parent/child angles, caliber and taper to their local fragment, withholding cross-core global orders and full-mother density at cyclic boundaries.
- Added local/global scope to JSON/CSV and result-page labels. Pure rings, competing bases, unresolved cores and global confirmation remain explicit.
- Phenotyping method identity is now root-phenotyping-v002; annotation input format remains compatible.

## 0.3.0 — 2026-10-29

- Added arc-weighted diameter classes and diameter-by-depth estimates with qualified single-branch radius coverage and explicit missing lengths.
- Added exact unique-segment inclination, undirected azimuth, position bearing, depth-by-direction/bearing and horizontal radial length distributions; position-based distributions require an explicit reference point.
- Added measurement-bound root-base, continuation, root-type and endpoint annotations, candidate/reviewed continuous axes, developmental order and per-edge Strahler order. Cycles, disconnected components and unresolved continuations remain explicit.
- Added lateral counts, whole-axis/branched-zone density, interbranch distances, unbranched zones, local and whole-axis PCA angles, diameter profiles, taper and offset-window child/parent diameter ratios.
- Added the `phenotype` CLI/Python API and local browser panel with node picking, axis highlighting, annotation import, confirmation, snapshot saving/restoration, JSON and five CSV exports. Original measurements and materials remain hash-verified and unchanged.
- Added analytic, missing-data, subdivision, annotation and public/HTTP integration controls. Human identity confirmation is separate from measurement accuracy; mass, soil, time-series and maize-whorl inference remain outside this release.

## 0.2.0 — 2026-10-29

- Added 17 original-solid/centerline morphology descriptors and exact segment/cube depth profiles, with configurable vertical axis/downward sign (default Z negative).
- Added continuous-segment 3D finite-scale box counting over a fixed scale range and eight grid phases; retained counts, R², grid phase slopes and unavailable states.
- Added result-page orientation controls and morphology/depth-profile exports, parameter integration, batch summary columns and an additive `morphology` CLI/Python API for SHA-bound historical measurements.
- Kept partial-axis, mesh-volume, missing-scale and biological interpretation limits explicit. Axis extraction, length, radius and accuracy acceptance algorithms are unchanged.

## 0.1.0 — 2026-10-29

- Extracted the actual dependency closure of the existing measurer into an installable standalone tool.
- Added a public CLI/Python API, external writable workspace, generated demonstrations and packaged browser assets.
- Retained original finite-solid measurements, explicit user root intervals, area-equivalent connected sections, SHA-bound evidence and strict failure states.
- Declared Shapely as an installed dependency and added cross-platform worker memory monitoring.
- Kept historical precision evidence and unresolved PLY/contact cases explicitly scoped.

## Unreleased — 2026-10-30 memory repair

- Preserve guarded Lee and topology exactly while processing disconnected solids in separate bounded domains; retain every source cell and the350m per-component guard.
-301 source tests and four real formerly guarded inputs verified; details in experiments/memory_scaling_20261030/analysis/memory_scaling_report_run001.md. No wheel published.
