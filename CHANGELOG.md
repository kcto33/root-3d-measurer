# Changelog

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
