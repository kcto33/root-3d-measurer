# Root 3D Measurer maintenance log

## 2026-11-26 Review repairs

Plan: from commit `087b672`, preserve TIFF volume axes, move volume import into the existing supervised measurement worker, and include section subprocesses in the existing memory budget and cancellation. This is software maintenance; historical experiments and scientific parameters remain unchanged. Add regression tests and retain test artifacts under `experiments/pr_validation_20261029/results/repair_20261126/`.

Result: declare/use tifffile directly; pass raw volume settings from API/web to the existing supervised worker. Its existing memory polling now includes descendants, and cancellation/budget stops include section processes. No additional watchdog or scientific estimator was introduced. The separate `import-volume` conversion command retains its existing unsupervised semantics.

Validation: all 368 tests passed in 41.25 s with `python -B -m pytest --basetemp build/review_20261126/run001 -q --tb=short`. This includes eight asymmetric TIFF axis cases, real two-process volume measurement, real low-budget rejection, API/HTTP deferred-import dispatch and process-tree accounting/termination controls. Initial focused suite: 19 passed. An earlier full run had 366 passes and one Windows long-path fixture error; the complete rerun used the shorter build path and passed. Both earlier experiment test artifacts and the shorter build artifacts are retained; no existing data was modified.

Changed files: `src/root_measurer/api.py`, `src/root_measurer/scripts/measure_model_delivery.py`, `src/root_measurer/src/root_model_delivery_service.py`, `src/root_measurer/src/root_voxel_import.py`, `tests/test_voxel_volume_import.py`, `pyproject.toml`, `README.md`, `docs/EXPERIMENT_LOG.md`.

## 2026-11-26 Measurement workflow improvements

Plan: improve only the measurer: label selection, volume metadata and explicit isotropic resampling; mixed native/volume/mesh batches with per-sample options, full parameter/status exports and restored batch history; complete the existing root-axis review roundtrip and verify explicit partial results on complex controls. Preserve source data, original automatic graphs, numerical estimators and unresolved statuses. No generator changes or new scientific experiment. Retain tests under newly allocated `build/functional_20261126/` runs.

Result: label selection, XYZ spacing/units/axis metadata from NPZ, OME-TIFF or explicit JSON, and explicit nearest-neighbor isotropic resampling are available through API/CLI/web. Source files remain unchanged; source/derived volumes are recorded. Calibrated imported packages retain their units on later native measurement; scanner origin is recorded separately from local measurement coordinates.

Result: mixed batches accept per-sample options and continue after failed inputs. Complete parameter/status/coverage rows, cancellations and original denominators are exported; history and unfinished progress are recoverable for viewing without automatically restarting jobs. Root-axis review now shares a checklist between UI, JSON and CSV; original graphs and unavailable values remain intact.

Validation: 380 tests passed in 47.34 s with `python -B -m pytest --basetemp build/functional_20261126/run005 -q --tb=short`; Python and JavaScript syntax/diff checks passed. The prior full run had 379 passes and one Windows TCP abort on an unauthorized POST; the bounded body is now drained before returning 403, and repeated-denial plus full regression passed. No additional watchdog or numerical retuning.

UI validation: actual local volume measurement with pitch 0.1 mm and label 1 gave 2.1 mm axis length and 0.189 mm^3 volume; root-base selection, confirmation, save and restore retained the reviewed axis and zero pending review counts. A two-input batch displayed one successful axis and one failed input, with all-parameter CSV available; no browser console errors. The temporary tab and local server were closed. Artifacts remain in `build/functional_20261126/run001` through `run005` and `ui001`; source datasets and generator files were not modified this round.

Changed files: `README.md`, `docs/input_format.md`, `docs/phenotyping.md`, `docs/roadmap.md`, `docs/EXPERIMENT_LOG.md`; `src/root_measurer/{__init__,api,cli}.py`; `scripts/{measure_model_batch,measure_model_delivery}.py` inside the package; packaged `src/{root_measurement_inputs,root_voxel_import,root_native_entity_measurement,root_model_parameters,root_model_phenotyping,root_model_delivery_service}.py`, `model_delivery_workspace.html`, `model_delivery_phenotyping.js`; `tests/{test_voxel_volume_import,test_batch_workflow,test_phenotyping_integration}.py`. Earlier review repairs and dependency declaration remain in the workspace.
