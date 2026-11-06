# Agent Development Guide

## Project Scope and Isolation

**The project root is `.` relative to this `AGENTS.md`: the directory containing this file. This directory is the exclusive project boundary for every agentic development session.**

- Read, search, analyze, and modify only files contained within this root. Keep commands scoped to it; do not search parent directories, sibling projects, home directories, or other external locations.
- Do not follow symbolic links, Git worktree pointers, submodule references, object alternates, configuration includes, or other references that lead outside the root. Check containment before following paths; a repository-relative path containing `..` may escape the boundary.
- Derive project architecture, requirements, and conventions only from verified contents of this repository. Do not import assumptions from other projects or previous sessions.
- If required information, data, credentials, configuration, or tooling resides outside the repository, report the dependency rather than accessing it. External links in documentation identify dependencies; they do not authorize inspecting external repositories or services.
- Keep development data, temporary files, caches, and generated outputs within the root. Commands documented below describe existing workflows and do not override these restrictions. Report workflows requiring external access or side effects before attempting them.
- These instructions apply throughout the repository. Re-read relevant local files when making changes; distinguish observed behavior from proposed behavior and unverified assumptions.

## Purpose and Repository Map

Arotake is a Python toolkit and CLI for evaluating NOAA HRRR meteorological forecasts against hourly station observations, including LCD v2 and ISD-Lite datasets. It calculates spatial statistics over time and produces NetCDF results and diagnostic PNGs. See [README.md](README.md) for usage, CLI arguments, input file layout, and the public API.

| Path | Responsibility |
| --- | --- |
| `src/arotake/analyze_hrrr_vs_observations.py` | `run_analysis` orchestration, argument parsing, and installed CLI entry point. |
| `src/arotake/observations.py` | `TimeSeries` loads observation NetCDF into xarray and adds timezone-aware `UTC` timestamps. |
| `src/arotake/interpolation.py` | Bilinear xESMF interpolation of model fields to station coordinates, with a cached regridder. |
| `src/arotake/stats_utils.py` | Regional and model-versus-observation statistics containers, consistency checks, pandas/xarray exports, and NetCDF serialization. |
| `src/arotake/plotting.py` | Headless Matplotlib/Cartopy station maps, time series, and Gaussian KDE probability-density plots. |
| `src/arotake/__init__.py` | Package version from distribution metadata, with `0.0.0+local` fallback. |
| `tests/test_test.py`, `src/arotake/test.py` | One test for an addition helper; no scientific regression suite is supplied. |
| `notebooks/ArotakeDemo.ipynb` | Download-and-analysis example using external HRRR/LCD packages; includes saved outputs. |
| `docs/`, `plots/` | Existing example images. |
| `scripts/`, `Makefile` | Environment setup, an exploratory data script, formatting/linting/testing, and Anaconda publishing. |
| `pyproject.toml`, `environment.yml` | Packaging, dependencies, tool configuration, and Conda environment. |
| `.github/workflows/ci.yml`, `.pre-commit-config.yaml` | CI checks and local hooks. |

## Architecture and Data Contracts

`run_analysis` loads observations once, iterates inclusively over daily HRRR initialization dates at a fixed initialization hour and forecast lead, opens each forecast file, and extends a `ModelVsObs2DStatisticsTimeSeries`. That container selects the observation at the forecast valid time, interpolates the model to station locations, and calculates bias, RMSE, correlation, means, standard deviations, variances, and reporting-station diagnostics. Results are exported to NetCDF; optional plots go into `<results_stem>.plots/`, with a station map for `conus`. The API returns the results path and plot paths.

`Model2DRegionalStatisticsTimeSeries` separately computes model-only regional statistics using GeoPandas geometry and regionmask. Both statistics containers accumulate values and metadata in `TimeList` objects and expose `extend`, `test_internal_consistency`, `DataFrame`, `DataSet`, and `write_ds2netcdf`. Preserve these existing public names.

Preserve these input/output contracts unless the task explicitly changes them:

- HRRR input layout: `hrrr.<YYYYMMDD>/<TAG>/hrrr.t<II>z.wrfsfcf<FF>.nc`; latitude/longitude names in the analysis module are `gridlat_0`/`gridlon_0`.
- Model variable metadata includes `units`, `initial_time` formatted as `%m/%d/%Y (%H:%M)`, `forecast_time`, and `forecast_time_units == 'hours'`. Statistics derive forecast timing from this metadata, so it must agree with the requested files.
- Observations use a `time` coordinate, station variables `LAT`/`LON`, a global `region` attribute, and variable `units`. The statistics implementation indexes observations as `(time, station)`. Maps also use global `name` and `long_name` attributes.
- Each forecast valid time must have exactly one observation timestamp. The observation range must cover initialization time plus lead time, including dates beyond the initialization range.
- Interpolation expects two final horizontal model dimensions with matching 2-D latitude/longitude arrays; target coordinates are matching 1-D arrays without NaNs. Model and target longitudes must lie in `[-180, 180]`. The interpolation helper supports leading dimensions, but the statistics consumer assumes station-indexable results; verify broader shape support before claiming it.
- Unit strings must match exactly at comparison time. Missing observations are identified using NaNs; no-reporting timesteps produce NaN comparison metrics. Bias is model minus observation.
- In-memory timestamps are UTC-aware; xarray output timestamps are timezone-naive `datetime64[ns]` representing UTC. Forecast lead values are serialized in hours, despite some comments/docstrings referring to seconds.
- Preserve result variable names, global attributes, and per-variable `units`, `long_name`, and `description`. NetCDF plotting discovers variables via `description == 'time series'`.
- DataFrame plotting prefers `df.attrs['column_attrs']` and falls back to legacy `Series.attrs`. Preserve both compatibility paths and global metadata when changing exports.

## Environment, Setup, and External Interfaces

[pyproject.toml](pyproject.toml) requires Python **3.12 or newer**, setuptools >=77, and wheel. CI exercises Python 3.12. Ruff's `py311` target is a separate tool setting, not the supported Python minimum.

Runtime dependencies include NumPy, pandas, xarray, SciPy, NetCDF backends, Matplotlib, Cartopy, GeoPandas, pyproj, regionmask >=0.13, and xESMF >=0.9. [environment.yml](environment.yml) additionally lists native ESMF/ESMPy, PROJ, libnetcdf, the external `hrrr-data` and `lcd-v2-data` packages, development tools, and Jupyter support. There is no dependency lockfile in the supplied tree.

Repository-defined setup workflows, subject to the isolation rules:

- `make setup-dev-env` invokes [scripts/setup-dev-env.sh](scripts/setup-dev-env.sh). It requires Conda, prefers Mamba when available, writes `environment-dev.yml`, removes/recreates `arotake-dev`, installs `.[dev]` editable with `--no-deps`, replaces the user Jupyter kernel, and installs pre-commit hooks. **It has destructive environment-replacement and outside-root effects; do not run it as a routine validation step under the isolation boundary.** Its final activation instruction is `conda activate arotake-dev`.
- CI installs with `pip install -e ".[dev]"`. Native geospatial dependency availability through that pip-only path has not been verified by this guide; do not infer successful installation from the workflow definition.
- Packaging uses `setuptools.build_meta` and the `src/` layout. There is no dedicated Makefile build target. `make upload-anaconda` builds and publishes through [scripts/upload2anaconda.sh](scripts/upload2anaconda.sh); it is a release operation, not a local build check.

The core analysis consumes local NetCDF files; downloading/preprocessing is delegated to the notebook's external packages. Their source and service behavior are outside this repository. The notebook uses `ESMFMKFILE` and removes it when it points to a missing file; investigate stale ESMF configuration only within the permitted boundary. No application `.env` loader or central configuration file is present; analysis settings are function/CLI arguments and module constants.

The release script uses `ANACONDA_API_TOKEN` and optional overrides such as `ANACONDA_USER_NAME`, `LICENSE_ID`, `PYPROJECT`, `IMPORT_NAME_OVERRIDE`, `BUILD_ENV_PREFIX`, and `BUILD_ENV_NAME`. It creates external temporary directories and Conda environments and uploads to anaconda.org. Report these dependencies instead of executing the release workflow under project isolation; never record credentials in source or logs.

## Coding Conventions

- Follow [.editorconfig](.editorconfig): UTF-8, LF, final newline, four-space indentation, and two spaces for YAML/JSON/TOML. [.gitattributes](.gitattributes) also specifies LF.
- Follow Ruff settings in `pyproject.toml`: line length 128, preserved quote style, rules `E`, `F`, `I`, `UP`, `B`, and `SIM`, with `E501` and `SIM108` ignored. Respect its excluded data/documentation/notebook/output directories.
- Match nearby patterns: snake_case functions/variables, typed public interfaces where already used, `pathlib.Path` filesystem arguments, NumPy-style scientific docstrings, and xarray/pandas metadata propagation.
- Statistics and analysis currently use assertions for consistency checks. Preserve checks and assess error behavior when modifying them; optimized Python can disable assertions.
- Regridding and regional masks are cached using function/class attributes. Changes must account for grid, coordinate, and region changes across successive calls. Do not assume these shared caches are safe for concurrent use.
- Keep plotting headless with the existing Agg backend and close created figures. Preserve the BSD-3-Clause licensing in [LICENSE](LICENSE).

## Validation Commands

Run from the project root in an already available, permitted environment. Report missing dependencies instead of installing them automatically. Keep any generated files within the project boundary.

| Command | Verified behavior |
| --- | --- |
| `python -m ruff check .` | CI lint check without autofixes. |
| `pytest` or `make test` | Runs tests under `tests/`; no coverage options are configured. |
| `make fmt` | Runs `python -m ruff format`; modifies files. |
| `make lint` | Runs `python -m ruff check --fix`; modifies files. |
| `make check` | Runs both formatting and autofixing lint targets; modifies files. |

The pre-commit configuration adds YAML/JSON checks, newline/whitespace fixes, and Ruff hooks. Running it may fetch external hook repositories. README claims about `make type`, strict mypy, and automatic coverage are inconsistent with the Makefile/configuration: there is no `type` target, mypy is commented out, and `make test` is plain pytest. CI runs lint and pytest; it does not establish scientific correctness or a successful end-to-end analysis.

## Change and Review Workflow

1. Read the relevant module and its callers, plus README/API examples. Keep changes focused on the requested behavior and avoid incidental notebook-output or image changes.
2. Preserve CLI positional arguments, public signatures, statistics definitions, UTC/lead-time handling, metadata schemas, and output naming unless a change is requested. Update local examples/documentation when behavior changes.
3. For scientific changes, add focused tests with small synthetic grids/station datasets covering the affected behavior. Useful cases include missing or duplicate timestamps, all-missing observations, incompatible units, changed interpolation grids/targets, metadata retention, and NetCDF round trips. Avoid requiring downloads in unit tests.
4. Run applicable lint/tests when the environment and boundary permit. Review any formatter/autofix diff. Use a small local integration dataset only if supplied inside the root; keep its outputs isolated within the root because writers overwrite existing result files.
5. Review for accidental generated artifacts, unrelated edits, credential exposure, and compatibility changes. Report what changed, checks performed, failures, and checks blocked by missing data or external dependencies.

## Known Pitfalls and Manual Verification

These are source observations or explicit verification gaps, not validated runtime outcomes:

- Temperature conversion in `run_analysis` converts observations from `C`/`degC` to `K` before asserting equality with HRRR units; model conversion happens after that assertion. Celsius model inputs can therefore fail despite README's general conversion claim. Test ordering and units explicitly when changing this path.
- `arg_parse(argv=None)` calls `parser.parse_args()` without passing `argv`, and `main(argv=None)` uses `sys.argv`. Do not assume supplied argument lists are honored.
- `generate_plots=False` returns `[None]` for plot paths although the return annotation says `list[Path]`. Confirm the intended compatibility contract before changing it.
- NetCDF writers overwrite existing files and remove conflicting time attributes from the passed dataset before serialization. Serialization changes require a round-trip check of timestamps, lead-time units, metadata, and numeric precision.
- PDF plotting passes finite samples directly to `gaussian_kde`; empty, insufficient, or constant samples have no explicit guard. Verify edge cases when modifying PDF generation.
- No scientific fixtures or input NetCDF datasets are supplied, and the sole test covers only `add(2, 2)`. Scientific accuracy, native dependency installation, CLI execution, and notebook reproducibility require manual validation with appropriate in-root inputs.
- The notebook declares `lcd-v2-data` but imports `lcd_data`; compatibility depends on that external distribution and cannot be established from this repository. Its download cells use network access and write data; do not execute them merely to validate an unrelated change.
- `scripts/test.py` is an exploratory script, not the pytest suite. It imports `isd_lite_data` (absent from declared dependencies) and refers to `../data/ISD-LITE/data` when run from the root. Do not follow that outside-root data dependency; report it.
