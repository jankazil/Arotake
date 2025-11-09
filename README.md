# Arotake

**Arotake** is a Python toolkit for evaluating meteorological forecast models against observations. It provides a command-line tool and an application programming interface that evaluates NOAA High-Resolution Rapid Refresh (HRRR) surface forecasts vs. observations.

## Installation

```bash
mamba install -c jan.kazil -c conda-forge Arotake
```

## Overview

The toolkit provides a top level command line (CLI) tool and application programming interface (API) for evaluating NOAA High-Resolution Rapid Refresh (HRRR) contiguous United States surface forecasts against surface observations:

- **analyze-hrrr-vs-observations**  
  Loads HRRR forecasts and observations, computes HRRR-vs-observation statistics time series, plots the results, and saves the results as a NetCDF file.

## Workflow (using API)

The workflow is shown in the [Arotake Demo](https://github.com/jankazil/Arotake/blob/main/notebooks/ArotakeDemo.ipynb) Jupyter notebook. The notebook evaluates High Resolution Rapid Refresh (HRRR) forecast data against Local Climatological Data (LCD) v2 observations. The workflow is

- Specify settings and working directories  
- Specify HRRR forecast parameters (analysis period, forecast initialization time, forecast valid time, ...)  
- Download and process HRRR forecast data
- Specify LCD v2 observation parameters (analysis period, analysis regions)  
- Download and process LCD v2 data
- Evaluate the HRRR forecast against the LCD v2 observations. This creates and saves time series as netCDF files and creates plots, for each analysis region, of the model bias, rmse, correlation with observations, etc.  
- Create a plot with the map of the LCD v2 observation locations.  
- Create time series and probability density functions plots of the model bias, rmse, correlation with observations, etc., for all analysis regions concurrently.  

Requirements: Jupyter kernel with the following packages installed:  

- [`Arotake`](https://github.com/jankazil/Arotake)  
- [`hrrr-data`](https://github.com/jankazil/hrrr-data)  
- [`lcd-data`](https://github.com/jankazil/lcd-data)

## Command-line interface (CLI)

### `analyze-hrrr-vs-observations`

**Description**

Compares HRRR forecasts with hourly surface observations, computes time series of model-observation statistics, and generates plots summarizing the results.  

This tool supports both ISD-Lite and LCDv2 observation datasets, provided as NetCDF files (generated, e.g., using [`isd-lite-data`](https://github.com/jankazil/isd-lite-data) or [`lcd-data`](https://github.com/jankazil/lcd-data)).

**Usage**

```bash
analyze-hrrr-vs-observations   <start_year> <start_month> <start_day>   <end_year> <end_month> <end_day>   <forecast_init_hour> <forecast_lead_hour>   <hrrr_region> <hrrr_data_dir> <hrrr_var>   <obs_file> <obs_var> <out_dir>
```

| Argument | Description |
|:----------|:------------|
| `start_year`, `start_month`, `start_day` | Start date (UTC) of the HRRR initialization range. |
| `end_year`, `end_month`, `end_day` | End date (UTC) of the HRRR initialization range. |
| `forecast_init_hour` | HRRR forecast initialization hour (UTC). |
| `forecast_lead_hour` | HRRR forecast lead time in hours. |
| `hrrr_region` | HRRR region tag (e.g. `conus`). |
| `hrrr_data_dir` | Directory containing HRRR forecast NetCDF files. |
| `hrrr_var` | Name of the variable in HRRR files (e.g. `TMP_P0_L103_GLC0` for 2-m temperature). |
| `obs_file` | Path to the NetCDF file with full-hourly UTC observations. |
| `obs_var` | Variable name in the observation file (e.g. `T`). |
| `out_dir` | Output directory for NetCDF and plots (created if absent). |

**HRRR directory structure**

HRRR forecast files must be available in subdirectories of `hrrr_data_dir` with the following path/name convention:

```
hrrr.<YYYYMMDD>/<TAG>/hrrr.t<II>z.wrfsfcf<FF>.nc
```

where:
- `YYYYMMDD` = initialization date  
- `TAG` = HRRR region tag (e.g. `conus`)  
- `II` = initialization hour (UTC)  
- `FF` = forecast lead time (hours)

**Example**

```bash
analyze-hrrr-vs-observations 2020 1 1 2020 12 31 0 12 conus data/ TMP_P0_L103_GLC0 ERCOT.2024-2024.nc T results/

```

**Output**

- A NetCDF file containing time series of HRRR-versus-observation statistics at observation locations.  
- Time series plots of bias, RMSE, correlation, and related metrics saved in a `.plots` subdirectory.  
- Map of observation station locations.

**Notes**

- HRRR and observation variables must share consistent units. Temperatures in Celsius are automatically converted to Kelvin.  
- The observation dataset must fully cover the HRRR forecast validation period.  
- This command can be run either as the installed entry point `analyze-hrrr-vs-observations` or directly via the Python script.

## Public application programming interface (API)

### Module `arotake.analyze_hrrr_vs_observations`  

####  Functions

- `run_analysis`: API function that performs the same HRRR-versus-observation evaluation as the CLI `analyze-hrrr-vs-observations`, but callable directly from Python code.

    ***Purpose***

    Computes forecast-observation statistics (bias, RMSE, correlation, etc.) for HRRR surface fields relative to hourly observations     and generates corresponding plots and NetCDF outputs.

    ***Returns***

    - NetCDF file containing HRRR-versus-observation statistics time series.  
    - Diagnostic plots (bias, RMSE, correlation) saved in the output directory.  
    - Optional in-memory DataFrame for further analysis.

### Module `arotake.interpolation`  

####  Functions

  - `model_2D_interpolate`: Interpolates a 2-D model field to arbitrary lat/lon points using bilinear regridding with xESMF, with internal caching to improve performance.

### Module `arotake.plotting`  

####  Functions

  - `plot_df_timeseries`: Plots one figure per statistic by overlaying multiple DataFrames and saves PNG files with consistent naming conventions.
  - `plot_locations_conus`: Creates a map of observation locations over the contiguous U.S. in Lambert Conformal projection.

## Development

### Code Quality and Testing Commands

- `make fmt` - Runs `ruff format` to auto-format Python files.
- `make lint` - Runs `ruff check --fix` to lint and autofix style issues.
- `make check` - Runs both formatting and linting.
- `make type` - Runs `mypy` type checker in strict mode.
- `make test` - Runs `pytest` with coverage reporting.

## Name

Arotake is a word in te reo Māori meaning “review” or “evaluate.” The name reflects the project’s focus on the evaluation of meteorological models.

## Author

Jan Kazil - jan.kazil.dev@gmail.com - [jankazil.com](https://jankazil.com)

## License

BSD 3-Clause
