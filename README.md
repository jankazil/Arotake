# Arotake

**Arotake** is a Python toolkit for evaluating meteorological forecast models against observations. It provides

- A command-line tool that evaluates NOAA High-Resolution Rapid Refresh (HRRR) surface forecasts vs. observations

## Installation (Linux / macOS)

```bash
mamba install -c jan.kazil -c conda-forge Arotake
```

## Overview

The toolkit provides a top level command line tool for evaluating NOAA High-Resolution Rapid Refresh (HRRR) contiguous United States surface forecasts against surface observations:

- **analyze-hrrr-vs-observations**  
  Loads HRRR forecasts and observations, computes HRRR-vs-observation statistics time series, plots the results, and saves the results as a NetCDF file.

## Workflow

The workflow is

1. Download HRRR contiguous United States surface forecast data with [hrrr-data](https://github.com/jankazil/hrrr-data)
2. Download and construct full-hourly UTC time series in a given U.S. state, territory, RTO/ISO region, or individual observation locations using one of the following tools:
   - [`isd-lite-data`](https://github.com/jankazil/isd-lite-data)
   - [`lcd-data`](https://github.com/jankazil/lcd-data)
3. Calculate the HRRR-vs-observation statistics time series using this package.

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

    ***Notes***

    - The function mirrors the CLI’s workflow and argument structure.  

    ***Example***

```python
from arotake.analyze_hrrr_vs_observations import run_analysis

start_date = datetime(year=2021, month=1, day=1, tzinfo=timezone.utc)
end_date = datetime(year=2021, month=12, day=31, tzinfo=timezone.utc)
forecast_init_hour = 12
forecast_lead_hour = 32
hrrr_region = 'conus'
hrrr_data_dir = Path('HRRR') / Path('data')
hrrr_var = 'TMP_P0_L103_GLC0'
obs_file = Path('LCD') / Path('data') / Path('CO.2020-2025.nc')
obs_var = 'T'
out_dir = Path('results')

results_file, plot_files = run_analysis(
    start_date,
    end_date,
    forecast_init_hour,
    forecast_lead_hour,
    hrrr_region,
    hrrr_data_dir,
    hrrr_var,
    obs_file,
    obs_var,
    out_dir,
    verbose=True,
    generate_plots=True
)
```

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
