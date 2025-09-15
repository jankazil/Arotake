# Arotake

**Arotake** is a Python toolkit for evaluating meteorological forecast models against observations. It currently offers a comparison of select NOAA HRRR surface forecast variables with NCEI ISD-Lite surface observations, stratified by U.S. Regional Transmission Organization (RTO) and Independent System Operator (ISO) regions.

## Installation (Linux / macOS)

```bash
wget -O environment.yml https://raw.githubusercontent.com/jankazil/Arotake/main/environment.yml  
mamba env create -y -f environment.yml  
conda activate arotake
```

## Overview

The toolkit provides analysis scripts for working with NOAA ISD-Lite surface observations and NOAA HRRR surface forecasts, organized by U.S. Regional Transmission Organization (RTO) and Independent System Operator (ISO) regions.

- **Construct_ISDLite_data_netCDF_for_RTO_ISO_regions.py**  
  Downloads ISD-Lite station metadata, filters stations by location and data availability, assigns stations to RTO/ISO regions, generates a regional station map, downloads the observations (optionally in parallel), and saves region-specific NetCDF files with observations. This needs to be done only once for the given time range.

- **Analyze_HRRR_vs_ISDLite_time_series_by_RTO_ISO_region.py**  
  Loads region-specific ISD-Lite observation NetCDF files and HRRR forecast files, computes regional HRRR-only statistics and HRRR-vs-observation diagnostics, saves results as NetCDF time series, and generates plots of the statistics for each region.

## Usage

The scripts accept on the command line arguments specifying the time range, input file locations, and output directories.

### 1. Construct RTO/ISO region ISD-Lite Datasets

```bash
Construct_ISDLite_data_netCDF_for_RTO_ISO_regions.py <start_year> <start_month> <start_day> <end_year> <end_month> <end_day> <geojson_file> <isdlite_data_dir> [-n <n_jobs>]
```

Example:

```bash
Construct_ISDLite_data_netCDF_for_RTO_ISO_regions.py 2021 1 1 2021 12 31 data/RTO_ISO_regions.geojson data/ISD-LITE/ -n 8
```

This will:  

- Load RTO/ISO region geometries.
- Download and filter ISD-Lite station metadata.  
- Assign stations to RTO/ISO regions and save station lists.  
- Plot a map of the ISD stations in each region (`ISD_stations_map.png`).  
- Download observations (parallelized if `-n` > 1).  
- Save region-specific NetCDF files with the observations.  

### 2. Calculate RTO/ISO region HRRR vs ISD-Lite Statistics Time Series

```bash
Analyze_HRRR_vs_ISDLite_time_series_by_RTO_ISO_region.py <start_year> <start_month> <start_day> <end_year> <end_month> <end_day> <forecast_init_hour> <forecast_lead_hour> <geojson_file> <isdlite_data_dir> <hrrr_data_dir> <out_dir>
```

<hrrr_data_dir\> is the parent directory of the HRRR data directory, which has the following structure (see Section ''Data Preparation''):  

hrrr.<YYYYMMDD\>/conus/hrrr.t<II\>z.wrfsfcf<FF\>_select_vars.nc

Example:

```bash
Analyze_HRRR_vs_ISDLite_time_series_by_RTO_ISO_region.py 2021 1 1 2021 12 30 12 6 data/RTO_ISO_regions.geojson data/ISD-LITE/ data/HRRR/ results/
```

This will:  

- Load RTO/ISO region geometries and ISD-Lite observations.
- Load HRRR forecasts for the specified initialization/lead hours.
- Compute and save time series of HRRR regional statistics and of HRRR vs ISD-Lite observation statistics.
- Generate plots of the time series per region in `./results/plots/<region>/`.

## Data Preparation

- The definitions of the RTO/ISO regions in GEOJson format can be downloaded from https://atlas.eia.gov/datasets/rto-regions.

- The HRRR forecast data for the contiguous United States, the time period, forecast initialization time, and forecast lead time can be downloaded and converted to netCDF files with the [HRRR-data](https://github.com/jankazil/hrrr-data) toolkit.
    - Use the script    [DownloadHRRRSurfaceForecast.py](https://github.com/jankazil/hrrr-data/blob/main/scripts/DownloadHRRRSurfaceForecast.py) to download HRRR surface forecast files in GRIB format for a given time range from the Amazon S3 HRRR bucket, which also extracts select variables into netCDF files:

        hrrr.<YYYYMMDD\>/conus/hrrr.t<II\>z.wrfsfcf<FF\>_select_vars.nc

## Public API

### Modules

- **`arotake.interpolation`**  
  Function:
  - `model_2D_interpolate`: Interpolates a 2-D model field to arbitrary lat/lon points using bilinear regridding with xESMF, with internal caching to improve performance.

- **`arotake.rto_iso`**  
  Functions:
  - `regions`: Loads all RTO/ISO regions from a GeoJSON file and returns a `GeoDataFrame`.
  - `region`: Loads a single RTO/ISO region by name and returns a `GeoDataFrame`.

- **`arotake.statistics`**  
  Classes:
  - `Model2DRegionalStatisticsTimeSeries`: Computes and accumulates regional statistics (mean, std, var) of a model variable over a polygonal region, exports results as DataFrames or Datasets.
  - `ModelVsObs2DStatisticsTimeSeries`: Computes model-vs-obs metrics (bias, RMSE, correlation, coverage) by interpolating model data to station locations and comparing with observations.

- **`arotake.plotting`**  
  Functions:
  - `plot_df_timeseries`: Plots one figure per statistic by overlaying multiple DataFrames and saves PNG files with consistent naming conventions.

## Development

### Code Quality and Testing Commands

- `make fmt` – Runs `ruff format` to auto-format Python files.
- `make lint` – Runs `ruff check --fix` to lint and autofix style issues.
- `make check` – Runs both formatting and linting.
- `make type` – Runs `mypy` type checker in strict mode.
- `make test` – Runs `pytest` with coverage reporting.

### Notes

Arotake uses [xESMF](https://xesmf.readthedocs.io]) (Universal Regridder for Geospatial Data), which uses [ESMPy](https://earthsystemmodeling.org/esmpy) (ESMF Python Regridding Interface) as backend. These packages are as of writing available as conda packages from conda-forge. This necessitates installing a conda environment to operate Arotake.

## Name

Arotake is a word in te reo Māori meaning “review” or “evaluate.” The name reflects the project’s focus on the evaluation of meteorological models.

## Author

Jan Kazil – jan.kazil.dev@gmail.com – [jankazil.com](https://jankazil.com)

## License

BSD 3-Clause
