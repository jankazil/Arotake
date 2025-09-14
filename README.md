# Arotake

**Arotake** is a Python toolkit for evaluating meteorological forecast models against observations. It currently offers a comparison of select NOAA HRRR surface forecasts with NCEI ISD-Lite surface observations, stratified by U.S. Regional Transmission Organization (RTO) and Independent System Operator (ISO) regions.

## Overview

The repository contains:

- **Analysis scripts** (`scripts/`):
  - `Construct_ISDLite_data_netCDF_for_RTO_ISO_regions.py`  

    Builds region-specific ISD-Lite datasets 2020-12-02 -- 2025-06-30:
    - Loads ISD station metadata, filters by geographic bounds and data availability.
    - Assigns stations to regions, saves per-region metadata, generates a map, and writes consolidated NetCDF observation files.

  - `Analyze_HRRR_vs_ISDLite_time_series_by_RTO_ISO_region.py`  

    Performs HRRR vs ISD-Lite comparison 2020-12-03 -- 2021-12-31:
    - Loads HRRR forecast files and ISD-Lite station NetCDF data for each region.
    - Computes regional HRRR-only statistics and model-vs-obs metrics.
    - Aggregates results into pandas DataFrames, saves NetCDF time series, and generates time-series plots.

## Usage

Scripts are executed as follows

1. **Construct ISD-Lite region datasets**  

   ```cd scripts  ```  
   ```./Construct_ISDLite_data_netCDF_for_RTO_ISO_regions.py```

   Produces:
   - A map of ISD station locations in each RTO/ISO region (`RTO_ISO_regions_ISD_stations_map.png`)
   - ISD metadata text files for each RTO/ISO region (`data/RTO_ISO_regions_ISD_stations/*.txt`)
   - netCDF with ISD-Lite observations for each region (`data/RTO_ISO_regions_ISD_stations/*.nc`)

   This script needs not to be run again unless additional ISD-Lite datasets need to be produced.

2. **Run HRRR vs ISD-Lite analysis**  

   ```cd scripts  ```  
   ```./Analyze_HRRR_vs_ISDLite_time_series_by_RTO_ISO_region.py```

   Produces:
   - Plots of HRRR vs ISD-Lite statistics time series for each RTO/ISO region (`./plots/<regions>/`)
   - NetCDF files with the HRRR vs ISD-Lite statistics time series for each RTO/ISO region (`./data/RTO_ISO_regions_timeseries/*.nc`)

## Data Preparation

The above scripts operate on HRRR forecast data and ISD-Lite observations stored in netCDF files. The data can be downloaded and converted to netCDF files with the following toolkits:

- [HRRR-data](https://github.com/jankazil/hrrr-data)
  - Use the script [DownloadHRRRv4Data.py](https://github.com/jankazil/hrrr-data/blob/main/scripts/DownloadHRRRv4Data.py) to download HRRR forecast files in GRIB format for a given time range from the Amazon S3 HRRR bucket, and to extract select variables into netCDF files. Name/place the files

    ./data/HRRR/data/hrrr.<YYYYMMDD\>/conus/hrrr.t<II\>z.wrfsfcf<FF\>_select_vars.nc

    The download script creates the directory structure with the correct directory/file names automatically. <YYYYMMDD\> are the year, month, and day, <II\> the forecast initialization time in hours, and <FF\> the forecast time in hours.
- [ISD-Lite-data](https://github.com/jankazil/isd-lite-data)  
Use the regional download scripts for the Contiguous United States (CONUS) to download the ISD-Lite observations and create a netCDF file covering all CONUS stations:

  - [CONUS_Metadata_Download.py](https://github.com/jankazil/isd-lite-data/blob/main/scripts/CONUS_Metadata_Download.py)  
    Downloads the ISDLite station metadata and produces the file conus_stations.2020-2025.txt  
  - [CONUS_Observations_Download.py](https://github.com/jankazil/isd-lite-data/blob/main/scripts/CONUS_Observations_Download.py)  
    Downloads the ISDLite station observation files <USAF_ID\>-<WBAN_ID\>-<YYYY\>.gz where <USAF_ID\> is the USAF station ID, <WBAN_ID\> the WBAN ID, and <YYYY\> the year. Place the files in

    ./data/ISD-LITE/data/

- Download the definitions of the RTO/ISO regions in GEOJson format from https://atlas.eia.gov/datasets/rto-regions, and name/place the file in  

    ./data/RTO_ISO_regions.geojson

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

## Name

Arotake is a word in te reo Māori meaning “review” or “evaluate.” The name reflects the project’s focus on the evaluation of meteorological models.

## Author

Jan Kazil – jan.kazil.dev@gmail.com – [jankazil.com](https://jankazil.com)

## License

BSD 3-Clause
