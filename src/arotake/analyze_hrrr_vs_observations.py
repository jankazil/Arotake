#!/usr/bin/env python

"""
HRRR Forecast vs. Observations Time Series Analysis and Plotting Tool
=====================================================================

This script compares HRRR (High-Resolution Rapid Refresh) model forecasts
against station-based observation time series and generates summary statistics
and plots for a specified time period and region.

Overview
--------
Given:
  - a range of HRRR initialization dates,
  - an initialization hour,
  - a forecast lead time (in hours),
  - an HRRR region tag (e.g., "conus"),
  - a directory containing HRRR forecast NetCDF files,
  - a variable name in the HRRR files (e.g., "TMP_P0_L103_GLC0" for 2 m temperature),
  - a NetCDF file with hourly observations (e.g., produced by isd-lite-data or lcd-data),
  - a variable name in the observation file (e.g., "T" for 2 m temperature),
  - and an output directory,

the program performs the following steps:

  1. Loads the observation dataset into a `TimeSeries` object.
  2. Iterates over HRRR forecast initialization dates within the specified range.
  3. For each forecast, loads the corresponding HRRR NetCDF file and verifies unit consistency.
  4. Computes model-versus-observation spatial statistics using
     `stats_utils.ModelVsObs2DStatisticsTimeSeries`.
  5. Aggregates results across all forecast cycles.
  6. Saves the resulting statistics as a NetCDF file.
  7. Produces:
       - time series plots of model–observation statistics, and
       - (if the region is "conus") a map showing the locations of observation stations.

File Naming and Directory Structure
-----------------------------------
HRRR forecast files must be organized as:

    hrrr.<YYYYMMDD>/<TAG>/hrrr.t<II>z.wrfsfcf<FF>.nc

where:
    YYYYMMDD = initialization date,
    TAG       = HRRR region tag (e.g., "conus"),
    II        = initialization hour (UTC),
    FF        = forecast lead time (hours).

Observation file:
    A NetCDF file containing hourly UTC observations with variables for latitude,
    longitude, and the observed field of interest.

Output:
    - A NetCDF file with aggregated HRRR vs. observation statistics.
    - One or more PNG plots saved under a `.plots` subdirectory.

Usage
-----
Run from the command line:

    ./analyze_hrrr_vs_observations.py START_YEAR START_MONTH START_DAY \
        END_YEAR END_MONTH END_DAY FORECAST_INIT_HOUR FORECAST_LEAD_HOUR \
        HRRR_REGION HRRR_DATA_DIR HRRR_VAR OBS_FILE OBS_VAR OUT_DIR

Example:
    ./analyze_hrrr_vs_observations.py 2021 1 1 2021 1 3 0 1 conus \
        /data/hrrr TMP_P0_L103_GLC0 /data/obs/isd_conus.nc T /data/output

Notes
-----
- The observation and HRRR datasets must use consistent physical units.
  Temperature variables in degrees Celsius are automatically converted to Kelvin.
- The observation period must fully cover the HRRR forecast validation times.

Entry Points
------------
Functions:
    arg_parse(argv=None)
        Parses and validates command-line arguments, returning structured inputs.

    main(argv=None)
        Main program entry point. Executes the complete HRRR vs. observation comparison
        workflow, saves results, and produces plots.

"""

import argparse
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import xarray as xr

from arotake import observations, plotting, stats_utils

#
# Settings
#

# Keep attributes of xarray Datasets upon numerical operations (otherwise they are lost)

xr.set_options(keep_attrs=True)

# Name of latitude and longitude variables in the HRRR forecastfiles
hrrr_lat_name = 'gridlat_0'
hrrr_lon_name = 'gridlon_0'

# Name of latitude and longitude variables in the observation files
obs_lat_name = 'LAT'
obs_lon_name = 'LON'


def arg_parse(argv=None):
    '''
    Argument parser which returns the parsed values given as arguments.
    '''

    code_description = (
        "Computes, saves as a netCDF file, and plots the time series of the comparison of HRRR forecasts against observations.\n\n"
        "Given:\n\n"
        "  - a range of HRRR initialization dates,\n"
        "  - an HRRR initialization hour,\n"
        "  - the HRRR forecast lead time in hours,\n"
        "  - an HRRR region tag (e.g. 'conus'),\n"
        "  - a path to a directory with HRRR forecast files in netCDF format,\n"
        "  - a name of variable in the HRRR files (e.g. 'TMP_P0_L103_GLC0' for 2 m temperature)\n"
        "  - a path to a file in netCDF format containing full-hourly UTC time series of observations at locations within the HRRR domain,\n"
        "  - a name of variable in the observations file (e.g. 'T' for 2 m temperature),\n"
        "  - a path to an output directory,\n\n"
        "this executable loads the forecast data and observations, computes model-versus-observation statistics time series, writes results "
        "to NetCDF, and produces time series plots.\n\n"
        "The given HRRR forecast data directory must contain files with the path/names\n\n"
        "    hrrr.<YYYYMMDD>/<TAG>/hrrr.t<II>z.wrfsfcf<FF>.nc\n\n"
        "\n\n"
        "where\n\n"
        "    YYYYMMDD is the year, month, and day,\n"
        "    TAG is the HRRR region tag \n"
        "    II is the initialization hour,\n"
        "    FF is the forecast lead time in hours\n\n"
        "These files can be downloaded into the given HRRR forecast data directory using the Python package https://github.com/jankazil/hrrr-data.\n"
        "\n"
        "The time range of HRRR forecast times (not HRRR initialization times) must be fully contained by the time range of the netCDF file holding the observations.\n"
        "\n"
        "The observation files can be generated using the Python packages https://github.com/jankazil/isd-lite-data or https://github.com/jankazil/lcd-data.\n"
    )

    parser = argparse.ArgumentParser(
        description=code_description,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Mandatory arguments
    parser.add_argument('start_year', type=int, help='HRRR initialization start year.')
    parser.add_argument('start_month', type=int, help='HRRR initialization start month.')
    parser.add_argument('start_day', type=int, help='HRRR initialization start day.')
    parser.add_argument('end_year', type=int, help='HRRR initialization end year.')
    parser.add_argument('end_month', type=int, help='HRRR initialization end month.')
    parser.add_argument('end_day', type=int, help='HRRR initialization end day.')
    parser.add_argument('forecast_init_hour', type=int, help='HRRR forecast initialization hour (UTC).')
    parser.add_argument('forecast_lead_hour', type=int, help='HRRR forecast lead time in hours.')
    parser.add_argument('hrrr_region', type=str, help='HRRR region tag (e.g. "conus").')
    parser.add_argument('hrrr_data_dir', type=str, help='Directory where HRRR data files in netCDF format are located.')
    parser.add_argument('hrrr_var', type=str, help='Name of a variable in the HRRR files.')
    parser.add_argument('obs_file', type=str, help='path to netCDF file with time series of observables.')
    parser.add_argument('obs_var', type=str, help='Name of a variable in the observations files.')
    parser.add_argument(
        'out_dir', type=str, help='Directory where results will be saved. Will be created if it does not exist.'
    )

    # Optional arguments
    # parser.add_argument('-x','--xxx', type=str, help='HELP STRING HERE')

    args = parser.parse_args()

    start_date = datetime(year=args.start_year, month=args.start_month, day=args.start_day, tzinfo=timezone.utc)
    end_date = datetime(year=args.end_year, month=args.end_month, day=args.end_day, tzinfo=timezone.utc)
    forecast_init_hour = args.forecast_init_hour
    forecast_lead_hour = args.forecast_lead_hour
    hrrr_region = args.hrrr_region
    hrrr_data_dir = Path(args.hrrr_data_dir)
    hrrr_var = args.hrrr_var
    obs_file = Path(args.obs_file)
    obs_var = args.obs_var
    out_dir = Path(args.out_dir)

    return (
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
    )


def main(argv=None):
    '''
    Command line interface entry point.
    '''

    (
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
    ) = arg_parse(sys.argv[1:])

    #
    # Load observations
    #

    data_obs = observations.TimeSeries(obs_file)

    obs_data_units = data_obs.time_series[obs_var].attrs['units']

    # Convert T to Kelvin
    if obs_data_units == 'C' or obs_data_units == 'degC':
        data_obs.time_series[obs_var] = data_obs.time_series[obs_var] + 273.15
        data_obs.time_series[obs_var].attrs['units'] = 'K'
        obs_data_units = 'K'

    #
    # Loop over HRRR forecasts
    #

    # Initialize empty model evaluation statistics object

    hrrr_vs_obs_statistics = stats_utils.ModelVsObs2DStatisticsTimeSeries(None, None, None, None, None, None, None, None)

    # Timing

    forecast_init_time_delta = timedelta(hours=forecast_init_hour)
    forecast_lead_time_delta = timedelta(hours=forecast_lead_hour)

    # Loop variable

    forecast_creation_date = start_date

    # Loop increment

    time_step = timedelta(days=1)

    while forecast_creation_date <= end_date:
        # Forecast initilization time
        forecast_init_time = forecast_creation_date + forecast_init_time_delta

        # Forecast valid time
        forecast_valid_time = forecast_creation_date + forecast_init_time_delta + forecast_lead_time_delta

        print(
            'Evaluating for the country/state/territory/region: '
            + data_obs.time_series.attrs['region']
            + ', forecast init time '
            + str(forecast_init_time)
            + ' , forecast valid time '
            + str(forecast_valid_time)
        )

        #
        # HRRR data
        #

        # Construct path of HRRR file (in HRRR data directory)

        hrrr_dir_name = (
            'hrrr.'
            + str(forecast_init_time.year).zfill(4)
            + str(forecast_init_time.month).zfill(2)
            + str(forecast_init_time.day).zfill(2)
        )

        hrrr_file_name = 'hrrr.t' + str(forecast_init_hour).zfill(2) + 'z.wrfsfcf' + str(forecast_lead_hour).zfill(2) + '.nc'

        hrrr_file = hrrr_data_dir / hrrr_dir_name / hrrr_region / hrrr_file_name

        # Open HRRR file and load data
        hrrr_ds = xr.open_dataset(hrrr_file)

        # Units

        hrrr_data_units = hrrr_ds[hrrr_var].attrs['units']

        assert hrrr_data_units == obs_data_units, 'Units mismatch between HRRR forecast and observations.'

        # Convert T to Kelvin
        if hrrr_data_units == 'C' or hrrr_data_units == 'degC':
            hrrr_ds[hrrr_var] = hrrr_ds[hrrr_var] + 273.15
            hrrr_ds[hrrr_var].attrs['units'] = 'K'
            hrrr_data_units = 'K'

        #
        # Statistics of HRRR forecast vs observations
        #

        hrrr_vs_obs_statistics.extend(
            stats_utils.ModelVsObs2DStatisticsTimeSeries(
                hrrr_ds,
                hrrr_var,
                hrrr_lat_name,
                hrrr_lon_name,
                data_obs,
                obs_var,
                obs_lat_name,
                obs_lon_name,
            )
        )

        # Close HRRR file

        hrrr_ds.close()

        # Iterate date

        forecast_creation_date += time_step

    # Save results

    results_file_path = hrrr_vs_obs_statistics.write_ds2netcdf(hrrr_vs_obs_statistics.DataSet(), out_dir)

    print()
    print('Saved results in', str(results_file_path))

    #
    # Plotting
    #

    plot_dir = results_file_path.parent / Path(str(results_file_path.stem) + '.plots')

    plot_file_paths = plotting.plot_df_timeseries(
        plot_dir, [hrrr_vs_obs_statistics.DataFrame()], [data_obs.time_series.attrs['region']]
    )

    if hrrr_region == 'conus':
        plot_file_path = plotting.plot_locations_conus(
            data_obs,
            start_date + forecast_init_time_delta + forecast_lead_time_delta,
            end_date + forecast_init_time_delta + forecast_lead_time_delta,
            obs_lat_name,
            obs_lon_name,
            plot_dir,
        )
        plot_file_paths.append(plot_file_path)

    print()
    print('Created the following plots:')
    print()

    for plot in plot_file_paths:
        print(plot)


if __name__ == '__main__':
    main()
