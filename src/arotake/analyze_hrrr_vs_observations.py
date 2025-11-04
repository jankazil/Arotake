#!/usr/bin/env python

'''
HRRR Forecast vs. Observations: Time-Series Analysis and Plotting
=================================================================

This module provides a command-line interface (CLI) and an application
programming interface (API) for comparing HRRR (High-Resolution Rapid
Refresh) model forecasts against station-based observation time series
over a requested period and region.

Functionality
-------------

Given
  - a range of HRRR *initialization* dates,
  - a single initialization hour (UTC),
  - a forecast lead time (hours),
  - an HRRR region tag (for example, "conus"),
  - a directory that contains HRRR forecast NetCDF files,
  - the target HRRR variable name,
  - a NetCDF file with hourly observations for stations within the HRRR domain,
  - the corresponding observation variable name, and
  - an output directory,

the tool:

  1. Loads the observation dataset into an `observations.TimeSeries` object.
  2. Iterates day-by-day through HRRR initialization dates, opening the matching
     HRRR file for the specified init hour and lead time.
  3. Verifies unit consistency between model and observations. Temperature
     variables provided in degrees Celsius may be converted to Kelvin as
     applicable.
  4. Computes model-versus-observation spatial statistics for each forecast
     valid time using `stats_utils.ModelVsObs2DStatisticsTimeSeries` and
     aggregates those results across the period.
  5. Writes the aggregated statistics to a NetCDF file in the output directory.
  6. Optionally generates time-series plots of the computed statistics, and, if
     the region is "conus", a map of station locations.

HRRR file layout
----------------

HRRR forecast files are expected in the given HRRR forecast directory at:

    hrrr.<YYYYMMDD>/<TAG>/hrrr.t<II>z.wrfsfcf<FF>.nc

where:
    YYYYMMDD = initialization date,
    TAG      = HRRR region tag (for example, "conus"),
    II       = initialization hour (UTC, zero-padded),
    FF       = forecast lead time in hours (zero-padded).

Outputs
-------

  - A NetCDF file with aggregated HRRR-vs-observation statistics.
  - If enabled, one or more PNG plots written under a sibling `.plots` directory.

Usage
-----

Command line:

    ./analyze_hrrr_vs_observations.py START_YEAR START_MONTH START_DAY \
        END_YEAR END_MONTH END_DAY FORECAST_INIT_HOUR FORECAST_LEAD_HOUR \
        HRRR_REGION HRRR_DATA_DIR HRRR_VAR OBS_FILE OBS_VAR OUT_DIR

Programmatic:

    from arotake.analyze_hrrr_vs_observations import run_analysis
    results_file, plot_files = run_analysis(
        forecast_init_start_date, forecast_init_end_date, forecast_init_hour, forecast_lead_hour,
        hrrr_region, hrrr_data_dir, hrrr_var, obs_file, obs_var, out_dir,
        verbose=True, generate_plots=True
    )

Notes
-----

- Observation and HRRR variables must be in compatible physical units.
- The observation record must fully cover the requested HRRR *forecast valid*
  times implied by the initialization times plus the specified lead time.

Entry points
------------

- `run_analysis(...)` performs the model-vs-observations comparison and produces outputs.
- `main(argv=None)` implements the command-line interface and calls `run_analysis`.

'''

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


def run_analysis(
    forecast_init_start_date: datetime,
    forecast_init_end_date: datetime,
    forecast_init_hour: int,
    forecast_lead_hour: int,
    hrrr_region: str,
    hrrr_data_dir: Path,
    hrrr_var: str,
    obs_file: Path,
    obs_var: str,
    out_dir: Path,
    verbose: bool = True,
    generate_plots: bool = True,
) -> tuple[Path, list[Path]]:
    '''
    Compute and, if requested, plot HRRR-versus-observation time-series statistics for a date range.

    This routine loads hourly station observations, iterates over daily HRRR
    initialization dates at a fixed initialization hour and forecast lead time,
    computes spatial statistics comparing HRRR fields to colocated station
    observations for each forecast valid time, aggregates the results across
    the requested period, writes the aggregated statistics to NetCDF, and
    optionally generates diagnostic plots. If the region is "conus", an
    additional plot showing station locations is produced.

    Parameters
    ----------
    forecast_init_start_date : datetime.datetime (timezone-aware, UTC)
        First HRRR *initialization* date to include. Time portion is irrelevant.
    forecast_init_end_date : datetime.datetime (timezone-aware, UTC)
        Last HRRR *initialization* date to include (inclusive). Time portion is irrelevant.
    forecast_init_hour : int
        HRRR initialization hour (UTC) to process, 0–23.
    forecast_lead_hour : int
        Forecast lead time in hours.
    hrrr_region : str
        HRRR region tag (for example, "conus").
    hrrr_data_dir : pathlib.Path
        Root directory containing HRRR NetCDF files organized as
        `hrrr.<YYYYMMDD>/<TAG>/hrrr.t<II>z.wrfsfcf<FF>.nc`.
    hrrr_var : str
        Name of the target variable in the HRRR NetCDF files.
    obs_file : pathlib.Path
        NetCDF file with full-hourly UTC station observations covering all
        forecast valid times.
    obs_var : str
        Name of the observed variable in `obs_file`.
    out_dir : pathlib.Path
        Output directory to receive the aggregated NetCDF file and the `.plots`
        directory with generated figures. The directory is created if missing.
    verbose : bool, optional
        If True, print progress information to stdout. Defaults to True.
    generate_plots : bool, optional
        If True, produce time-series plots (and a CONUS station-location plot
        when applicable). Defaults to True.

    Behavior and assumptions
    ------------------------
    - Units must be consistent between HRRR and observations. Temperature
      variables provided in degrees Celsius may be converted to Kelvin as
      applicable.
    - Observation metadata must include station latitude and longitude
      variables compatible with those expected by the statistics routine.
    - HRRR latitude and longitude variable names are expected to be
      `gridlat_0` and `gridlon_0` in the forecast files used here.
    - If `verbose` is True, prints progress information to stdout.

    Results
    -------
    - Writes an aggregated statistics NetCDF file to `out_dir`.
    - If `generate_plots` is True, writes one or more PNG plots to
      `out_dir/<results_basename>.plots/`.

    Returns
    -------
    tuple[pathlib.Path, list[pathlib.Path]]
        A tuple `(results_file_path, plot_file_paths)` where:

          - `results_file_path` is the path to the aggregated statistics NetCDF.
          - `plot_file_paths` is a list of paths to the generated plot files, or
            the list `[None]` when `generate_plots` is False.

    Raises
    ------
    AssertionError
        If HRRR and observation units are incompatible after any automatic
        Celsius-to-Kelvin conversion.
    FileNotFoundError
        If an expected HRRR file is missing.
    xarray.backends.plugins.BackendEntrypointError or OSError
        If an HRRR file cannot be opened as a NetCDF dataset.
    KeyError
        If `hrrr_var` or `obs_var` is not present in the provided datasets.

    Examples
    --------
    >>> from pathlib import Path
    >>> from datetime import datetime, timezone
    >>> run_analysis(
    ...     forecast_init_start_date=datetime(2021, 1, 1, tzinfo=timezone.utc),
    ...     forecast_init_end_date=datetime(2021, 1, 3, tzinfo=timezone.utc),
    ...     forecast_init_hour=0,
    ...     forecast_lead_hour=1,
    ...     hrrr_region="conus",
    ...     hrrr_data_dir=Path("/data/hrrr"),
    ...     hrrr_var="TMP_P0_L103_GLC0",
    ...     obs_file=Path("/data/obs/isd_conus.nc"),
    ...     obs_var="T",
    ...     out_dir=Path("/data/output"),
    ...     verbose=True,
    ...     generate_plots=True,
    ... )
    '''

    # Set the time portion of the first and last HRRR initialization date to midnight
    forecast_init_start_date_midnight = forecast_init_start_date.replace(hour=0, minute=0, second=0, microsecond=0)
    forecast_init_end_date_midnight = forecast_init_end_date.replace(hour=0, minute=0, second=0, microsecond=0)

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

    forecast_creation_date = forecast_init_start_date_midnight

    # Loop increment

    time_step = timedelta(days=1)

    while forecast_creation_date <= forecast_init_end_date_midnight:
        # Forecast initilization time
        forecast_init_time = forecast_creation_date + forecast_init_time_delta

        # Forecast valid time
        forecast_valid_time = forecast_creation_date + forecast_init_time_delta + forecast_lead_time_delta

        if verbose:
            print(
                'Evaluating for the country/state/territory/region: '
                + data_obs.time_series.attrs['region']
                + ', forecast init time '
                + str(forecast_init_time)
                + ', forecast valid time '
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

    if verbose:
        print()
        print('Saved results in', str(results_file_path))

    #
    # Plotting
    #

    plot_file_paths = [None]

    if generate_plots:
        plot_dir = results_file_path.parent / Path(str(results_file_path.stem) + '.plots')

        plot_file_paths = plotting.plot_df_timeseries(
            plot_dir, [hrrr_vs_obs_statistics.DataFrame()], [data_obs.time_series.attrs['region']]
        )

        if hrrr_region == 'conus':
            plot_file_path = plotting.plot_locations_conus(
                data_obs,
                forecast_init_start_date_midnight + forecast_init_time_delta + forecast_lead_time_delta,
                forecast_init_end_date_midnight + forecast_init_time_delta + forecast_lead_time_delta,
                obs_lat_name,
                obs_lon_name,
                plot_dir,
            )
            plot_file_paths.append(plot_file_path)

        if verbose:
            print()
            print('Created the following plots:')
            print()

            for plot in plot_file_paths:
                print(plot)

    return results_file_path, plot_file_paths


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

    forecast_init_start_date = datetime(year=args.start_year, month=args.start_month, day=args.start_day, tzinfo=timezone.utc)
    forecast_init_end_date = datetime(year=args.end_year, month=args.end_month, day=args.end_day, tzinfo=timezone.utc)
    forecast_init_hour = args.forecast_init_hour
    forecast_lead_hour = args.forecast_lead_hour
    hrrr_region = args.hrrr_region
    hrrr_data_dir = Path(args.hrrr_data_dir)
    hrrr_var = args.hrrr_var
    obs_file = Path(args.obs_file)
    obs_var = args.obs_var
    out_dir = Path(args.out_dir)

    return (
        forecast_init_start_date,
        forecast_init_end_date,
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
        forecast_init_start_date,
        forecast_init_end_date,
        forecast_init_hour,
        forecast_lead_hour,
        hrrr_region,
        hrrr_data_dir,
        hrrr_var,
        obs_file,
        obs_var,
        out_dir,
    ) = arg_parse(sys.argv[1:])

    results_file, plot_files = run_analysis(
        forecast_init_start_date,
        forecast_init_end_date,
        forecast_init_hour,
        forecast_lead_hour,
        hrrr_region,
        hrrr_data_dir,
        hrrr_var,
        obs_file,
        obs_var,
        out_dir,
    )


if __name__ == '__main__':
    main()
