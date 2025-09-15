#!/usr/bin/env python

import argparse
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import xarray as xr
from isd_lite_data import stations

from arotake import plotting, rto_iso, statistics

#
# Settings
#

# Keep attributes of xarray Datasets upon numerical operations (otherwise they are lost)

xr.set_options(keep_attrs=True)

#
# Command line arguments
#


def arg_parse(argv=None):
    '''
    Argument parser which returns the parsed values given as arguments.
    '''

    code_description = (
        'Compute HRRR versus ISD-Lite temperature statistics over U.S. RTO/ISO regions. '
        'Given a start and end date, an HRRR initialization hour and forecast lead time, and paths to the region GeoJSON, '
        'ISD-Lite data, HRRR data, and an output directory, the script loads region geometries, forecast data, '
        'and observations, and computes regional aggregates and model-versus-observation diagnostics, writes results'
        'to NetCDF, and produces timeseries plots.'
    )

    parser = argparse.ArgumentParser(description=code_description)

    # Mandatory arguments
    parser.add_argument('start_year', type=int, help='Start year of time range.')
    parser.add_argument('start_month', type=int, help='Start month of time range.')
    parser.add_argument('start_day', type=int, help='Start day of time range.')
    parser.add_argument('end_year', type=int, help='End year of time range.')
    parser.add_argument('end_month', type=int, help='End month of time range.')
    parser.add_argument('end_day', type=int, help='End day of time range.')
    parser.add_argument('forecast_init_hour', type=int, help='HRRR forecast initialization hour.')
    parser.add_argument('forecast_lead_hour', type=int, help='HRRR forecast lead time in hours.')
    parser.add_argument(
        'path_to_geojson',
        type=str,
        help='Path to RTO/ISO geometries file (available from https://atlas.eia.gov/datasets/rto-regions)',
    )
    parser.add_argument('isdlite_data_dir', type=str, help='Directory where ISDLite data are located.')
    parser.add_argument('hrrr_data_dir', type=str, help='Directory where HRRR data are located.')
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
    path_to_geojson = Path(args.path_to_geojson)
    isdlite_data_dir = Path(args.isdlite_data_dir)
    hrrr_data_dir = Path(args.hrrr_data_dir)
    out_dir = Path(args.out_dir)

    return (
        start_date,
        end_date,
        forecast_init_hour,
        forecast_lead_hour,
        path_to_geojson,
        isdlite_data_dir,
        hrrr_data_dir,
        out_dir,
    )


(start_date, end_date, forecast_init_hour, forecast_lead_hour, path_to_geojson, isdlite_data_dir, hrrr_data_dir, out_dir) = (
    arg_parse(sys.argv[1:])
)

#
# Variables to compare
#

hrrr_variable = 'TMP_P0_L103_GLC0'
hrrr_lat_name = 'gridlat_0'
hrrr_lon_name = 'gridlon_0'

isdlite_variable = 'T'
isd_lat_name = 'LAT'
isd_lon_name = 'LON'

#
# Loop over RTO / ISO regions
#

region_names = ['ERCOT', 'CAISO', 'ISONE', 'NYISO', 'MISO', 'SPP', 'PJM']

dfs_hrrr_regional_statistics = []
dfs_hrrr_vs_isdlite_statistics = []

region_name_last = None

for region_name in region_names:
    # Read RTO/ISO region geometry

    region_gdf = rto_iso.region(path_to_geojson, region_name)

    # Load ISDLite observations

    isdlite_file_name = region_name + '.' + str(start_date.year) + '-' + str(end_date.year) + '_ISD_Lite_observations.nc'

    isdlite_file = isdlite_data_dir / isdlite_file_name

    isdlite_stations = stations.Stations.from_netcdf(isdlite_file)

    isdlite_data_units = isdlite_stations.observations[isdlite_variable].attrs['units']

    # Convert T to Kelvin
    if isdlite_data_units == 'C' or isdlite_data_units == 'degC':
        isdlite_stations.observations[isdlite_variable] = isdlite_stations.observations[isdlite_variable] + 273.15
        isdlite_stations.observations[isdlite_variable].attrs['units'] = 'K'
        isdlite_data_units = 'K'

    #
    # Loop over HRRR forecasts
    #

    # Initialize empty HRRR regional statistics object

    hrrr_regional_statistics = statistics.Model2DRegionalStatisticsTimeSeries(None, None, None, None, None)

    # Initialize empty HRRR vs ISDLite evaluation statistics object

    hrrr_vs_isdlite_statistics = statistics.ModelVsObs2DStatisticsTimeSeries(None, None, None, None, None, None, None, None)

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
            region_name
            + ' region, forecast init time '
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

        hrrr_file_name = (
            'hrrr.t' + str(forecast_init_hour).zfill(2) + 'z.wrfsfcf' + str(forecast_lead_hour).zfill(2) + '_select_vars.nc'
        )

        hrrr_file = hrrr_data_dir / hrrr_dir_name / 'conus' / hrrr_file_name

        # Open HRRR file and load data
        hrrr_ds = xr.open_dataset(hrrr_file)

        # Units

        hrrr_data_units = hrrr_ds[hrrr_variable].attrs['units']

        assert hrrr_data_units == isdlite_data_units, 'Units mismatch between HRRR forecast and ISDLite observations.'

        # Convert T to Kelvin
        if hrrr_data_units == 'C' or hrrr_data_units == 'degC':
            hrrr_ds[hrrr_variable] = hrrr_ds[hrrr_variable] + 273.15
            hrrr_ds[hrrr_variable].attrs['units'] = 'K'
            hrrr_data_units = 'K'

        #
        # Statistics of HRRR forecast in region
        #

        if region_name != region_name_last:
            run_region_change_test = True
            region_name_last = region_name
        else:
            run_region_change_test = False

        hrrr_regional_statistics.extend(
            statistics.Model2DRegionalStatisticsTimeSeries(
                hrrr_ds, hrrr_variable, hrrr_lat_name, hrrr_lon_name, region_gdf, run_region_change_test=run_region_change_test
            )
        )

        #
        # Statistics of HRRR forecast vs ISDLite observations in region
        #

        hrrr_vs_isdlite_statistics.extend(
            statistics.ModelVsObs2DStatisticsTimeSeries(
                hrrr_ds,
                hrrr_variable,
                hrrr_lat_name,
                hrrr_lon_name,
                isdlite_stations,
                isdlite_variable,
                isd_lat_name,
                isd_lon_name,
            )
        )

        # Close HRRR file

        hrrr_ds.close()

        # Iterate date

        forecast_creation_date += time_step

    # Save results

    hrrr_regional_statistics.write_ds2netcdf(hrrr_regional_statistics.DataSet(), out_dir)
    hrrr_vs_isdlite_statistics.write_ds2netcdf(hrrr_vs_isdlite_statistics.DataSet(), out_dir)

    # Construct/append dataframes

    dfs_hrrr_regional_statistics.append(hrrr_regional_statistics.DataFrame())
    dfs_hrrr_vs_isdlite_statistics.append(hrrr_vs_isdlite_statistics.DataFrame())

#
# Plotting
#

# Each RTO/ISO individually (one region per plot)

for region_name, df in zip(region_names, dfs_hrrr_regional_statistics, strict=False):
    plot_dir = out_dir / 'plots' / region_name

    plot_file_paths = plotting.plot_df_timeseries(plot_dir, [df], [region_name])

for region_name, df in zip(region_names, dfs_hrrr_vs_isdlite_statistics, strict=False):
    plot_dir = out_dir / 'plots' / region_name

    plot_file_paths = plotting.plot_df_timeseries(plot_dir, [df], [region_name])

plot_dir = out_dir / 'plots' / '_'.join(region_names)

# All RTO/ISO together (all regions per plot)

plot_file_paths = plotting.plot_df_timeseries(plot_dir, dfs_hrrr_regional_statistics, region_names)
plot_file_paths = plotting.plot_df_timeseries(plot_dir, dfs_hrrr_vs_isdlite_statistics, region_names)
