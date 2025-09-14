#!/usr/bin/env python

'''
Analyze HRRR forecasts versus ISDLite observations by RTO/ISO region.

This script compares 2-D HRRR forecast fields against ISDLite point
observations across specified U.S. RTO/ISO regions (ERCOT, CAISO, SPP, etc.).
For each region and forecast time, the script:

1. Loads the polygonal region geometry from a GeoJSON file.
2. Loads ISDLite observations for stations in that region from a NetCDF file.
   - Converts observation units to Kelvin if stored in Celsius/degC.
3. Iterates through daily HRRR forecasts for the analysis period:
   - Constructs the HRRR file path based on initialization time and lead time.
   - Opens the HRRR dataset and checks units against ISDLite observations.
   - Converts HRRR variable to Kelvin if needed.
   - Computes HRRR-only regional statistics (mean, std, variance) over the region
     using Model2DRegionalStatisticsTimeSeries.
   - Computes HRRR vs ISDLite statistics (bias, RMSE, correlation, coverage) at
     station locations using ModelVsObs2DStatisticsTimeSeries.
   - Extends the accumulating time series containers with results for each forecast.
   - Closes the HRRR dataset before proceeding to the next day.
4. After looping through all forecasts, converts accumulated results into pandas
   DataFrames for both HRRR-only and HRRR-vs-ISDLite statistics.
5. Produces plots of all statistics versus time using plotting.plot_df_timeseries,
   saving one figure per statistic in a region-specific directory under ../plots/.
6. Saves the statistics time series as netCDF files '../data/RTO_ISO_regions_timeseries.

Configuration
-------------

- Forecast lead time and forecast hour are set near the top of the script
- Analysis period is controlled by start_time and end_time.
- HRRR input files are expected under ../data/HRRR/data/YYYYMMDD/conus/
  with filenames of the form
    hrrr.t{HH}z.wrfsfcf{LL}_select_vars.nc
  where HH = initialization hour, LL = lead time in hours.
- ISDLite input files are expected under
  ../data/RTO_ISO_regions_ISD_stations/{REGION}.{YEARS}_ISD_Lite_observations.nc
- Region geometries are read from ../data/RTO_ISO_regions.geojson.

Outputs
-------

- dfs_hrrr_regional_statistics: list of pandas.DataFrame objects with HRRR-only
  regional statistics for each region.
- dfs_hrrr_vs_isdlite_statistics: list of pandas.DataFrame objects with HRRR vs
  ISDLite evaluation metrics for each region.
- Plots: PNG files saved under ../plots/<region_names_combined>/, one per statistic,
  showing timeseries overlays across the analyzed regions.
'''

from datetime import datetime, timedelta, timezone
from pathlib import Path

import xarray as xr
from isd_lite_data import stations

from arotake import plotting, rto_iso, statistics

#
# Keep attributes of xarray Datasets upon numerical operations (otherwise they are lost)
#

xr.set_options(keep_attrs=True)

#
# HRRR forecast time parameters
#

forecast_hour = 20
forecast_lead_time_hours = 32

#
# Period to analyze
#

start_time = datetime(2020, 12, 3, forecast_hour, 0, 0, tzinfo=timezone.utc)
end_time = datetime(2021, 12, 31, forecast_hour, 0, 0, tzinfo=timezone.utc)

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
    # Obtain RTO/ISO region geometry

    path_to_geojson = Path('') / '..' / 'data' / 'RTO_ISO_regions.geojson'

    region_gdf = rto_iso.region(path_to_geojson, region_name)

    # ISDLite observations

    isdlite_file_name = region_name + '.2020-2025_ISD_Lite_observations.nc'

    isdlite_file = Path('..') / 'data' / 'RTO_ISO_regions_ISD_stations' / isdlite_file_name

    isdlite_stations = stations.Stations.from_netcdf(isdlite_file)

    isdlite_data_units = isdlite_stations.observations[isdlite_variable].attrs['units']

    # Convert T to Kelvin
    if isdlite_data_units == 'C' or isdlite_data_units == 'degC':
        isdlite_stations.observations[isdlite_variable] = isdlite_stations.observations[isdlite_variable] + 273.15
        isdlite_stations.observations[isdlite_variable].attrs['units'] = 'K'
        isdlite_data_units = 'K'

    # Station names
    # isd_station_names = isdlite_stations.observations.data_vars['STATION_NAME'].values
    # isd_station_isd = isdlite_stations.observations.data_vars['STATION_ID'].values

    #
    # Loop over HRRR forecasts
    #

    # HRRR forecast timing

    forecast_lead_time_delta = timedelta(hours=forecast_lead_time_hours)

    time_step = timedelta(days=1)

    forecast_time = start_time

    # Initialize empty HRRR regional statistics object

    hrrr_regional_statistics = statistics.Model2DRegionalStatisticsTimeSeries(None, None, None, None, None)

    # Initialize empty HRRR vs ISDLite evaluation statistics object

    hrrr_vs_isdlite_statistics = statistics.ModelVsObs2DStatisticsTimeSeries(None, None, None, None, None, None, None, None)

    while forecast_time <= end_time:
        print(region_name + ' region, forecast time ' + str(forecast_time))

        #
        # HRRR data
        #

        # Forecast initilization time
        forecast_init_time = forecast_time - forecast_lead_time_delta

        # Construct path of HRRR file (in HRRR data directory)
        hrrr_data_dir = Path('..') / 'data' / 'HRRR' / 'data'
        hrrr_dir_name = (
            'hrrr.'
            + str(forecast_init_time.year).zfill(4)
            + str(forecast_init_time.month).zfill(2)
            + str(forecast_init_time.day).zfill(2)
        )
        hrrr_file_name = (
            'hrrr.t'
            + str(forecast_init_time.hour).zfill(2)
            + 'z.wrfsfcf'
            + str(forecast_lead_time_hours).zfill(2)
            + '_select_vars.nc'
        )
        hrrr_file = hrrr_data_dir / hrrr_dir_name / 'conus' / hrrr_file_name

        # Open HRRR file
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

        # Iterate forecast time

        forecast_time += time_step

    # Save results

    data_dir = Path('..') / 'data' / 'RTO_ISO_regions_timeseries'
    hrrr_regional_statistics.write_ds2netcdf(hrrr_regional_statistics.DataSet(), data_dir)
    hrrr_vs_isdlite_statistics.write_ds2netcdf(hrrr_vs_isdlite_statistics.DataSet(), data_dir)

    # Construct/append dataframes

    dfs_hrrr_regional_statistics.append(hrrr_regional_statistics.DataFrame())
    dfs_hrrr_vs_isdlite_statistics.append(hrrr_vs_isdlite_statistics.DataFrame())

#
# Plotting
#

# Each RTO/ISO individually (one region per plot)

for region_name, df in zip(region_names, dfs_hrrr_regional_statistics, strict=False):
    plot_dir = Path('..') / 'plots' / region_name

    plot_file_paths = plotting.plot_df_timeseries(plot_dir, [df], [region_name])

for region_name, df in zip(region_names, dfs_hrrr_vs_isdlite_statistics, strict=False):
    plot_dir = Path('..') / 'plots' / region_name

    plot_file_paths = plotting.plot_df_timeseries(plot_dir, [df], [region_name])

plot_dir = Path('..') / 'plots' / '_'.join(region_names)

# All RTO/ISO together (all regions per plot)

plot_file_paths = plotting.plot_df_timeseries(plot_dir, dfs_hrrr_regional_statistics, region_names)
plot_file_paths = plotting.plot_df_timeseries(plot_dir, dfs_hrrr_vs_isdlite_statistics, region_names)
