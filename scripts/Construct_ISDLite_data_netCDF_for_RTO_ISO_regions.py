#!/usr/bin/env python

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib.lines as lines
import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LatitudeFormatter, LatitudeLocator, LongitudeFormatter, LongitudeLocator
from isd_lite_data import ncei, stations
from shapely.geometry import Point

from arotake import rto_iso

#
# Command line arguments
#


def arg_parse(argv=None):
    '''
    Argument parser which returns the parsed values given as arguments.
    '''

    code_description = (
        'Download, filter, and visualize NOAA ISD-Lite station data for a specified date range, '
        'grouping stations by Regional Transmission Organization (RTO) / Independent System Operator (ISO) regions. '
        'Given a start and end date, a path to the RTO/ISO GeoJSON, and a target data directory, '
        'the script retrieves station metadata, filters stations to the contiguous U.S., checks data availability, '
        'associates stations with RTO/ISO regions, saves region-specific station lists, and plots station locations '
        'on a map using Cartopy. It then downloads the ISD-Lite observations (optionally in parallel), loads them '
        'into memory, assigns region attributes, and writes the results as region-specific netCDF files.'
    )

    parser = argparse.ArgumentParser(description=code_description)

    # Mandatory arguments
    parser.add_argument('start_year', type=int, help='Start year of time range.')
    parser.add_argument('start_month', type=int, help='Start month of time range.')
    parser.add_argument('start_day', type=int, help='Start day of time range.')
    parser.add_argument('end_year', type=int, help='End year of time range.')
    parser.add_argument('end_month', type=int, help='End month of time range.')
    parser.add_argument('end_day', type=int, help='End day of time range.')
    parser.add_argument(
        'path_to_geojson',
        type=str,
        help='Path to RTO/ISO geometries file (available from https://atlas.eia.gov/datasets/rto-regions)',
    )
    parser.add_argument(
        'isdlite_data_dir',
        type=str,
        help='Directory where ISDLite data are located/will be downloaded to. Will be created if it does not exist.',
    )

    # Optional arguments
    parser.add_argument(
        '-n',
        '--n',
        type=int,
        help='Number of parallel download processes. n > 1 accelerates downloads significantly, but can result in network errors or in the server refusing to cooperate. In case of errors, set to 1',
    )

    args = parser.parse_args()

    start_date = datetime(year=args.start_year, month=args.start_month, day=args.start_day, tzinfo=timezone.utc)
    end_date = datetime(year=args.end_year, month=args.end_month, day=args.end_day, tzinfo=timezone.utc)
    path_to_geojson = Path(args.path_to_geojson)
    isdlite_data_dir = Path(args.isdlite_data_dir)
    n_jobs = args.n

    return (start_date, end_date, path_to_geojson, isdlite_data_dir, n_jobs)


(start_date, end_date, path_to_geojson, isdlite_data_dir, n_jobs) = arg_parse(sys.argv[1:])

#
# Read RTO/ISO geometries from file
#

regions_gdf = rto_iso.regions(path_to_geojson)

#
# Download ISDLite station meta data
#

isdlite_history_file = Path(isdlite_data_dir) / 'isd-history.txt'

# ncei.download_file(isd_lite_stations_url,isdlite_history_file, refresh = True, verbose = True)

ncei.download_stations(isdlite_history_file)

# Create stations object
ISD_stations = stations.Stations.from_file(isdlite_history_file)

# Filter by coordinates - for the contiguous US - this will reduce the number of stations to process

min_lat = 24
max_lat = 50
min_lon = -125
max_lon = -65

ISD_stations = ISD_stations.filter_by_coordinates(min_lat, max_lat, min_lon, max_lon)

# The following might take a little time as we will check whether data files are actually available for download
ISD_stations = ISD_stations.filter_by_data_availability(start_date, end_date, verbose=True)

# Create pandas dataframe with stations
isd_stations_df = ISD_stations.meta_data

# Create geopandas dataframe with stations
isd_stations_gdf = gpd.GeoDataFrame(
    isd_stations_df,
    geometry=[Point(coordinates) for coordinates in zip(isd_stations_df.LON, isd_stations_df.LAT, strict=False)],
    crs="EPSG:4326",  # lon/lat in degrees
)

# Identify which RTO/ISO region each ISD station is in

isd_stations_regions_gdf = gpd.sjoin(
    isd_stations_gdf, regions_gdf[["name", "geometry"]], how="left", predicate="covered_by"
).rename(columns={"name": "rto_iso_region"})

# Drop the index column in the result

if "index_right" in isd_stations_regions_gdf.columns:
    isd_stations_regions_gdf = isd_stations_regions_gdf.drop(columns="index_right")

# RTO/ISO regions

# for region in isd_stations_regions_gdf['rto_iso_region'].unique():
#    print(region)

# Define RTO/ISO regions as list with a specific order for convenience

regions = ['MISO', 'PJM', 'ERCOT', 'CAISO', 'SPP', 'ISONE', 'NYISO']

# Save the station metadata for each RTO/ISO region in files in ISD station "history" file format

for region in regions:
    mask = isd_stations_regions_gdf['rto_iso_region'] == region

    isd_stations_meta_data = stations.Stations(isd_stations_regions_gdf[mask])

    isd_stations_meta_data_file_path = region + '.' + str(start_date.year) + '-' + str(end_date.year) + '_ISD_stations.txt'

    isd_stations_meta_data_file_path = isdlite_data_dir / isd_stations_meta_data_file_path

    isd_stations_meta_data_file_path.parent.mkdir(parents=True, exist_ok=True)

    isd_stations_meta_data.save(
        'Integrated Surface Dataset (ISD) stations with observations in the period '
        + start_date.isoformat()
        + ' - '
        + end_date.isoformat()
        + ' in the RTO/ISO region '
        + region,
        isd_stations_meta_data_file_path,
    )

#
# Plot ISD stations in the RTO/ISO regions
#

regions_gdf = regions_gdf.to_crs('EPSG:4326')

# Define the target map projection (Lambert Conformal as in HRRR)
hrrr_proj = ccrs.LambertConformal(
    central_longitude=-97.5,
    central_latitude=38.5,
    standard_parallels=(38.5, 38.5),
)

# Crate plot - here, passing a Cartopy projection will make 'ax' being a cartopy.mpl.geoaxes.GeoAxes object, and Cartopy methods can be used with 'ax'

fig, ax = plt.subplots(figsize=(18, 9), subplot_kw={"projection": hrrr_proj})

# Set the map extent (using the PlateCarree coordinate reference system which just means lon/lat)
ax.set_extent([-123, -71, 25, 50], crs=ccrs.PlateCarree())

# Coastlines, national borders, state lines

ax.coastlines(resolution="50m", linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.25)
ax.add_feature(cfeature.STATES, linewidth=0.25)

# Colors for regions
colors = {
    'MISO': '#9467bd',
    'PJM': '#8c564b',
    'ERCOT': '#0000ff',
    'CAISO': '#1f77b4',
    'SPP': '#ff0000',
    'ISONE': '#2ca02c',
    'NYISO': '#ff7f0e',
}

# Size of location markers
marker_sizes = {
    "MISO": 5,
    "PJM": 3,
    "ERCOT": 3,
    "CAISO": 2,
    "SPP": 1,
    "ISONE": 3,
    "NYISO": 2,
}

# Draw locations

for region in regions:
    mask = isd_stations_regions_gdf['rto_iso_region'] == region

    ax.plot(
        isd_stations_regions_gdf.loc[mask, 'LON'],
        isd_stations_regions_gdf.loc[mask, 'LAT'],
        'o',
        markersize=marker_sizes.get(region, 2),
        color=colors.get(region, '#000000'),
        transform=ccrs.PlateCarree(),  # tell Cartopy these are lat/lon coords
    )

# Coordinate gridlines

gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='black', alpha=1.0, linestyle='--')

# Grid labels on the axes (not inline on the map)

gl.x_inline = False
gl.y_inline = False
gl.bottom_labels = True
gl.left_labels = True
gl.top_labels = False
gl.right_labels = False

# Choose where lines go

gl.xlocator = LongitudeLocator(nbins=6)
gl.ylocator = LatitudeLocator(nbins=6)

# Nice degree formatting

gl.xformatter = LongitudeFormatter(number_format='.0f', degree_symbol='°')
gl.yformatter = LatitudeFormatter(number_format='.0f', degree_symbol='°')

# Force horizontal labels

gl.xlabel_style = {'rotation': 0, 'ha': 'center', 'va': 'top', 'size': 12}
gl.ylabel_style = {'rotation': 0, 'ha': 'right', 'va': 'center', 'size': 12}

# Build a legend

handles = [
    lines.Line2D(
        [],
        [],
        color=colors[key],
        marker='o',
        linestyle='None',  # no connecting line
        markersize=marker_sizes[key],
        label=key,
    )
    for key in regions
    if key in colors
]

legend = ax.legend(handles=handles, loc="center right")

title = (
    'Integrated Surface Dataset (ISD) stations by Regional Transmission Organizations (RTOs) / Independend System Operators (ISOs)\n'
    + 'Stations shown have observations in the period '
    + start_date.isoformat()
    + ' - '
    + end_date.isoformat()
    + '\n'
    + 'ISD stations: National Centers for Environmental Information (NCEI), https://www.ncei.noaa.gov (accessed '
    + datetime.today().date().strftime("%Y-%m-%d")
    + ')\n'
    'RTO/ISO regions: https://atlas.eia.gov/datasets/rto-regions'
)

title_object = ax.set_title(title, fontsize=12)

plot_path = isdlite_data_dir / 'plots' / 'RTO_ISO_regions_ISD_stations_map.png'

plot_path.parent.mkdir(parents=True, exist_ok=True)

path = fig.savefig(plot_path, bbox_inches="tight", dpi=600)

#
# Download ISDLite data and save them as netCDF files
#

# Directory where the ISD station history files for the RTO/ISO regions are located,
# and where the netCDF files with the ISDLite observations for these regions will be saved

for region in regions:
    # Load stations in the region for the period of interest

    isd_stations_meta_data_file_path = region + '.' + str(start_date.year) + '-' + str(end_date.year) + '_ISD_stations.txt'

    isd_stations_meta_data_file_path = isdlite_data_dir / isd_stations_meta_data_file_path

    region_stations = stations.Stations.from_file(isd_stations_meta_data_file_path)

    # Download IDSLite station observations - this will download the observations
    # files only if they have changed online, e.g. the current year's files that
    # have been updated with observations since the last download.

    local_files = ncei.download_many(
        start_date.year,
        end_date.year,
        region_stations.id(),
        isdlite_data_dir,
        n_jobs=n_jobs,
        refresh=False,
        verbose=True,
    )

    # Load the observations from the local files

    region_stations.load_observations(isdlite_data_dir, start_date.year, end_date.year, verbose=True)

    # Set region global attribute

    region_stations.observations.attrs['region'] = region

    # Save observations as netCDF file

    nc_file_name = region + '.' + str(start_date.year) + '-' + str(end_date.year) + '_ISD_Lite_observations.nc'

    nc_file_path = isdlite_data_dir / nc_file_name

    region_stations.write_observations2netcdf(nc_file_path)
