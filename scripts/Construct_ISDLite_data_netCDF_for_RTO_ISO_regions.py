#!/usr/bin/env python

'''
Construct ISD-Lite station metadata, a station map, and per-region ISD-Lite
netCDF files for U.S. RTO/ISO regions over a configurable period.

Overview
--------
This script performs an end-to-end workflow to:
1) Load RTO/ISO geographic regions from GeoJSON.
2) Load the NCEI ISD station "history" metadata, restrict to a bounding box
   over the CONUS, and further filter stations by observed data availability
   within the requested time window.
3) Spatially join stations to their containing RTO/ISO region.
4) For each target region, write a region-specific ISD "history" file that
   preserves the ISD text format.
5) Render and save a Cartopy map of all retained stations, colored by region.
6) For each region, download ISD-Lite observation files for all retained
   stations for the requested years (skipping unchanged local files),
   load observations, attach region metadata, and save a consolidated
   netCDF file of ISD-Lite observations.

Inputs (on disk)
----------------

- RTO/ISO regions GeoJSON:

    ../data/RTO_ISO_regions.geojson

  This Regional Transmission Organization / Independent System Operators map can be obtained from

    https://atlas.eia.gov/datasets/rto-regions

- ISD station "history" file (text):

    ../data/ISD-LITE/data/isd-history.txt

  This can be downloaded with https://github.com/jankazil/isd-lite-data/blob/main/demos/demo_ncei_download_stations.py

- Local ISD-Lite data directory

    ../data/ISD-LITE/data/

  These data can be downloaded with the "CONUS" scripts in https://github.com/jankazil/isd-lite-data/tree/main/scripts

Key Configuration (edit in code)
--------------------------------

- Time window:
    start_year, end_year
    start_date = datetime(start_year, 12, 2)   # HRRR v4 operational start
    end_date   = datetime(end_year, 6, 30, 23, 59, 59)
- Spatial filter (CONUS bounding box):
    min_lat = 24, max_lat = 50, min_lon = -125, max_lon = -65
- Regions to process and their display order:
    regions = ['MISO', 'PJM', 'ERCOT', 'CAISO', 'SPP', 'ISONE', 'NYISO']

Outputs (on disk)
-----------------

- Region-specific ISD station metadata files (ISD "history" text format):
    ../data/RTO_ISO_regions_ISD_stations/{REGION}.{YYYY}-{YYYY}_ISD_stations.txt
- Station map PNG (Cartopy Lambert Conformal projection approximating HRRR):
    ../data/RTO_ISO_regions_ISD_stations/RTO_ISO_regions_ISD_stations_map.png
- Region-specific ISD-Lite observations in netCDF:
    ../data/RTO_ISO_regions_ISD_stations/{REGION}.{YYYY}-{YYYY}_ISD_Lite_observations.nc

Data Processing Details
-----------------------

- Station metadata are loaded via isd_lite_data.stations.Stations.from_file(...)
  and filtered by:
    a) geographic bounding box (EPSG:4326 coordinates)
    b) observed data availability within [start_date, end_date]
- Regions are assigned with a GeoPandas spatial join using predicate "covered_by".
- Map rendering uses Cartopy with coastlines, national borders, U.S. states,
  labeled graticules, and a legend keyed by region with per-region marker sizes.
- ISD-Lite files are downloaded from NCEI using ncei.download_many(...) with
  change detection (skips unchanged files). Observations are then loaded into
  each region's Stations object and written to netCDF via
  write_observations2netcdf(...). A global attribute "region" is set on the
  output.

Assumptions and Conventions
---------------------------

- Coordinate reference system for vector data is EPSG:4326 (lat/lon degrees).
- Region geometries are expected to cover stations ("covered_by" predicate).
- The time window is inclusive of station-level availability filtering and is
  year-spanning for downloads (start_year..end_year).
- File system layout matches the relative paths referenced above.

Notes
-----

- Filtering by data availability can be time-consuming; n_jobs is configurable.
- The map title documents access dates for the ISD and RTO/ISO region sources.
- Output folders are created if they do not exist.
'''

from datetime import datetime
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
# Time period - HRRR v4 became operational December 2 2020
#

start_year = 2020
end_year = 2025

start_date = datetime(start_year, 12, 2)
end_date = datetime(end_year, 6, 30, 23, 59, 59)

#
# RTO/ISO geometries file
#

path_to_geojson = Path('') / '..' / 'data' / 'RTO_ISO_regions.geojson'

regions_gdf = rto_iso.regions(path_to_geojson)

#
# ISDLite station meta data
#

# Directory where ISDLite data is located/will be downloaded to
isdlite_data_dir = Path('..') / 'data' / 'ISD-LITE' / 'data'

# Open ISDLite "history" (station metadata) file
isdlite_history_file = isdlite_data_dir / 'isd-history.txt'

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

isd_rto_iso_dir = Path('..') / 'data' / 'RTO_ISO_regions_ISD_stations'

for region in regions:
    mask = isd_stations_regions_gdf['rto_iso_region'] == region

    isd_stations_meta_data = stations.Stations(isd_stations_regions_gdf[mask])

    isd_stations_meta_data_file_path = region + '.' + str(start_date.year) + '-' + str(end_date.year) + '_ISD_stations.txt'

    isd_stations_meta_data_file_path = isd_rto_iso_dir / isd_stations_meta_data_file_path

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
    + 'ISD stations: National Centers for Environmental Information (NCEI), https://www.ncei.noaa.gov (accessed 2025-08-15)\n'
    'RTO/ISO regions: https://atlas.eia.gov/datasets/rto-regions (accessed 2025-08-26)'
)

title_object = ax.set_title(title, fontsize=12)

plot_path = Path('..') / 'data' / 'RTO_ISO_regions_ISD_stations' / 'RTO_ISO_regions_ISD_stations_map.png'

plot_path.parent.mkdir(parents=True, exist_ok=True)

path = fig.savefig(plot_path, bbox_inches="tight", dpi=600)

#
# Download ISDLite data and save them as netCDF files
#

# Directory where the ISD station history files for the RTO/ISO regions are located,
# and where the netCDF files with the ISDLite observations for these regions will be saved

rto_iso_isdlite_data_dir = Path('..') / 'data' / 'RTO_ISO_regions_ISD_stations'

for region in regions:
    # Load stations in the region for the period of interest

    isd_stations_meta_data_file_path = region + '.' + str(start_date.year) + '-' + str(end_date.year) + '_ISD_stations.txt'

    isd_stations_meta_data_file_path = isd_rto_iso_dir / isd_stations_meta_data_file_path

    region_stations = stations.Stations.from_file(isd_stations_meta_data_file_path)

    # Download IDSLite station observations - this will download the observations
    # files only if they have changed online, e.g. the current year's files that
    # have been updated with observations since the last download.

    local_files = ncei.download_many(
        start_date.year,
        end_date.year,
        region_stations.id(),
        isdlite_data_dir,
        n_jobs=32,
        refresh=False,
        verbose=True,
    )

    # Load the observations from the local files

    region_stations.load_observations(isdlite_data_dir, start_year, end_year, verbose=True)

    # Set region global attribute

    region_stations.observations.attrs['region'] = region

    # Save observations as netCDF file

    nc_file_name = region + '.' + str(start_date.year) + '-' + str(end_date.year) + '_ISD_Lite_observations.nc'

    nc_file_path = rto_iso_isdlite_data_dir / nc_file_name

    region_stations.write_observations2netcdf(nc_file_path)
