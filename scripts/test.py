from pathlib import Path

from isd_lite_data import stations

# Load IDS Lite observations

data_dir = Path('..') / 'data' / 'ISD-LITE' / 'data'

nc_file_name = 'tx_stations.2021-2025.nc'

nc_file_path = data_dir / nc_file_name

tx_stations = stations.Stations.from_netcdf(nc_file_path)
