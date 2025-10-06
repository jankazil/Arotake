from pathlib import Path

import pandas as pd
import xarray as xr


class TimeSeries:
    """
    Represents a time series.

    Functionalities:

      - Loads time series from a NetCDF file when initialized.

    """

    def __init__(self, file_path: Path):
        """
        Initializes the TimeSeries object with data from a NetCDF file.

        After loading, a new variable `'UTC'` is added to the dataset containing
        timezone-aware datetime objects representing observation times in UTC.

        Args:
            file_path (Path): Path to the NetCDF file containing station observations.
        """

        # Load the NetCDF dataset into an xarray Dataset
        self.time_series = xr.load_dataset(file_path)

        # Convert the 'time' coordinate to pandas datetime objects
        times = pd.to_datetime(self.time_series.coords['time'].values)

        # Ensure all times are timezone-aware and converted to UTC
        if times.tz is None:
            times = times.tz_localize('UTC')
        else:
            times = times.tz_convert('UTC')

        # Add a new variable 'UTC' containing the UTC timestamps as Python datetime objects
        self.time_series['UTC'] = ('time', times.to_pydatetime())

        return
