'''
Statistics utilities for 2-D gridded model data and point observations.

This module provides two time-series containers:

- Model2DRegionalStatisticsTimeSeries:
  Computes single-timestep regional statistics of a 2-D model field over a
  polygonal region (union of all polygons in a GeoDataFrame), caches region and
  grid to avoid recomputation, and stores results as extendable list-derived
  time series. Also exports results to a pandas.DataFrame with per-column metadata.

- ModelVsObs2DStatisticsTimeSeries:
  Evaluates a 2-D model field against point observations at a single valid time by
  interpolating the model to observation locations and computing summary and
  comparison metrics. Stores results as extendable list-derived time series and
  exports them to a pandas.DataFrame with per-column metadata.

Design notes
------------
- Time series are represented by lightweight list subclasses (TimeList) to avoid
  repeated DataFrame.append() costs during incremental accumulation.
- For model-only regional stats, the unioned region mask and model grid are cached
  at the class level to accelerate subsequent instances using the same inputs.
- Per-series metadata (units, description, var_name) are carried alongside values
  and are propagated into DataFrame Series.attrs when exporting.

Dependencies and expectations
-----------------------------
- xarray model datasets must provide:
  * A 2-D data variable (y, x) with attributes:
      'initial_time' (str, "%m/%d/%Y (%H:%M)", UTC),
      'forecast_time_units' == "hours",
      'forecast_time' convertible to int,
      'long_name', 'units'.
  * 2-D latitude/longitude variables (y, x) named by caller.
  * Optional dataset attribute 'model' for labeling.
- Observations containers must expose `observations` with:
  * Coordinate/variable 'UTC' containing timezone-aware datetimes in UTC.
  * 1-D arrays of lat/lon at the same index dimension as observation values.
  * A 2-D variable (time, station) for the evaluated quantity, NaN where missing.
- Spatial masking uses `regionmask.from_geopandas` on the unary union of geometries
  in the provided GeoDataFrame. CRS is aligned to the union GeoDataFrame when needed.
'''

from datetime import datetime, timedelta, timezone
from pathlib import Path

import geopandas as gp
import numpy as np
import pandas as pd
import regionmask
import xarray as xr

from arotake import interpolation


class TimeList(list):
    '''
    TimeList is a list that supports attributes (unlike lists). Inherits from list (obviously).
    '''

    pass


class Model2DRegionalStatisticsTimeSeries:
    '''
    Compute and hold a time series of regional statistics for a 2-D model field.

    Each instance adds one time step consisting of summary statistics of the
    selected model variable over the union of all polygons in a GeoDataFrame.
    Values and their metadata are stored in list-derived series that can be
    concatenated by calling "extend" with another instance.

    Caching
    -------
    Class-level caches speed up repeated usage with unchanged inputs:
    - _cached_model_lats, _cached_model_lons: the last model grid (lat/lon).
    - _cached_gdf: the last unioned region GeoDataFrame.
    - _cached_model_region_index: numpy indices selecting grid cells inside the region.

    The region index is recomputed when:
    - The model grid differs from the cached grid (array-equality check), or
    - The region (after unary union) is topologically different from the cached region,
      or no region has been cached yet. You may skip the topology test by setting
      "run_region_change_test=False" to save time when you know the region is unchanged.
    '''

    #
    # Class attributes - hold cached values to accelerate computation when
    # neither the model grid nor the region changes
    #

    # Cached model grid
    _cached_model_lats = None
    _cached_model_lons = None

    # Cached region
    _cached_gdf = None

    # Index array with indices in model data arrays representing the region
    _cached_model_region_index = None

    def __init__(
        self,
        model_ds: xr.Dataset,
        model_variable: str,
        model_lat_name: str,
        model_lon_name: str,
        region_gdf: gp.GeoDataFrame,
        run_region_change_test: bool = True,
    ):
        '''
        Initialize the object with one timestep of regional statistics.

        Behavior
        --------
        - Builds a single geometry by unary-union of all polygons in "region_gdf".
        - Masks the 2-D model grid with that union using "regionmask" and extracts
          all in-region values of "model_variable" to compute:
            mean, standard deviation, and variance.
        - Parses forecast initialization time (UTC), lead time (hours), and computes
          valid time from model variable attributes.
        - Stores values and per-series metadata in list-derived fields. If any
          required parameter is "None", initializes an empty container.

        Parameters
        ----------
        model_ds : xarray.Dataset
            Dataset holding the 2-D field, 2-D latitude and longitude, and required
            timing attributes on the target variable.
        model_variable : str
            Name of the model variable to summarize.
        model_lat_name : str
            Name of the latitude variable in "model_ds" (2-D, y×x).
        model_lon_name : str
            Name of the longitude variable in "model_ds" (2-D, y×x).
        region_gdf : geopandas.GeoDataFrame
            Polygons/multipolygons to be unioned and used as the regional mask.
        run_region_change_test : bool, default True
            If True, tests topological equality of the new unioned region against
            the cached one to decide whether to recompute indices. If False,
            assumes the region is unchanged and skips this test.

        Raises
        ------
        AssertionError
            If "model_variable" has "forecast_time_units" not equal to "hours".

        Notes
        -----
        - The cached region is stored in the CRS of the provided "region_gdf";
          if a cached region exists with a different CRS, it is reprojected for
          topological comparison.
        - Units for variance are recorded as "({units})^2".
        '''

        # If any of the parameters is none, initialize an empty object.

        if any(obj is None for obj in [model_ds, model_variable, model_lat_name, model_lon_name, region_gdf]):
            # Make sure to initialize lists for time series values and units

            self.model_name = []

            self.model_region_name = []

            self.model_variable = TimeList()
            self.model_variable_long_name = TimeList()

            self.model_forecast_init_time = TimeList()
            self.model_forecast_lead_time = TimeList()
            self.model_forecast_time = TimeList()

            self.model_forecast_init_time.units = []
            self.model_forecast_lead_time.units = []
            self.model_forecast_time.units = []

            self.model_forecast_init_time.var_name = None
            self.model_forecast_lead_time.var_name = None
            self.model_forecast_time.var_name = None

            self.model_forecast_init_time.description = None
            self.model_forecast_lead_time.description = None
            self.model_forecast_time.description = None

            self.model_region_mean = TimeList()
            self.model_region_std = TimeList()
            self.model_region_var = TimeList()

            self.model_region_mean.units = []
            self.model_region_std.units = []
            self.model_region_var.units = []

            self.model_region_mean.var_name = None
            self.model_region_std.var_name = None
            self.model_region_var.var_name = None

            self.model_region_mean.description = None
            self.model_region_std.description = None
            self.model_region_var.description = None

            return

        # Update cached model grid if necessary

        cached_model_lats = getattr(Model2DRegionalStatisticsTimeSeries, '_cached_model_lats', None)
        cached_model_lons = getattr(Model2DRegionalStatisticsTimeSeries, '_cached_model_lons', None)

        if (
            cached_model_lats is None
            or cached_model_lons is None
            or not np.array_equal(cached_model_lats, model_ds[model_lat_name])
            or not np.array_equal(cached_model_lons, model_ds[model_lon_name])
        ):
            Model2DRegionalStatisticsTimeSeries._cached_model_lats = model_ds[model_lat_name].copy(deep=True)
            Model2DRegionalStatisticsTimeSeries._cached_model_lons = model_ds[model_lon_name].copy(deep=True)
            model_grid_changed = True
        else:
            model_grid_changed = False

        # Update cached region if necessary

        cached_gdf = getattr(Model2DRegionalStatisticsTimeSeries, '_cached_gdf', None)

        if cached_gdf is None or run_region_change_test:
            # Build a single-geometry union from the provided region_gdf

            union_geom = region_gdf.unary_union
            union_gdf = gp.GeoDataFrame(
                {'name': ['-'.join(region_gdf['name'].tolist())]}, geometry=[union_geom], crs=region_gdf.crs
            )

            # Compare cached vs provided unioned geometry by topology

            if cached_gdf is not None:
                # Ensure same coordinate reference system first
                if union_gdf.crs != cached_gdf.crs:
                    cached_gdf = cached_gdf.to_crs(union_gdf.crs)

                # Topological equality (ignores vertex order and ring orientation)
                region_changed = not union_gdf.geometry.iloc[0].equals(cached_gdf.geometry.iloc[0])
            else:
                region_changed = True

        else:
            region_changed = False

        if region_changed:
            Model2DRegionalStatisticsTimeSeries._cached_gdf = union_gdf.copy(deep=True)

        # Recalculate the indices in the model grid representing the provided unioned region
        # and update cached values, if necessary

        if model_grid_changed or region_changed:
            # Build a mask for the unioned region in the model grid
            model_region_mask = regionmask.from_geopandas(union_gdf, names='name')

            # For a single unioned region, the mask is 2D (y, x)
            model_region_mask_da = model_region_mask.mask(
                lon_or_obj=model_ds[model_lon_name],
                lat=model_ds[model_lat_name],
            )

            # Update cached values
            Model2DRegionalStatisticsTimeSeries._cached_model_region_index = np.where(~np.isnan(model_region_mask_da))

        model_region_index = Model2DRegionalStatisticsTimeSeries._cached_model_region_index

        # Identify the model initialization time, forecast lead time, and forecast time

        self.model_forecast_init_time = datetime.strptime(model_ds[model_variable].attrs['initial_time'], '%m/%d/%Y (%H:%M)')
        self.model_forecast_init_time = self.model_forecast_init_time.replace(tzinfo=timezone.utc)

        assert model_ds[model_variable].attrs['forecast_time_units'] == 'hours', 'model forecast lead time units must be hours'

        self.model_forecast_lead_time_hours = int(model_ds[model_variable].attrs['forecast_time'])
        self.model_forecast_lead_time = timedelta(hours=self.model_forecast_lead_time_hours)

        self.model_forecast_time = self.model_forecast_init_time + self.model_forecast_lead_time

        # Model name

        if 'model' in model_ds.attrs:
            self.model_name = [model_ds.attrs['model']]
        else:
            self.model_name = ['model']

        # Calculate statistics over the unioned region in the model data

        self.model_variable = TimeList([model_variable])
        self.model_variable_long_name = TimeList([model_ds[model_variable].attrs['long_name']])

        values = model_ds[model_variable].values[model_region_index]

        self.model_region_mean = TimeList([np.mean(values)])
        self.model_region_std = TimeList([np.std(values)])
        self.model_region_var = TimeList([np.var(values)])

        self.model_region_mean.units = [model_ds[model_variable].attrs['units']]
        self.model_region_std.units = [model_ds[model_variable].attrs['units']]
        self.model_region_var.units = ['(' + model_ds[model_variable].attrs['units'] + ')^2']

        self.model_region_mean.var_name = 'model_mean'
        self.model_region_std.var_name = 'model_std'
        self.model_region_var.var_name = 'model_var'

        self.model_region_mean.description = self.model_name[0] + ' mean'
        self.model_region_std.description = self.model_name[0] + ' std. dev.'
        self.model_region_var.description = self.model_name[0] + ' variance'

        # Set other instance variables

        self.model_region_name = ['-'.join(Model2DRegionalStatisticsTimeSeries._cached_gdf['name'].tolist())]

        self.model_forecast_init_time = TimeList([self.model_forecast_init_time])
        self.model_forecast_lead_time = TimeList([self.model_forecast_lead_time])
        self.model_forecast_time = TimeList([self.model_forecast_time])

        self.model_forecast_init_time.units = ['UTC']
        self.model_forecast_lead_time.units = ['h']
        self.model_forecast_time.units = ['UTC']

        self.model_forecast_init_time.var_name = 'model_forecast_init_time'
        self.model_forecast_lead_time.var_name = 'model_forecast_lead_time'
        self.model_forecast_time.var_name = 'model_forecast_time'

        self.model_forecast_init_time.description = 'forecast initialization time (UTC)'
        self.model_forecast_lead_time.description = 'forecast lead time'
        self.model_forecast_time.description = 'forecast time'

        return

    def test_internal_consistency(self):
        '''
        Validate internal consistency across all accumulated timesteps.

        Checks that the region name, model variable identifiers, forecast times
        (by HH:MM:SS), lead times, and all recorded units are constant across
        the series. Raises AssertionError on inconsistency. Does nothing for
        an empty container.
        '''

        if self.model_forecast_init_time:
            assert len(set(self.model_name)) <= 1, 'Time series covers different models'

            assert len(set(self.model_region_name)) <= 1, 'Time series covers different regions'

            assert len(set(self.model_variable)) <= 1, 'Time series is internally inconsistent in the model variables'

            assert (
                len(set(self.model_variable_long_name)) <= 1
            ), 'Time series is internally inconsistent in the model variable long names'

            assert len(set([string.strftime("%H:%M:%S") for string in self.model_forecast_init_time[:]])) <= 1, (
                'Time series is internally inconsistent in the ' + self.model_forecast_init_time.description
            )

            assert len(set(self.model_forecast_lead_time)) <= 1, (
                'Time series is internally inconsistent in the ' + self.model_forecast_lead_time.description
            )

            assert len(set(self.model_forecast_init_time.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_forecast_init_time.description
            )
            assert len(set(self.model_forecast_time.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_forecast_time.description
            )
            assert len(set(self.model_forecast_lead_time.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_forecast_lead_time.description
            )

            assert len(set(self.model_region_mean.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_region_mean.description
            )
            assert len(set(self.model_region_std.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_region_std.description
            )
            assert len(set(self.model_region_var.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_region_var.description
            )

    def extend(self, other: 'Model2DRegionalStatisticsTimeSeries'):
        '''
        Append another instance's time series to this instance's time series.

        Extends all value, unit, description, and name lists in place from
        "other" and then validates internal consistency.

        Parameters
        ----------
        other : Model2DRegionalStatisticsTimeSeries
            Another container produced with the same configuration whose contents
            should be appended as the next timestep.
        '''

        # Append

        self.model_name.extend(other.model_name)
        self.model_region_name.extend(other.model_region_name)
        self.model_variable.extend(other.model_variable)
        self.model_variable_long_name.extend(other.model_variable_long_name)

        self.model_forecast_init_time.extend(other.model_forecast_init_time)
        self.model_forecast_lead_time.extend(other.model_forecast_lead_time)
        self.model_forecast_time.extend(other.model_forecast_time)

        self.model_forecast_init_time.units.extend(other.model_forecast_init_time.units)
        self.model_forecast_lead_time.units.extend(other.model_forecast_lead_time.units)
        self.model_forecast_time.units.extend(other.model_forecast_time.units)

        self.model_forecast_init_time.var_name = other.model_forecast_init_time.var_name
        self.model_forecast_lead_time.var_name = other.model_forecast_lead_time.var_name
        self.model_forecast_time.var_name = other.model_forecast_time.var_name

        self.model_forecast_init_time.description = other.model_forecast_init_time.description
        self.model_forecast_lead_time.description = other.model_forecast_lead_time.description
        self.model_forecast_time.description = other.model_forecast_time.description

        self.model_region_mean.extend(other.model_region_mean)
        self.model_region_std.extend(other.model_region_std)
        self.model_region_var.extend(other.model_region_var)

        self.model_region_mean.units.extend(other.model_region_mean.units)
        self.model_region_std.units.extend(other.model_region_std.units)
        self.model_region_var.units.extend(other.model_region_var.units)

        self.model_region_mean.var_name = other.model_region_mean.var_name
        self.model_region_std.var_name = other.model_region_std.var_name
        self.model_region_var.var_name = other.model_region_var.var_name

        self.model_region_mean.description = other.model_region_mean.description
        self.model_region_std.description = other.model_region_std.description
        self.model_region_var.description = other.model_region_var.description

        # Test if the time series are internally consistent

        self.test_internal_consistency()

        return

    def DataFrame(self) -> pd.DataFrame:
        '''
        Export the accumulated series to a pandas DataFrame.

        Columns
        -------
        - model_forecast_time (datetime, timezone-aware UTC)
        - model_forecast_init_time (datetime, timezone-aware UTC)
        - model_forecast_lead_time (datetime.timedelta)
        - model_mean / model_std / model_var (floats; regional stats)

        Metadata
        --------
        - Per-column: stored in "Series.attrs" as "{'units', 'description'}".
        - Global: "DataFrame.attrs" includes
          "model_variable", "model_variable_long_name", "model_region_name".

        Returns
        -------
        pandas.DataFrame
            A new DataFrame; arrays do not alias the internal lists.

        Raises
        ------
        AssertionError
            If internal consistency checks fail.
        '''

        # Test if the time series are internally consistent

        self.test_internal_consistency()

        # Create the data frame with one time series per column:

        df = pd.DataFrame(
            {
                self.model_forecast_time.var_name: self.model_forecast_time,
                self.model_forecast_init_time.var_name: self.model_forecast_init_time,
                self.model_forecast_lead_time.var_name: self.model_forecast_lead_time,
                self.model_region_mean.var_name: self.model_region_mean,
                self.model_region_std.var_name: self.model_region_std,
                self.model_region_var.var_name: self.model_region_var,
            }
        )

        # Store units and description as per-column attributes

        df[self.model_forecast_time.var_name].attrs = {
            'units': self.model_forecast_time.units[0],
            'description': self.model_forecast_time.description,
        }
        df[self.model_forecast_init_time.var_name].attrs = {
            'units': self.model_forecast_init_time.units[0],
            'description': self.model_forecast_init_time.description,
        }
        df[self.model_forecast_lead_time.var_name].attrs = {
            'units': self.model_forecast_lead_time.units[0],
            'description': self.model_forecast_lead_time.description,
        }
        df[self.model_region_mean.var_name].attrs = {
            'units': self.model_region_mean.units[0],
            'description': self.model_region_mean.description,
        }
        df[self.model_region_std.var_name].attrs = {
            'units': self.model_region_std.units[0],
            'description': self.model_region_std.description,
        }
        df[self.model_region_var.var_name].attrs = {
            'units': self.model_region_var.units[0],
            'description': self.model_region_var.description,
        }

        # Global attributes

        df.attrs['model'] = self.model_name[0]
        df.attrs['model_variable'] = self.model_variable[0]
        df.attrs['model_variable_long_name'] = self.model_variable_long_name[0]
        df.attrs['model_region_name'] = self.model_region_name[0]

        return df

    def DataSet(self) -> xr.Dataset:
        '''
        Export the accumulated series to an Xarray Dataset.

        Dimensions
        ----------
        - time
            Shared dimension for all time-dependent series, based on `model_forecast_time`.

        Coordinates
        -----------
        - time : datetime64[ns], naive UTC (represents UTC)

        Data Variables
        --------------
        - model_forecast_init_time : naive UTC (represents UTC) (datetime64[ns])
        - model_forecast_lead_time : forecast lead time in seconds
        - model_mean / model_std / model_var : Regional model statistics

        Returns
        -------
        xarray.Dataset
            A new Dataset object. Variables and coordinates are copies, not views.

        Raises
        ------
        AssertionError
            If internal consistency checks fail.
        '''

        self.test_internal_consistency()

        # Coordinates

        coordinates = {}

        # Ensure time coordinate is datetime64[ns] and timezone-naive UTC
        coordinates['time'] = np.array(
            [dt.astimezone(timezone.utc).replace(tzinfo=None) for dt in self.model_forecast_time], dtype='datetime64[ns]'
        )

        # Variables

        data_vars = {}

        # Ensure model initialization time is datetime64[ns] and timezone-naive UTC
        init_time = np.array(
            [dt.astimezone(timezone.utc).replace(tzinfo=None) for dt in self.model_forecast_init_time], dtype='datetime64[ns]'
        )
        data_vars[self.model_forecast_init_time.var_name] = ('time', init_time)

        # Model lead time in hours
        lead_td = np.array(self.model_forecast_lead_time, dtype='timedelta64[ns]')
        lead_seconds = lead_td / np.timedelta64(1, 'h')
        data_vars[self.model_forecast_lead_time.var_name] = ('time', lead_seconds)

        data_vars[self.model_region_mean.var_name] = ('time', self.model_region_mean)
        data_vars[self.model_region_std.var_name] = ('time', self.model_region_std)
        data_vars[self.model_region_var.var_name] = ('time', self.model_region_var)

        ds = xr.Dataset(data_vars=data_vars, coords=coordinates)

        ds['time'].attrs['long_name'] = 'Forecast valid time (UTC)'
        ds['time'].attrs['units'] = 'seconds since 1970-01-01T00:00:00Z'
        ds['time'].attrs['calendar'] = 'proleptic_gregorian'

        ds[self.model_forecast_init_time.var_name].attrs['long_name'] = self.model_forecast_init_time.description
        ds[self.model_forecast_init_time.var_name].attrs['units'] = 'seconds since 1970-01-01T00:00:00Z'
        ds[self.model_forecast_init_time.var_name].attrs['calendar'] = 'proleptic_gregorian'

        # Units explicitly "seconds" for converted lead time (a duration, not a date)
        ds[self.model_forecast_lead_time.var_name].attrs['units'] = 'hours'
        ds[self.model_forecast_lead_time.var_name].attrs['long_name'] = self.model_forecast_lead_time.description

        ds[self.model_region_mean.var_name].attrs['units'] = self.model_region_mean.units[0]
        ds[self.model_region_std.var_name].attrs['units'] = self.model_region_std.units[0]
        ds[self.model_region_var.var_name].attrs['units'] = self.model_region_var.units[0]

        ds[self.model_region_mean.var_name].attrs['long_name'] = self.model_region_mean.description
        ds[self.model_region_std.var_name].attrs['long_name'] = self.model_region_std.description
        ds[self.model_region_var.var_name].attrs['long_name'] = self.model_region_var.description

        ds.attrs['model'] = self.model_name[0]
        ds.attrs['model_variable'] = self.model_variable[0]
        ds.attrs['model_variable_long_name'] = self.model_variable_long_name[0]
        ds.attrs['model_region_name'] = self.model_region_name[0]

        # Ensure variables and coordinates are copies, not views
        ds = ds.copy(deep=True)

        return ds

    def write_ds2netcdf(self, ds: xr.Dataset, dir_path: Path):
        '''
        Writes an xarray ds holding the time series data of this class into a netCDF file.

        Args:
            ds (xarray.Dataset): xarray holding the time series data of this class
            dir_path (Path): Directory where the netCDF file will be saved.
                             The directory will be created if it does not exist.
                             The file will be overwritten if it exists.
        '''

        # Construct file name

        tvals = ds['time'].values
        if np.issubdtype(tvals.dtype, np.datetime64):
            t64 = tvals.astype('datetime64[ns]')
        else:
            idx = pd.to_datetime(tvals, utc=True)
            t64 = idx.tz_localize(None).values.astype('datetime64[ns]')

        min_iso = np.datetime_as_string(t64.min(), unit='D')
        max_iso = np.datetime_as_string(t64.max(), unit='D')

        file_prefix = ds.attrs['model'] + '.' + min_iso + '-' + max_iso
        file_suffix = ds.attrs['model_variable'] + '.' + self.model_region_name[0]
        nc_file_name = Path(file_prefix + '.' + file_suffix + '.nc')

        # File path

        file_path = dir_path / nc_file_name

        # Remove attribute/encoding conflicts with any units set previously on time-like variables

        init_name = self.model_forecast_init_time.var_name

        for name in ('time', init_name):
            if name in ds.variables:
                # Remove conflicting 'units' or 'calendar' attributes if present
                for key in ('units', 'calendar'):
                    if key in ds[name].attrs:
                        del ds[name].attrs[key]

        # Build encoding

        encoding = {
            'time': {
                'dtype': 'float64',
                'units': 'seconds since 1970-01-01T00:00:00Z',
                'calendar': 'proleptic_gregorian',
            },
            init_name: {
                'dtype': 'float64',
                'units': 'seconds since 1970-01-01T00:00:00Z',
                'calendar': 'proleptic_gregorian',
            },
            **{
                var: {'dtype': 'float32'}
                for var in ds.data_vars
                if np.issubdtype(ds[var].dtype, np.floating) and var not in (init_name)
            },
        }

        # Identify any variables in the xarray that hold objects (such as datatime objects)
        object_vars = [name for name, var in ds.data_vars.items() if var.dtype == object]

        # Save to netCDF skipping any variables that are objects

        file_path.parent.mkdir(parents=True, exist_ok=True)

        if object_vars:
            ds.drop_vars(object_vars).to_netcdf(file_path, encoding=encoding)
        else:
            ds.to_netcdf(file_path, encoding=encoding)
        return


class ModelVsObs2DStatisticsTimeSeries:
    '''
    Compute and hold model-vs-observation evaluation statistics for 2-D fields.

    Upon initialization with the default constructor, the instance will hold
    a one-element time series (one-element list) of each statistical quantity
    obtaind by:

    1) Parsing model initialization time and forecast lead time to obtain the forecast valid time
    2) Matching that valid time to a unique timestamp in the observations,
    3) Interpolating the model 2-D field to observation locations,
    4) Computing model and observation summary statistics at reporting sites,
    5) Computing comparison metrics: bias, RMSE, Pearson correlation,
       and reporting coverage metrics.

    Time series are represented with list objects for performance, rather than with
    Pandas DataFrames, as the DataFrame.append() operation is slow.

    The time series can be expanded by calling the 'extend' instance method with
    the default constructor as an argument.
    '''

    def __init__(
        self,
        model_ds: xr.Dataset,
        model_variable: str,
        model_lat_name: str,
        model_lon_name: str,
        obs: object,
        obs_variable: str,
        obs_lat_name: str,
        obs_lon_name: str,
    ):
        '''
        Initialize the object with one timestep of model-vs-obs statistics.

        Behavior
        --------
        - Asserts matching physical units between model variable and observations.
        - Extracts model valid time from variable attributes and finds exactly one
          matching time in "obs.observations['UTC']" (UTC, timezone-aware).
        - Interpolates the model field to observation locations using
          "arotake.interpolation.model_2D_interpolate".
        - Identifies reporting vs non-reporting locations via NaN screening of the
          observation field at the matched time.
        - Computes:
            * Model stats at reporting locations: mean, std, var.
            * Observation stats at reporting locations: mean, std, var (NaN if none report).
            * Model-obs metrics at reporting locations: bias, RMSE, Pearson r
              (NaN if none report).
            * Coverage: counts of reporting/not-reporting locations and fraction not reporting.
        - Stores values and metadata in list-derived series. If any required parameter
          is "None", initializes an empty container.

        Parameters
        ----------
        model_ds : xarray.Dataset
            Dataset with 2-D model variable and 2-D lat/lon plus required timing attrs.
        model_variable : str
            Name of the model variable to evaluate.
        model_lat_name : str
            Name of the model latitude variable (2-D, y×x).
        model_lon_name : str
            Name of the model longitude variable (2-D, y×x).
        obs : object
            Container with attribute "observations" exposing the observation dataset.
        obs_variable : str
            Name of the observation variable to compare against.
        obs_lat_name : str
            Name of the observation latitude array variable.
        obs_lon_name : str
            Name of the observation longitude array variable.

        Raises
        ------
        AssertionError
            If model/observation units differ,
            if model forecast time units are not "hours",
            if the model valid time is absent or appears more than once in observations.
        '''

        # If any of the parameters is None, initialize an empty object

        if any(obj is None for obj in [model_ds, model_variable, model_lat_name, model_lon_name, obs, obs_variable]):
            # Make sure to initialize lists for time series values and units

            self.model_name = []

            self.model_region_name = []

            self.model_variable = TimeList()
            self.model_variable_long_name = TimeList()

            self.obs_name = []

            self.obs_variable = TimeList()
            self.obs_variable_long_name = TimeList()

            self.model_forecast_init_time = TimeList()
            self.model_forecast_lead_time = TimeList()
            self.model_forecast_time = TimeList()
            self.model_mean = TimeList()
            self.model_std = TimeList()
            self.model_var = TimeList()
            self.obs_mean = TimeList()
            self.obs_std = TimeList()
            self.obs_var = TimeList()
            self.model_bias = TimeList()
            self.model_rmse = TimeList()
            self.model_r_corr = TimeList()
            self.loc_rep_n = TimeList()
            self.loc_notrep_n = TimeList()
            self.loc_notrep_frac = TimeList()

            self.model_forecast_init_time.units = []
            self.model_forecast_lead_time.units = []
            self.model_forecast_time.units = []
            self.model_mean.units = []
            self.model_std.units = []
            self.model_var.units = []
            self.obs_mean.units = []
            self.obs_std.units = []
            self.obs_var.units = []
            self.model_bias.units = []
            self.model_rmse.units = []
            self.model_r_corr.units = []
            self.loc_rep_n.units = []
            self.loc_notrep_n.units = []
            self.loc_notrep_frac.units = []

            self.model_forecast_init_time.var_name = None
            self.model_forecast_lead_time.var_name = None
            self.model_forecast_time.var_name = None
            self.model_mean.var_name = None
            self.model_std.var_name = None
            self.model_var.var_name = None
            self.obs_mean.var_name = None
            self.obs_std.var_name = None
            self.obs_var.var_name = None
            self.model_bias.var_name = None
            self.model_rmse.var_name = None
            self.model_r_corr.var_name = None
            self.loc_rep_n.var_name = None
            self.loc_notrep_n.var_name = None
            self.loc_notrep_frac.var_name = None

            self.model_forecast_init_time.description = None
            self.model_forecast_lead_time.description = None
            self.model_forecast_time.description = None
            self.model_mean.description = None
            self.model_std.description = None
            self.model_var.description = None
            self.obs_mean.description = None
            self.obs_std.description = None
            self.obs_var.description = None
            self.model_bias.description = None
            self.model_rmse.description = None
            self.model_r_corr.description = None
            self.loc_rep_n.description = None
            self.loc_notrep_n.description = None
            self.loc_notrep_frac.description = None

            return

        # Model name

        if 'model' in model_ds.attrs:
            self.model_name = [model_ds.attrs['model']]
        else:
            self.model_name = ['model']

        # Isolate model variable

        model_data = model_ds[model_variable]
        model_data_units = model_data.attrs['units']

        # Data name

        if 'name' in obs.observations.attrs:
            self.obs_name = [obs.observations.attrs['name']]
        else:
            self.obs_name = ['observations']

        # Isolate observations variable

        obs_data = obs.observations.data_vars[obs_variable].values
        obs_data_units = obs.observations[obs_variable].attrs['units']

        # Check units
        assert model_data_units == obs_data_units, 'Mismatching units between model and observations'

        # Identify the model initialization time, forecast lead time, and forecast time

        self.model_forecast_init_time = datetime.strptime(model_data.attrs['initial_time'], '%m/%d/%Y (%H:%M)')
        self.model_forecast_init_time = self.model_forecast_init_time.replace(tzinfo=timezone.utc)

        assert model_data.attrs['forecast_time_units'] == 'hours', 'model forecast lead time units must be hours'

        self.model_forecast_lead_time_hours = int(model_data.attrs['forecast_time'])
        self.model_forecast_lead_time = timedelta(hours=self.model_forecast_lead_time_hours)

        self.model_forecast_time = self.model_forecast_init_time + self.model_forecast_lead_time

        # Identify the observation time matching the model forecast

        obs_time_index = np.where(obs.observations['UTC'].values[:] == self.model_forecast_time)[0]
        assert len(obs_time_index) > 0, 'model forecast time not found in observations data'
        assert len(obs_time_index) == 1, 'model forecast time present more than once in observations data'
        obs_time_index = obs_time_index[0]

        # obs_time = obs.observations['UTC'].values[obs_time_index]

        # Interpolate the model data to the observation locations

        model_data = interpolation.model_2D_interpolate(
            model_ds,
            model_variable,
            model_lat_name,
            model_lon_name,
            obs.observations.data_vars[obs_lat_name].values,
            obs.observations.data_vars[obs_lon_name].values,
        )

        # Identify observation locations that are/are not reporting

        loc_rep_index = np.where(~np.isnan(obs_data[obs_time_index, :]))[0]
        loc_notrep_index = np.where(np.isnan(obs_data[obs_time_index, :]))[0]

        self.loc_rep_n = len(loc_rep_index)
        self.loc_notrep_n = len(loc_notrep_index)

        loc_n = self.loc_rep_n + self.loc_notrep_n

        if loc_n > 0:
            self.loc_notrep_frac = self.loc_notrep_n / loc_n
        else:
            self.loc_notrep_frac = 0

        # Model statistics at observation locations
        self.model_mean = np.mean(model_data[loc_rep_index])
        self.model_std = np.std(model_data[loc_rep_index])
        self.model_var = np.var(model_data[loc_rep_index])

        if self.loc_rep_n > 0:
            # Observations statistics at observation locations
            self.obs_mean = np.mean(obs_data[obs_time_index, loc_rep_index])
            self.obs_std = np.std(obs_data[obs_time_index, loc_rep_index])
            self.obs_var = np.var(obs_data[obs_time_index, loc_rep_index])

            # Model vs observations bias
            self.model_bias = self.model_mean - self.obs_mean

            # Model vs observations root mean square error
            self.model_rmse = rmse(model_data[loc_rep_index], obs_data[obs_time_index, loc_rep_index])

            # Model vs observations correlation coefficient
            self.model_r_corr = np.corrcoef(model_data[loc_rep_index], obs_data[obs_time_index, loc_rep_index])[0, 1]

        else:
            self.obs_mean = np.nan
            self.obs_std = np.nan
            self.obs_var = np.nan
            self.model_bias = np.nan
            self.model_rmse = np.nan
            self.model_r_corr = np.nan

        # Make sure to initialize lists for time series values and units

        self.model_variable = TimeList([model_variable])

        if 'long_name' in obs.observations.attrs:
            self.model_variable_long_name = TimeList([model_ds[model_variable].attrs['long_name']])
        else:
            self.model_variable_long_name = TimeList([model_variable])

        self.obs_variable = TimeList([obs_variable])

        if 'long_name' in obs.observations[obs_variable].attrs:
            self.obs_variable_long_name = TimeList([obs.observations[obs_variable].attrs['long_name']])
        else:
            self.obs_variable_long_name = TimeList([obs_variable])

        self.model_region_name = ['-'.join(Model2DRegionalStatisticsTimeSeries._cached_gdf['name'].tolist())]

        self.model_forecast_init_time = TimeList([self.model_forecast_init_time])
        self.model_forecast_lead_time = TimeList([self.model_forecast_lead_time])
        self.model_forecast_time = TimeList([self.model_forecast_time])
        self.model_mean = TimeList([self.model_mean])
        self.model_std = TimeList([self.model_std])
        self.model_var = TimeList([self.model_var])
        self.obs_mean = TimeList([self.obs_mean])
        self.obs_std = TimeList([self.obs_std])
        self.obs_var = TimeList([self.obs_var])
        self.model_bias = TimeList([self.model_bias])
        self.model_rmse = TimeList([self.model_rmse])
        self.model_r_corr = TimeList([self.model_r_corr])
        self.loc_rep_n = TimeList([self.loc_rep_n])
        self.loc_notrep_n = TimeList([self.loc_notrep_n])
        self.loc_notrep_frac = TimeList([self.loc_notrep_frac])

        self.model_forecast_init_time.units = ['UTC']
        self.model_forecast_lead_time.units = ['h']
        self.model_forecast_time.units = ['UTC']
        self.model_mean.units = [model_data_units]
        self.model_std.units = [model_data_units]
        self.model_var.units = ['(' + model_data_units + ')^2']
        self.obs_mean.units = [obs_data_units]
        self.obs_std.units = [obs_data_units]
        self.obs_var.units = ['(' + obs_data_units + ')^2']
        self.model_bias.units = [model_data_units]
        self.model_rmse.units = [model_data_units]
        self.model_r_corr.units = ['']
        self.loc_rep_n.units = ['']
        self.loc_notrep_n.units = ['']
        self.loc_notrep_frac.units = ['']

        self.model_forecast_init_time.var_name = 'model_forecast_init_time'
        self.model_forecast_lead_time.var_name = 'model_forecast_lead_time'
        self.model_forecast_time.var_name = 'model_forecast_time'
        self.model_mean.var_name = 'model_mean_at_obs_locs'
        self.model_std.var_name = 'model_std_at_obs_locs'
        self.model_var.var_name = 'model_var_at_obs_locs'
        self.obs_mean.var_name = 'obs_mean_at_obs_locs'
        self.obs_std.var_name = 'obs_std_at_obs_locs'
        self.obs_var.var_name = 'obs_var_at_obs_locs'
        self.model_bias.var_name = 'model_vs_obs_bias_at_obs_locs'
        self.model_rmse.var_name = 'model_vs_obs_rmse_at_obs_locs'
        self.model_r_corr.var_name = 'model_vs_obs_r_corr_at_obs_locs'
        self.loc_rep_n.var_name = 'obs_loc_rep_n'
        self.loc_notrep_n.var_name = 'obs_loc_notrep_n'
        self.loc_notrep_frac.var_name = 'obs_loc_notrep_frac'

        self.model_forecast_init_time.description = 'forecast initialization time (UTC)'
        self.model_forecast_lead_time.description = 'forecast lead time'
        self.model_forecast_time.description = 'forecast time'
        self.model_mean.description = self.model_name[0] + ' mean at reporting ' + self.obs_name[0] + ' locations'
        self.model_std.description = self.model_name[0] + ' std. dev. at reporting ' + self.obs_name[0] + ' locations'
        self.model_var.description = self.model_name[0] + ' variance at reporting ' + self.obs_name[0] + ' locations'
        self.obs_mean.description = self.obs_name[0] + ' mean at reporting locations'
        self.obs_std.description = self.obs_name[0] + ' std. dev. at reporting locations'
        self.obs_var.description = self.obs_name[0] + ' variance at reporting locations'
        self.model_bias.description = self.model_name[0] + '-' + self.obs_name[0] + ' bias at reporting locations'
        self.model_rmse.description = self.model_name[0] + '-' + self.obs_name[0] + ' RMSE at reporting locations'
        self.model_r_corr.description = self.model_name[0] + '-' + self.obs_name[0] + ' correlation at reporting locations'
        self.loc_rep_n.description = 'number of reporting ' + self.obs_name[0] + ' locations'
        self.loc_notrep_n.description = 'number of not reporting ' + self.obs_name[0] + ' locations'
        self.loc_notrep_frac.description = 'fraction of not reporting ' + self.obs_name[0] + ' locations'

        return

    def test_internal_consistency(self):
        '''
        Validate internal consistency across all accumulated timesteps.

        Checks that model and observation variable identifiers, forecast times
        (by HH:MM:SS), lead times, and all recorded units are constant across
        the series. Raises AssertionError on inconsistency. Does nothing for
        an empty container.
        '''

        if self.model_forecast_init_time:
            assert len(set(self.obs_name)) <= 1, 'Time series covers different observations'

            assert len(set(self.obs_variable)) <= 1, 'Time series is internally inconsistent in the observation variables'

            assert (
                len(set(self.obs_variable_long_name)) <= 1
            ), 'Time series is internally inconsistent in the observation variable long names'

            assert len(set(self.model_name)) <= 1, 'Time series covers different models'

            assert len(set(self.model_region_name)) <= 1, 'Time series covers different regions'

            assert len(set(self.model_variable)) <= 1, 'Time series is internally inconsistent in the model variables'

            assert (
                len(set(self.model_variable_long_name)) <= 1
            ), 'Time series is internally inconsistent in the model variable long names'

            assert len(set([string.strftime("%H:%M:%S") for string in self.model_forecast_init_time[:]])) <= 1, (
                'Time series is internally inconsistent in the ' + self.model_forecast_init_time.description
            )

            assert len(set(self.model_forecast_lead_time)) <= 1, (
                'Time series is internally inconsistent in the ' + self.model_forecast_lead_time.description
            )

            assert len(set(self.model_forecast_init_time.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_forecast_init_time.description
            )
            assert len(set(self.model_forecast_time.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_forecast_time.description
            )
            assert len(set(self.model_forecast_lead_time.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_forecast_lead_time.description
            )

            assert len(set(self.model_mean.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_mean.description
            )
            assert len(set(self.model_std.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_std.description
            )
            assert len(set(self.model_var.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_var.description
            )

            assert len(set(self.obs_mean.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.obs_mean.description
            )
            assert len(set(self.obs_std.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.obs_std.description
            )
            assert len(set(self.obs_var.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.obs_var.description
            )

            assert len(set(self.model_bias.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_bias.description
            )
            assert len(set(self.model_rmse.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_rmse.description
            )
            assert len(set(self.model_r_corr.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.model_r_corr.description
            )

            assert len(set(self.loc_rep_n.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.loc_rep_n.description
            )
            assert len(set(self.loc_notrep_n.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.loc_notrep_n.description
            )
            assert len(set(self.loc_notrep_frac.units)) <= 1, (
                'Time series is internally inconsistent in the units for ' + self.loc_notrep_frac.description
            )

    def extend(self, other: 'ModelVsObs2DStatisticsTimeSeries'):
        '''
        Append another instance's time series to this instance's time series.

        Extends all value, unit, description, and name lists in place from
        "other" and then validates internal consistency.

        Parameters
        ----------
        other : ModelVsObs2DStatisticsTimeSeries
            Another container produced with the same configuration whose contents
            should be appended as the next timestep.
        '''

        # Append

        self.model_name.extend(other.model_name)
        self.model_region_name.extend(other.model_region_name)
        self.model_variable.extend(other.model_variable)
        self.model_variable_long_name.extend(other.model_variable_long_name)

        self.obs_name.extend(other.obs_name)
        self.obs_variable.extend(other.obs_variable)
        self.obs_variable_long_name.extend(other.obs_variable_long_name)

        self.model_forecast_init_time.extend(other.model_forecast_init_time)
        self.model_forecast_lead_time.extend(other.model_forecast_lead_time)
        self.model_forecast_time.extend(other.model_forecast_time)
        self.model_mean.extend(other.model_mean)
        self.model_std.extend(other.model_std)
        self.model_var.extend(other.model_var)
        self.obs_mean.extend(other.obs_mean)
        self.obs_std.extend(other.obs_std)
        self.obs_var.extend(other.obs_var)
        self.model_bias.extend(other.model_bias)
        self.model_rmse.extend(other.model_rmse)
        self.model_r_corr.extend(other.model_r_corr)
        self.loc_rep_n.extend(other.loc_rep_n)
        self.loc_notrep_n.extend(other.loc_notrep_n)
        self.loc_notrep_frac.extend(other.loc_notrep_frac)

        self.model_forecast_init_time.units.extend(other.model_forecast_init_time.units)
        self.model_forecast_lead_time.units.extend(other.model_forecast_lead_time.units)
        self.model_forecast_time.units.extend(other.model_forecast_time.units)
        self.model_mean.units.extend(other.model_mean.units)
        self.model_std.units.extend(other.model_std.units)
        self.model_var.units.extend(other.model_var.units)
        self.obs_mean.units.extend(other.obs_mean.units)
        self.obs_std.units.extend(other.obs_std.units)
        self.obs_var.units.extend(other.obs_var.units)
        self.model_bias.units.extend(other.model_bias.units)
        self.model_rmse.units.extend(other.model_rmse.units)
        self.model_r_corr.units.extend(other.model_r_corr.units)
        self.loc_rep_n.units.extend(other.loc_rep_n.units)
        self.loc_notrep_n.units.extend(other.loc_notrep_n.units)
        self.loc_notrep_frac.units.extend(other.loc_notrep_frac.units)

        self.model_forecast_init_time.var_name = other.model_forecast_init_time.var_name
        self.model_forecast_lead_time.var_name = other.model_forecast_lead_time.var_name
        self.model_forecast_time.var_name = other.model_forecast_time.var_name
        self.model_mean.var_name = other.model_mean.var_name
        self.model_std.var_name = other.model_std.var_name
        self.model_var.var_name = other.model_var.var_name
        self.obs_mean.var_name = other.obs_mean.var_name
        self.obs_std.var_name = other.obs_std.var_name
        self.obs_var.var_name = other.obs_var.var_name
        self.model_bias.var_name = other.model_bias.var_name
        self.model_rmse.var_name = other.model_rmse.var_name
        self.model_r_corr.var_name = other.model_r_corr.var_name
        self.loc_rep_n.var_name = other.loc_rep_n.var_name
        self.loc_notrep_n.var_name = other.loc_notrep_n.var_name
        self.loc_notrep_frac.var_name = other.loc_notrep_frac.var_name

        self.model_forecast_init_time.description = other.model_forecast_init_time.description
        self.model_forecast_lead_time.description = other.model_forecast_lead_time.description
        self.model_forecast_time.description = other.model_forecast_time.description
        self.model_mean.description = other.model_mean.description
        self.model_std.description = other.model_std.description
        self.model_var.description = other.model_var.description
        self.obs_mean.description = other.obs_mean.description
        self.obs_std.description = other.obs_std.description
        self.obs_var.description = other.obs_var.description
        self.model_bias.description = other.model_bias.description
        self.model_rmse.description = other.model_rmse.description
        self.model_r_corr.description = other.model_r_corr.description
        self.loc_rep_n.description = other.loc_rep_n.description
        self.loc_notrep_n.description = other.loc_notrep_n.description
        self.loc_notrep_frac.description = other.loc_notrep_frac.description

        # Test if the time series are internally consistent

        self.test_internal_consistency()

        return

    def DataFrame(self) -> pd.DataFrame:
        '''
        Export the accumulated series to a pandas DataFrame.

        Columns
        -------
        - model_forecast_time (datetime, timezone-aware UTC)
        - model_forecast_init_time (datetime, timezone-aware UTC)
        - model_forecast_lead_time (datetime.timedelta)
        - model_mean_at_obs_locs / model_std_at_obs_locs / model_var_at_obs_locs
        - obs_mean_at_obs_locs / obs_std_at_obs_locs / obs_var_at_obs_locs
        - model_vs_obs_bias_at_obs_locs / model_vs_obs_rmse_at_obs_locs /
          model_vs_obs_r_corr_at_obs_locs
        - obs_loc_rep_n / obs_loc_notrep_n / obs_loc_notrep_frac

        Implementation detail
        ---------------------
        Datetime and timedelta columns are created with dtype='object' to preserve
        native Python datetime/timedelta objects and avoid implicit conversion to
        pandas dtypes.

        Metadata
        --------
        - Per-column: stored in "Series.attrs" as "{'units', 'description'}".
        - Global: "DataFrame.attrs" includes
          "model_variable", "model_variable_long_name",
          "obs_variable", "obs_variable_long_name".

        Returns
        -------
        pandas.DataFrame
            A new DataFrame; arrays do not alias the internal lists.

        Raises
        ------
        AssertionError
            If internal consistency checks fail.
        '''

        # Test if the time series are internally consistent

        self.test_internal_consistency()

        # Create the data frame with one time series per column:

        df = pd.DataFrame(
            {
                self.model_forecast_time.var_name: pd.Series(
                    self.model_forecast_time, dtype='object'
                ),  # Preserve the datetime type which Pandas would otherwise convert to pandas.datetime64[ns]
                self.model_forecast_init_time.var_name: pd.Series(
                    self.model_forecast_init_time, dtype='object'
                ),  # Preserve the datetime type which Pandas would otherwise convert to pandas.datetime64[ns]
                self.model_forecast_lead_time.var_name: pd.Series(
                    self.model_forecast_lead_time, dtype='object'
                ),  # Preserve the datetime type which Pandas would otherwise convert to pandas.datetime64[ns]
                self.model_mean.var_name: self.model_mean,
                self.model_std.var_name: self.model_std,
                self.model_var.var_name: self.model_var,
                self.obs_mean.var_name: self.obs_mean,
                self.obs_std.var_name: self.obs_std,
                self.obs_var.var_name: self.obs_var,
                self.model_bias.var_name: self.model_bias,
                self.model_rmse.var_name: self.model_rmse,
                self.model_r_corr.var_name: self.model_r_corr,
                self.loc_rep_n.var_name: self.loc_rep_n,
                self.loc_notrep_n.var_name: self.loc_notrep_n,
                self.loc_notrep_frac.var_name: self.loc_notrep_frac,
            }
        )

        # Store units and description as per-column attributes

        df[self.model_forecast_time.var_name].attrs = {
            'units': self.model_forecast_time.units[0],
            'description': self.model_forecast_time.description,
        }
        df[self.model_forecast_init_time.var_name].attrs = {
            'units': self.model_forecast_init_time.units[0],
            'description': self.model_forecast_init_time.description,
        }
        df[self.model_forecast_lead_time.var_name].attrs = {
            'units': self.model_forecast_lead_time.units[0],
            'description': self.model_forecast_lead_time.description,
        }
        df[self.model_mean.var_name].attrs = {'units': self.model_mean.units[0], 'description': self.model_mean.description}
        df[self.model_std.var_name].attrs = {'units': self.model_std.units[0], 'description': self.model_std.description}
        df[self.model_var.var_name].attrs = {'units': self.model_var.units[0], 'description': self.model_var.description}
        df[self.obs_mean.var_name].attrs = {'units': self.obs_mean.units[0], 'description': self.obs_mean.description}
        df[self.obs_std.var_name].attrs = {'units': self.obs_std.units[0], 'description': self.obs_std.description}
        df[self.obs_var.var_name].attrs = {'units': self.obs_var.units[0], 'description': self.obs_var.description}
        df[self.model_bias.var_name].attrs = {'units': self.model_bias.units[0], 'description': self.model_bias.description}
        df[self.model_rmse.var_name].attrs = {'units': self.model_rmse.units[0], 'description': self.model_rmse.description}
        df[self.model_r_corr.var_name].attrs = {
            'units': self.model_r_corr.units[0],
            'description': self.model_r_corr.description,
        }
        df[self.loc_rep_n.var_name].attrs = {'units': self.loc_rep_n.units[0], 'description': self.loc_rep_n.description}
        df[self.loc_notrep_n.var_name].attrs = {
            'units': self.loc_notrep_n.units[0],
            'description': self.loc_notrep_n.description,
        }
        df[self.loc_notrep_frac.var_name].attrs = {
            'units': self.loc_notrep_frac.units[0],
            'description': self.loc_notrep_frac.description,
        }

        # Global attributes

        df.attrs['model'] = self.model_name[0]
        df.attrs['model_variable'] = self.model_variable[0]
        df.attrs['model_variable_long_name'] = self.model_variable_long_name[0]
        df.attrs['model_region_name'] = self.model_region_name[0]
        df.attrs['obs_variable'] = self.obs_variable[0]
        df.attrs['obs_variable_long_name'] = self.obs_variable_long_name[0]

        return df

    def DataSet(self) -> xr.Dataset:
        '''
        Export the accumulated series to an Xarray Dataset.

        Dimensions
        ----------
        - time
            Shared dimension for all time-dependent series, based on `model_forecast_time`.

        Coordinates
        -----------
        - time : datetime64[ns], naive UTC (represents UTC)

        Data Variables
        --------------
        - model_forecast_init_time : naive UTC (represents UTC) (datetime64[ns])
        - model_forecast_lead_time : forecast lead time in seconds
        - model_mean / model_std / model_var : Regional model statistics
        - obs_mean / obs_std / obs_var : Regional observed statistics
        - model_bias / model_rmse / model_r_corr : Model-obs comparison metrics
        - loc_rep_n / loc_notrep_n / loc_notrep_frac : Location representativity diagnostics.

        Returns
        -------
        xarray.Dataset
            A new Dataset object. Variables and coordinates are copies, not views.

        Raises
        ------
        AssertionError
            If internal consistency checks fail.
        '''

        self.test_internal_consistency()

        # Coordinates

        coordinates = {}

        # Ensure time coordinate is datetime64[ns] and timezone-naive UTC
        coordinates['time'] = np.array(
            [dt.astimezone(timezone.utc).replace(tzinfo=None) for dt in self.model_forecast_time], dtype='datetime64[ns]'
        )

        # Variables

        data_vars = {}

        # Ensure model initialization time is datetime64[ns] and timezone-naive UTC
        init_time = np.array(
            [dt.astimezone(timezone.utc).replace(tzinfo=None) for dt in self.model_forecast_init_time], dtype='datetime64[ns]'
        )
        data_vars[self.model_forecast_init_time.var_name] = ('time', init_time)

        # Model lead time in hours
        lead_td = np.array(self.model_forecast_lead_time, dtype='timedelta64[ns]')
        lead_seconds = lead_td / np.timedelta64(1, 'h')
        data_vars[self.model_forecast_lead_time.var_name] = ('time', lead_seconds)

        data_vars[self.model_mean.var_name] = ('time', self.model_mean)
        data_vars[self.model_std.var_name] = ('time', self.model_std)
        data_vars[self.model_var.var_name] = ('time', self.model_var)
        data_vars[self.obs_mean.var_name] = ('time', self.obs_mean)
        data_vars[self.obs_std.var_name] = ('time', self.obs_std)
        data_vars[self.obs_var.var_name] = ('time', self.obs_var)
        data_vars[self.model_bias.var_name] = ('time', self.model_bias)
        data_vars[self.model_rmse.var_name] = ('time', self.model_rmse)
        data_vars[self.model_r_corr.var_name] = ('time', self.model_r_corr)
        data_vars[self.loc_rep_n.var_name] = ('time', self.loc_rep_n)
        data_vars[self.loc_notrep_n.var_name] = ('time', self.loc_notrep_n)
        data_vars[self.loc_notrep_frac.var_name] = ('time', self.loc_notrep_frac)

        ds = xr.Dataset(data_vars=data_vars, coords=coordinates)

        ds['time'].attrs['long_name'] = 'Forecast valid time (UTC)'
        ds['time'].attrs['units'] = 'seconds since 1970-01-01T00:00:00Z'
        ds['time'].attrs['calendar'] = 'proleptic_gregorian'

        ds[self.model_forecast_init_time.var_name].attrs['long_name'] = self.model_forecast_init_time.description
        ds[self.model_forecast_init_time.var_name].attrs['units'] = 'seconds since 1970-01-01T00:00:00Z'
        ds[self.model_forecast_init_time.var_name].attrs['calendar'] = 'proleptic_gregorian'

        # Units explicitly "seconds" for converted lead time (a duration, not a date)
        ds[self.model_forecast_lead_time.var_name].attrs['units'] = 'hours'
        ds[self.model_forecast_lead_time.var_name].attrs['long_name'] = self.model_forecast_lead_time.description

        ds[self.model_mean.var_name].attrs['units'] = self.model_mean.units[0]
        ds[self.model_std.var_name].attrs['units'] = self.model_std.units[0]
        ds[self.model_var.var_name].attrs['units'] = self.model_var.units[0]
        ds[self.obs_mean.var_name].attrs['units'] = self.obs_mean.units[0]
        ds[self.obs_std.var_name].attrs['units'] = self.obs_std.units[0]
        ds[self.obs_var.var_name].attrs['units'] = self.obs_var.units[0]
        ds[self.model_bias.var_name].attrs['units'] = self.model_bias.units[0]
        ds[self.model_rmse.var_name].attrs['units'] = self.model_rmse.units[0]
        ds[self.model_r_corr.var_name].attrs['units'] = self.model_r_corr.units[0]
        ds[self.loc_rep_n.var_name].attrs['units'] = self.loc_rep_n.units[0]
        ds[self.loc_notrep_n.var_name].attrs['units'] = self.loc_notrep_n.units[0]
        ds[self.loc_notrep_frac.var_name].attrs['units'] = self.loc_notrep_frac.units[0]

        ds[self.model_mean.var_name].attrs['long_name'] = self.model_mean.description
        ds[self.model_std.var_name].attrs['long_name'] = self.model_std.description
        ds[self.model_var.var_name].attrs['long_name'] = self.model_var.description
        ds[self.obs_mean.var_name].attrs['long_name'] = self.obs_mean.description
        ds[self.obs_std.var_name].attrs['long_name'] = self.obs_std.description
        ds[self.obs_var.var_name].attrs['long_name'] = self.obs_var.description
        ds[self.model_bias.var_name].attrs['long_name'] = self.model_bias.description
        ds[self.model_rmse.var_name].attrs['long_name'] = self.model_rmse.description
        ds[self.model_r_corr.var_name].attrs['long_name'] = self.model_r_corr.description
        ds[self.loc_rep_n.var_name].attrs['long_name'] = self.loc_rep_n.description
        ds[self.loc_notrep_n.var_name].attrs['long_name'] = self.loc_notrep_n.description
        ds[self.loc_notrep_frac.var_name].attrs['long_name'] = self.loc_notrep_frac.description

        ds.attrs['model'] = self.model_name[0]
        ds.attrs['model_variable'] = self.model_variable[0]
        ds.attrs['model_variable_long_name'] = self.model_variable_long_name[0]
        ds.attrs['model_region_name'] = self.model_region_name[0]

        ds.attrs['observations'] = self.obs_name[0]
        ds.attrs['obs_variable'] = self.obs_variable[0]
        ds.attrs['obs_variable_long_name'] = self.obs_variable_long_name[0]

        # Ensure variables and coordinates are copies, not views
        ds = ds.copy(deep=True)

        return ds

    def write_ds2netcdf(self, ds: xr.Dataset, dir_path: Path):
        '''
        Writes an xarray ds holding the time series data of this class into a netCDF file.

        Args:
            ds (xarray.Dataset): xarray holding the time series data of this class
            dir_path (Path): Directory where the netCDF file will be saved.
                             The directory will be created if it does not exist.
                             The file will be overwritten if it exists.
        '''

        # Construct file name

        tvals = ds['time'].values
        if np.issubdtype(tvals.dtype, np.datetime64):
            t64 = tvals.astype('datetime64[ns]')
        else:
            idx = pd.to_datetime(tvals, utc=True)
            t64 = idx.tz_localize(None).values.astype('datetime64[ns]')

        min_iso = np.datetime_as_string(t64.min(), unit='D')
        max_iso = np.datetime_as_string(t64.max(), unit='D')

        file_prefix = ds.attrs['model'] + '_vs_' + ds.attrs['observations'] + '.' + min_iso + '-' + max_iso
        file_suffix = ds.attrs['model_variable'] + '.' + self.model_region_name[0]
        nc_file_name = Path(file_prefix + '.' + file_suffix + '.nc')

        # File path

        file_path = dir_path / nc_file_name

        # Remove attribute/encoding conflicts with any units set previously on time-like variables

        init_name = self.model_forecast_init_time.var_name

        for name in ('time', init_name):
            if name in ds.variables:
                # Remove conflicting 'units' or 'calendar' attributes if present
                for key in ('units', 'calendar'):
                    if key in ds[name].attrs:
                        del ds[name].attrs[key]

        # Build encoding

        encoding = {
            'time': {
                'dtype': 'float64',
                'units': 'seconds since 1970-01-01T00:00:00Z',
                'calendar': 'proleptic_gregorian',
            },
            init_name: {
                'dtype': 'float64',
                'units': 'seconds since 1970-01-01T00:00:00Z',
                'calendar': 'proleptic_gregorian',
            },
            **{
                var: {'dtype': 'float32'}
                for var in ds.data_vars
                if np.issubdtype(ds[var].dtype, np.floating) and var not in (init_name)
            },
        }

        # Identify any variables in the xarray that hold objects (such as datatime objects)
        object_vars = [name for name, var in ds.data_vars.items() if var.dtype == object]

        # Save to netCDF skipping any variables that are objects

        file_path.parent.mkdir(parents=True, exist_ok=True)

        if object_vars:
            ds.drop_vars(object_vars).to_netcdf(file_path, encoding=encoding)
        else:
            ds.to_netcdf(file_path, encoding=encoding)
        return


def rmse(reference: np.ndarray, prediction: np.ndarray) -> float:
    '''
    Returns the root mean square error.
    '''

    return np.sqrt(np.mean((prediction - reference) ** 2))
