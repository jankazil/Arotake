'''
Interpolation utilities for geospatial model fields.

This module provides:
- model_2D_interpolate: Interpolates a 2-D variable from an xarray Dataset
  onto arbitrary latitude/longitude points using bilinear regridding via
  xESMF (xesmf.Regridder). To improve performance, the function caches the
  regridder and reuses it when inputs have not changed.

Design notes
------------
- Input lats/lons must be 1-D, equal length, non-NaN, and within [-180, 180].
- The model Dataset must provide a 2-D variable plus corresponding 2-D latitude
  and longitude arrays (WGS84, degrees North/East).
- Regridders are rebuilt only if the model grid or the interpolation locations
  differ from those cached in the last call.
'''

import numpy as np
import xarray as xr
import xesmf as xe  # Universal Regridder for Geospatial Data


def model_2D_interpolate(
    model_ds: xr.Dataset, model_variable: str, model_lat_name: str, model_lon_name: str, lats: np.ndarray, lons: np.ndarray
) -> np.ndarray:
    '''
    Interpolate a 2-D model variable to arbitrary latitude/longitude points
    using bilinear regridding with xESMF.

    Caching
    -------
    To reduce overhead from repeated calls, the function caches an
    "xesmf.Regridder" instance plus the model grid and target coordinates.
    A new regridder is constructed only if:
      - No regridder has been cached yet,
      - The model grid latitude/longitude arrays differ from the cached ones,
      - The requested interpolation coordinates (lats/lons) differ.

    Parameters
    ----------
    model_ds : xarray.Dataset
        Dataset containing:
        - 2-D data variable "model_variable" (y, x),
        - 2-D latitude variable "model_lat_name" (y, x),
        - 2-D longitude variable "model_lon_name" (y, x),
        with longitudes constrained to [-180, 180].
    model_variable : str
        Name of the variable in "model_ds" to interpolate.
    model_lat_name : str
        Name of the dataset variable holding latitude coordinates.
    model_lon_name : str
        Name of the dataset variable holding longitude coordinates.
    lats : numpy.ndarray
        One-dimensional array of latitude values (degrees North, WGS84) at which
        to interpolate.
    lons : numpy.ndarray
        One-dimensional array of longitude values (degrees East, WGS84, in [-180, 180])
        at which to interpolate. Must be the same length as "lats".

    Returns
    -------
    numpy.ndarray
        Interpolated values of the model variable at the requested coordinates.
        If the model variable has a time dimension, interpolation is performed
        independently for each time step, returning a 2-D array.

    Raises
    ------
    KeyError
        If the requested variable or latitude/longitude variables are not found
        in the dataset.
    ValueError
        If longitudes fall outside [-180, 180], if "lats" and "lons" differ
        in length, or if either contains NaNs.
    '''

    # Test input

    if model_variable not in model_ds:
        raise KeyError(f"Variable '{model_variable}' not found in dataset.")

    if model_lat_name not in model_ds or model_lon_name not in model_ds:
        raise KeyError("Latitude/longitude variable names not found in dataset.")

    if np.min(model_ds[model_lon_name].values) < -180 or np.max(model_ds[model_lon_name].values) > 180:
        raise ValueError("Longitude of interpolation model outside of [-180,180].")

    if lats.shape[0] != lons.shape[0]:
        raise ValueError("'lats' and 'lons' must have equal length.")

    if np.isnan(lats).any() or np.isnan(lons).any():
        raise ValueError("'lats' and 'lons' must not contain NaN values.")

    if np.min(lons) < -180 or np.max(lons) > 180:
        raise ValueError("Longitude of interpolation locations outside of [-180,180].")

    # Xarray Dataset holding the model grid with appropriately named latitude and longitude dimensions

    lat_name = 'lat'
    lon_name = 'lon'

    model_grid_ds = model_ds.rename({model_lat_name: lat_name, model_lon_name: lon_name}).set_coords([lat_name, lon_name])[
        [lat_name, lon_name]
    ]

    # Xarray Dataset holding the interpolation coordinates

    interpolation_locs_ds = xr.Dataset(
        coords={
            lat_name: (('locations',), lats),
            lon_name: (('locations',), lons),
        }
    )

    # Create regridder if necessary, otherwise use the old regridder

    cached_model_regridder = getattr(model_2D_interpolate, '_cached_model_regridder', None)
    cached_model_grid = getattr(model_2D_interpolate, '_cached_model_grid', None)
    cached_lats = getattr(model_2D_interpolate, '_cached_lats', None)
    cached_lons = getattr(model_2D_interpolate, '_cached_lons', None)

    if (
        cached_model_regridder is None
        or not np.array_equal(cached_model_grid[lat_name], model_grid_ds[lat_name].values)
        or not np.array_equal(cached_model_grid[lon_name], model_grid_ds[lon_name].values)
        or not np.array_equal(cached_lats, lats)
        or not np.array_equal(cached_lons, lons)
    ):
        cached_model_regridder = xe.Regridder(
            model_grid_ds, interpolation_locs_ds, method='bilinear', locstream_out=True, periodic=False
        )
        model_2D_interpolate._cached_model_regridder = cached_model_regridder
        model_2D_interpolate._cached_model_grid = {
            lat_name: model_grid_ds[lat_name].values,
            lon_name: model_grid_ds[lon_name].values,
        }
        model_2D_interpolate._cached_lats = np.asarray(lats)
        model_2D_interpolate._cached_lons = np.asarray(lons)

    # Interpolate the model data to the interpolation locations

    model_data_interpolated = cached_model_regridder(model_ds[model_variable]).values

    return model_data_interpolated
