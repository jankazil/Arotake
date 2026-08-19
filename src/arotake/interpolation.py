'''
Interpolation utilities for geospatial model fields.

This module provides:
- model_2D_interpolate: Interpolates a variable from an xarray Dataset onto
  arbitrary latitude/longitude points using bilinear regridding via xESMF
  (xesmf.Regridder). The variable may have additional leading dimensions
  before its two horizontal dimensions. To improve performance, the function
  caches the regridder and reuses it when inputs have not changed.

Design notes
------------
- Input lats/lons must be 1-D, equal length, non-NaN, and within [-180, 180].
- The model Dataset must provide a variable with dimensions (..., y, x) plus
  corresponding 2-D latitude and longitude arrays (y, x; WGS84, degrees
  North/East). Interpolation is applied independently over any leading
  dimensions.
- Regridders are rebuilt only if the model grid or the interpolation locations
  differ from those cached in the last call.
'''

import numpy as np
import xarray as xr
import xesmf as xe  # Universal Regridder for Geospatial Data


def model_2D_interpolate(
    model_ds: xr.Dataset, model_variable: str, model_lat_name: str, model_lon_name: str, lats: np.ndarray, lons: np.ndarray
) -> xr.DataArray:
    '''
    Interpolate a model variable to arbitrary latitude/longitude points using
    bilinear regridding with xESMF. The variable may have arbitrary leading
    dimensions before its two horizontal dimensions; interpolation is applied
    independently to every combination of the leading-dimension coordinates.

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
        - data variable "model_variable" (..., y, x),
        - 2-D latitude variable "model_lat_name" (y, x),
        - 2-D longitude variable "model_lon_name" (y, x),
        with longitudes constrained to [-180, 180].
        The horizontal dimensions must be the final two dimensions of the data
        variable.
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
    xarray.DataArray
        Interpolated values with dimensions (..., locations). Any leading
        dimensions and their coordinates are retained; the two horizontal input
        dimensions (y, x) are replaced by the output dimension "locations".

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

    # Interpolate the two rightmost, horizontal dimensions to the requested
    # locations. xESMF applies the same interpolation independently over any
    # additional leading dimensions and retains their coordinates.

    model_data_interpolated = cached_model_regridder(model_ds[model_variable])

    return model_data_interpolated
