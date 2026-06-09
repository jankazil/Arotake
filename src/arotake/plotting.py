'''

Plotting utilities

This module provides the following functions:

- plot_locations_conus: Creates a map of observation locations over the contiguous U.S.
  in Lambert Conformal projection.

- plot_evaluation_time_series: Loads one or more Arotake model-vs-observation NetCDF
  dataset files, generates time series plots, and saves them as PNGs.

- plot_evaluation_pdf: Loads one or more Arotake model-vs-observation NetCDF dataset files,
  computes kernel-density estimates (Gaussian KDE) of the probability density, generates
  probability density plots, and saves them as PNGs.

- plot_df_timeseries:
  Create one PNG per statistic by overlaying the same non-time column from multiple
  pandas.DataFrame objects that share schema and metadata. The first column in each
  DataFrame is interpreted as time. Figure titles and axis labels are derived from
  DataFrame/Series attributes. Output filenames encode the global date range, statistic
  name, and per-DataFrame suffixes based on forecast valid time and lead hour.

- line_plot_1d:
  Helper for 1-D multi-series line plots used by higher-level functions. Handles
  labels, limits, legend, optional saving, and always closes the figure.

- pandas_series_plottable:
  Heuristic check for whether a pandas.Series can be plotted directly by Matplotlib
  without pre-conversion. Supports numeric, datetime64, and timedelta64 dtypes; object
  dtype is accepted only when the first non-null value is numeric.

'''

from datetime import datetime, timedelta
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pandas.api.types as ptypes
import xarray as xr
from cartopy.mpl.ticker import LatitudeFormatter, LatitudeLocator, LongitudeFormatter, LongitudeLocator
from scipy.stats import gaussian_kde

matplotlib.use("Agg")  # Important to avoid runaway memory use upon creating plots repeatedly.


def dataframe_column_attr(df: pd.DataFrame, column: str, attr: str) -> object:
    '''
    Return metadata for a DataFrame column.

    Prefer DataFrame-level column metadata stored in
    ``df.attrs['column_attrs'][column][attr]``. Fall back to the older
    ``df[column].attrs[attr]`` storage used by earlier code.
    '''

    if 'column_attrs' in df.attrs:
        try:
            return df.attrs['column_attrs'][column][attr]
        except KeyError:
            pass

    try:
        return df[column].attrs[attr]
    except KeyError as exc:
        raise KeyError(
            'Missing DataFrame column metadata: '
            f"df.attrs['column_attrs'][{column!r}][{attr!r}], and "
            f'missing Series.attrs[{attr!r}] for DataFrame column {column!r}'
        ) from exc


def plot_df_timeseries(plot_dir: Path, dfs: list[pd.DataFrame], legend_strings: list[str] = None) -> list[Path]:
    '''
    Plot one figure per statistic (per non-time column) by overlaying that statistic across
    multiple DataFrames and save each figure as a PNG under "plot_dir".

    Behavior
    --------
    - Interprets the first column in each DataFrame as the time axis; all subsequent columns
      are treated as statistics to plot.
    - For every statistic present in the FIRST DataFrame (columns[1:]), creates a single
      figure and draws one line per DataFrame.
    - Skips a statistic if "pandas_series_plottable" returns False for that column in the
      FIRST DataFrame.
    - X-axis limits span the global min and max time across ALL DataFrames.
    - Figure title is taken from "df.attrs['model_variable_long_name']" (which must match
      across all DataFrames).
    - Y-axis label uses the plottable column’s "attrs['description']" and "attrs['units']".
    - Legend entries are built as:
        * Without "legend_strings":
          "<HH:MM:SS><space>time_units forecast valid time, <lead_hours> h lead time"
        * With "legend_strings[j]":
          "<legend_strings[j]>, <HH:MM:SS><space>time_units forecast valid time, <lead_hours> h lead time"
      where "lead_hours" is derived from 'model_forecast_lead_time' in hours.
    - Output filenames follow:
        "<start_date>-<end_date>.<column_label>.<suffix>.png"
      where:
        * <start_date> and <end_date> are ISO dates of the global time range,
        * <column_label> is the statistic name.

    Requirements and Assumptions
    ----------------------------
    - All DataFrames:
        - Have the same number of columns and column order, where corresponding columns
          represent the same statistics.
        - Share the same "attrs['model_variable_long_name']" string.
        - Include a 'model_forecast_lead_time' column with timedelta-like values.
    - The first column in each DataFrame contains datetime values, and its Series has
      "attrs['units']" set to a human-readable time unit string for labeling.
    - Each statistic Series used for labeling provides "attrs['description']" and "attrs['units']".
    - If provided, "legend_strings" length equals "len(dfs)".

    Parameters
    ----------
    plot_dir : pathlib.Path
        Directory where PNG files will be written. Created if it does not exist.
    dfs : list[pandas.DataFrame]
        List of DataFrames to overlay. The first column must be datetime; columns[1:] are
        candidate statistics. Non-plottable statistics are skipped.
    legend_strings : list[str], optional
        Optional labels used both in the legend and as components of the filename suffix.
        Must match "len(dfs)" if provided.

    Returns
    -------
    list[pathlib.Path]
        Full paths of the saved PNG files, one per plottable statistic.

    Raises
    ------
    AssertionError
        If "legend_strings" is provided with a length different from "len(dfs)".
        If the DataFrames do not share the same "attrs['model_variable_long_name']".

    Notes
    -----
    - Figures are created with size (14, 7) and DPI 600.
    - The function does not modify the input DataFrames.
    '''

    # Check input

    if legend_strings is not None:
        assert len(dfs) == len(legend_strings), 'Length of the iteratable arguments dfs and legend_strings differs.'

    model_variable_long_names = set([df.attrs['model_variable_long_name'] for df in dfs])

    assert len(model_variable_long_names) == 1, 'The Pandas DataFrames provided hold different model variables.'

    model_variable_long_name = list(model_variable_long_names)[0]

    # Plot title

    plot_title = model_variable_long_name

    # Plot time range and other time parameters

    min_times = []
    max_times = []
    forecast_valid_time_strings = []
    forecast_lead_hour_strings = []

    for df in dfs:
        time = df[df.columns[0]]
        min_times.append(min(time))
        max_times.append(max(time))

        forecast_valid_time_string = (
            time.iloc[0].strftime("%H:%M:%S") + ' ' + str(dataframe_column_attr(df, df.columns[0], 'units'))
        )

        forecast_valid_time_strings.append(forecast_valid_time_string)

        model_forecast_lead_hour = str(round(df['model_forecast_lead_time'].iloc[0].total_seconds() / 3600))

        forecast_lead_hour_strings.append(model_forecast_lead_hour)

    min_time = min(min_times)
    max_time = max(max_times)

    # Plot file prefix

    file_prefix = min_time.date().isoformat() + '-' + max_time.date().isoformat()

    # Plot file suffix

    file_suffixes = []

    if legend_strings is None:
        for jj, _df in enumerate(dfs):
            file_suffixes.append(forecast_valid_time_strings[jj].replace(' ', '') + '_f' + forecast_lead_hour_strings[jj] + 'h')
    else:
        for jj, _df in enumerate(dfs):
            file_suffixes.append(
                legend_strings[jj].replace(' ', '_')
                + '.'
                + forecast_valid_time_strings[jj].replace(' ', '')
                + 'f'
                + forecast_lead_hour_strings[jj]
                + 'h'
            )

    file_suffix = '_'.join(file_suffixes)

    # Create plots

    plot_file_paths = []

    for column_label in dfs[0].columns[1:]:
        time_series = dfs[0][column_label]

        # Skip columns whose elements cannot be plotted

        if not pandas_series_plottable(time_series):
            continue

        fig, ax = plt.subplots(figsize=(14, 7), dpi=600)

        for jj, df in enumerate(dfs):
            time = df[df.columns[0]]
            time_series = df[column_label]

            label = forecast_valid_time_strings[jj] + ' forecast valid time, ' + forecast_lead_hour_strings[jj] + ' h lead time'

            if legend_strings is not None:
                label = legend_strings[jj] + ', ' + label

            ax.plot(time, time_series, linewidth=1, label=label)

            ax.set_title(plot_title)

        ax.set_xlim(min_time, max_time)
        # ax.set_ylim(min(time_series),max(time_series))

        ax.set_xlabel('Forecast valid date')
        ax.set_ylabel(
            str(dataframe_column_attr(df, column_label, 'description'))
            + ' ('
            + str(dataframe_column_attr(df, column_label, 'units'))
            + ')'
        )

        ax.legend(title=None)

        # Save and close figures

        plot_file_name = file_prefix + '.' + column_label + '.' + file_suffix + '.png'

        plot_file_path = plot_dir / plot_file_name

        plot_file_path.parent.mkdir(parents=True, exist_ok=True)

        fig.savefig(plot_file_path, bbox_inches='tight')

        plot_file_paths.append(plot_file_path)

        plt.close(fig)

    return plot_file_paths


def pandas_series_plottable(col: pd.Series) -> bool:
    '''
    Return True if "col" can be plotted directly by Matplotlib without pre-conversion,
    otherwise return False.

    Rules
    -----
    - Returns True for numeric dtypes.
    - Returns True for any datetime64 dtype.
    - Returns True for timedelta64 dtype.
    - For object dtype:
        * If the first non-null element is an int or a float, return True.
        * If the first non-null element is a "datetime.datetime" or "datetime.timedelta",
          return False because explicit conversion (for example, "pd.to_datetime",
          "pd.to_timedelta", or ".total_seconds()") is expected before plotting.
        * Otherwise return False.
    - Returns False for empty Series (all values NA) and for unsupported dtypes.

    Parameters
    ----------
    col : pandas.Series
        Series to test for direct plot compatibility.

    Returns
    -------
    bool
        True if Matplotlib can plot the Series without additional conversion, else False.
    '''

    if ptypes.is_numeric_dtype(col):
        return True

    if ptypes.is_datetime64_any_dtype(col):
        return True

    if ptypes.is_timedelta64_dtype(col):
        return True

    if ptypes.is_object_dtype(col):
        # Inspect first non-null element
        if col.dropna().empty:
            return False
        first_non_null = col.dropna().iloc[0]
        if isinstance(first_non_null, (int | float)):
            return True
        if isinstance(first_non_null, datetime):
            return False  # requires conversion with pd.to_datetime
        if isinstance(first_non_null, timedelta):
            return False  # requires conversion with pd.to_timedelta or .total_seconds()
        return False

    return False


def plot_locations_conus(
    data_obss: list[object] | object,
    start_dates: list[datetime] | datetime,
    end_dates: list[datetime] | datetime,
    obs_lat_names: list[str] | str,
    obs_lon_names: list[str] | str,
    plot_dir: Path,
    colors: list[str] | str = None,
    markersizes: list[int] | int = None,
    legend_font_size: int = 10,
    verbose: bool = False,
) -> Path:
    '''
    Generates and saves a map showing the locations of observation stations within the
    contiguous United States (CONUS) using a Lambert Conformal projection consistent
    with the HRRR model domain.

    This function accepts one or more observation dataset objects, each containing a
    `.time_series` xarray Dataset with latitude and longitude coordinates for station
    locations. Each dataset is plotted as a separate layer with distinct marker sizes,
    allowing visual differentiation between datasets. The resulting figure includes
    coastlines, state and national borders, gridlines with labeled coordinates, and a
    legend describing each dataset and its associated observation period.

    Args:
        data_obss (list[object] | object): One or more objects containing a `.time_series`
            xarray Dataset with station observation data. Each dataset must include
            attributes `'long_name'`, `'region'`, and `'name'`, as well as coordinate
            variables for station latitude and longitude.
        start_dates (list[datetime] | datetime): One or more start dates corresponding
            to the observation periods for each dataset. Used in legend labeling.
        end_dates (list[datetime] | datetime): One or more end dates corresponding
            to the observation periods for each dataset. Used in legend labeling.
        obs_lat_names (list[str] | str): One or more names of the latitude variables
            within the `.time_series` Datasets.
        obs_lon_names (list[str] | str): One or more names of the longitude variables
            within the `.time_series` Datasets.
        plot_dir (Path): Directory where the generated map image will be saved. The
            directory will be created if it does not exist.
        colors (list[str] | str), optional: One or more colors for location markers.
        markersizes (list[int] | int), optional: One or more integers specifying the
            size of location markers.
        legend_font_size : int, optional
            Font size for figure legend.
        verbose : bool, optional
            If True, print progress information.

    Returns:
        Path: Path to the saved PNG image file showing the station map.

    Raises:
        AssertionError: If the provided input lists do not all have the same length.

    Notes:
        - The projection used is Lambert Conformal Conic, centered at 97.5°W and 38.5°N,
          matching the HRRR domain.
        - The map extent spans approximately 25°N–50°N and 123°W–71°W, covering the
          continental United States.
        - Each dataset is plotted with a distinct marker size (larger for earlier entries).
        - Gridlines include labeled parallels and meridians with degree symbols.
        - The output image is saved at 600 DPI with the filename constructed from the
          dataset regions and names, joined by periods, followed by `".stations_map.png"`.
        - A legend indicates each dataset’s long name, region, and date range.
    '''

    # Normalize all inputs to lists

    if not isinstance(data_obss, list):
        data_obss = [data_obss]
    if not isinstance(start_dates, list):
        start_dates = [start_dates]
    if not isinstance(end_dates, list):
        end_dates = [end_dates]
    if not isinstance(obs_lat_names, list):
        obs_lat_names = [obs_lat_names]
    if not isinstance(obs_lon_names, list):
        obs_lon_names = [obs_lon_names]
    if colors is not None:
        if not isinstance(colors, list):
            colors = [colors]
    else:
        colors = [None for data_obs in data_obss]
    if markersizes is not None:
        if not isinstance(markersizes, list):
            markersizes = [markersizes]
    else:
        markersizes = list(range(len(data_obss), 0, -1))

    # Make sure that all lists have the same length:

    assert (
        len(data_obss) == len(start_dates) == len(end_dates) == len(obs_lat_names) == len(obs_lon_names) == len(colors)
    ), 'Arguments are lists of different length. Aborting.'

    # Define the target map projection (Lambert Conformal as in HRRR)
    hrrr_proj = ccrs.LambertConformal(
        central_longitude=-97.5,
        central_latitude=38.5,
        standard_parallels=(38.5, 38.5),
    )

    # Crate plot - here, passing a Cartopy projection will make 'ax' being a cartopy.mpl.geoaxes.GeoAxes object, and Cartopy methods can be used with 'ax'

    fig, ax = plt.subplots(figsize=(16, 8), subplot_kw={"projection": hrrr_proj})

    # Set the map extent (using the PlateCarree coordinate reference system which just means lon/lat)
    ax.set_extent([-123, -71, 25, 50], crs=ccrs.PlateCarree())

    # Coastlines, national borders, state lines

    ax.coastlines(resolution="50m", linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.25)
    ax.add_feature(cfeature.STATES, linewidth=0.25)

    # Draw locations

    for data_obs, start_date, end_date, obs_lat_name, obs_lon_name, markersize, color in zip(
        data_obss, start_dates, end_dates, obs_lat_names, obs_lon_names, markersizes, colors, strict=True
    ):
        label = (
            data_obs.time_series.attrs['long_name']
            + ',\n'
            + data_obs.time_series.attrs['region']
            + ', '
            + start_date.date().isoformat()
            + '$-$'
            + end_date.date().isoformat()
        )

        ax.plot(
            data_obs.time_series[obs_lon_name],
            data_obs.time_series[obs_lat_name],
            'o',
            markersize=markersize,
            label=label,
            color=color,
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

    gl.xlabel_style = {'rotation': 0, 'ha': 'center', 'va': 'top', 'size': 14}
    gl.ylabel_style = {'rotation': 0, 'ha': 'right', 'va': 'center', 'size': 14}

    # Add legend
    ax.legend(fontsize=legend_font_size)

    title = ''

    ax.set_title(title, fontsize=12)

    obs_descriptors = [
        data_obs.time_series.attrs['region'] + '_' + data_obs.time_series.attrs['name'] for data_obs in data_obss
    ]
    obs_descriptor = '.'.join(obs_descriptors)

    plot_name = obs_descriptor + '.stations_map.png'

    plot_path = plot_dir / plot_name

    plot_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(plot_path, bbox_inches="tight", dpi=600)

    plt.close(fig)

    if verbose:
        print('Created plot', plot_path)

    return plot_path


def line_plot_1d(
    x_size,
    y_size,
    x_label,
    y_label,
    title,
    x_datas,
    y_datas,
    data_labels,
    linecolors,
    linewidths,
    linestyles,
    x_min=None,
    x_max=None,
    y_min=None,
    y_max=None,
    legend_font_size=10,
    plot_path=None,
):
    '''
    Generate and optionally save a one-dimensional line plot for one or more datasets.

    Parameters
    ----------
    x_size : float
        Width of the figure in inches.
    y_size : float
        Height of the figure in inches.
    x_label : str
        Label for the x-axis.
    y_label : str
        Label for the y-axis.
    title : str
        Title of the plot.
    x_datas : list of array-like or array-like
        One or more sequences of x-axis values, one per dataset.
    y_datas : list of array-like or array-like
        One or more sequences of y-axis values corresponding to each entry in `x_datas`.
    data_labels : list of str or str
        Labels for each data series, used in the legend.
    linecolors : list of str or str
        Colors for each plotted line (any valid Matplotlib color specification).
    linewidths : list of float or float
        Line widths for each dataset.
    linestyles : list of str or str
        Line styles for each dataset (e.g., "-", "--", "-.", ":").
    x_min : float, optional
        Minimum limit for the x-axis. If None, determined automatically.
    x_max : float, optional
        Maximum limit for the x-axis. If None, determined automatically.
    y_min : float, optional
        Minimum limit for the y-axis. If None, determined automatically.
    y_max : float, optional
        Maximum limit for the y-axis. If None, determined automatically.
    legend_font_size : int, optional
        Font size for figure legend.
    plot_path : pathlib.Path or str, optional
        Path to save the generated plot image. If None, the plot is not saved.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If any of the input lists differ in length, as enforced by ``zip(..., strict=True)``.

    Notes
    -----
    - Inputs that are provided as single objects (rather than lists) are automatically
      wrapped into lists for consistent iteration.
    - All iterable inputs (`x_datas`, `y_datas`, `data_labels`, `linecolors`,
      `linewidths`, and `linestyles`) must have equal length.
    - A Matplotlib figure is created and each dataset is plotted with the corresponding
      visual properties. Axis labels, title, limits, and legend are applied as specified.
    - If ``plot_path`` is given, the figure is saved to that location; otherwise, it is
      only created and closed.
    - The figure is always closed after saving to free memory.
    '''

    # Normalize inputs that are supposed to be lists to lists

    if not isinstance(x_datas, list):
        x_datas = [x_datas]
    if not isinstance(y_datas, list):
        y_datas = [y_datas]
    if not isinstance(data_labels, list):
        data_labels = [data_labels]
    if not isinstance(linecolors, list):
        linecolors = [linecolors]
    if not isinstance(linewidths, list):
        linewidths = [linewidths]
    if not isinstance(linestyles, list):
        linestyles = [linestyles]

    # Plot

    fig, ax = plt.subplots(figsize=(x_size, y_size), dpi=300)

    for x_data, y_data, data_label, linecolor, linewidth, linestyle in zip(
        x_datas, y_datas, data_labels, linecolors, linewidths, linestyles, strict=True
    ):
        ax.plot(
            x_data,
            y_data,
            linewidth=linewidth,
            color=linecolor,
            linestyle=linestyle,
            label=data_label,
        )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    ax.set_title(title)
    ax.legend(title=None, fontsize=legend_font_size)

    if x_min is not None and x_max is not None:
        ax.set_xlim(x_min, x_max)

    if y_min is not None and y_max is not None:
        ax.set_ylim(y_min, y_max)

    if plot_path:
        plot_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(plot_path, bbox_inches='tight')

    plt.close()


def plot_evaluation_time_series(
    model_vs_obs_files: list[Path],
    plot_dir: Path,
    title: 'str | None' = None,
    t_min: 'datetime | None' = None,
    t_max: 'datetime | None' = None,
    x_size: 'float | None' = None,
    y_size: 'float | None' = None,
    colors: 'list[object] | None' = None,
    linethicknesses: 'list[float] | None' = None,
    linestyles: 'list[str] | None' = None,
    legend_font_size: int = 10,
    verbose: bool = False,
) -> list[Path]:
    '''
    Create time series evaluation plots from one or more Arotake model-vs-observation NetCDF datasets.

    Parameters
    ----------
    model_vs_obs_files : list[pathlib.Path]
        Paths to Arotake-produced model-vs-observations evaluation datasets that contain
        time series variables marked with the attribute ``description == 'time series'``
        as well as global attributes used for labeling (e.g., model name, region, observations).

    plot_dir : pathlib.Path
        Output directory where PNG plots are written.

    title : str | None, optional
        Figure title applied to all plots.

    t_min, t_max : datetime.datetime | None, optional
        X-axis (time) limits passed through to the plotting backend.

    x_size, y_size : float | None, optional
        Figure size in inches for width and height, respectively.

    colors : list[object] | None, optional
        Per-series color specifications forwarded to the plotting backend.

    linethicknesses : list[float] | None, optional
        Per-series line widths.

    linestyles : list[str] | None, optional
        Per-series line styles.

    legend_font_size : int, optional
        Font size for figure legend.

    verbose : bool, optional
        If True, print progress information.

    Returns
    -------
    list[pathlib.Path]
        File paths to the created plot images, in the same order as variables are processed.

    Notes
    -----
    This function does not alter the datasets and forwards plotting parameters unchanged
    to ``plotting.line_plot_1d``. It assumes that all provided datasets describe the same
    base variable, forecast initialization hour, and lead time metadata used for labeling.
    '''

    #
    # Load datasets
    #

    dss = [xr.load_dataset(model_vs_obs_file, decode_timedelta=True) for model_vs_obs_file in model_vs_obs_files]

    #
    # Base (model) variable
    #

    model_variable_name_ = [ds.attrs['model_variable'] for ds in dss]

    if len(set(model_variable_name_)) == 1:
        model_variable_name = model_variable_name_[0]
    else:
        raise Exception('Datasets contain different model (base) variables:', model_variable_name_)

    #
    # Forecast initialization and lead times
    #

    model_forecast_init_hours = [
        ds['model_forecast_init_time'].values[0].astype('datetime64[h]').astype(int) % 24 for ds in dss
    ]
    model_forecast_lead_hours = [int(round(ds['model_forecast_lead_time'].values[0] / np.timedelta64(1, 'h'))) for ds in dss]

    #
    # Time series variables
    #

    var_names_ = []

    for ds in dss:
        var_names = [var_name for var_name in ds.data_vars if ds[var_name].attrs.get('description') == 'time series']
        var_names_.append(var_names)

    # Test if the time series variables are present in all model-vs-obs files:
    if all(lst == var_names_[0] for lst in var_names_):
        # The datasets have identical time series variables
        var_names = var_names_[0]
    else:
        # The datasets have different time series variables
        print('Time series variables differ between model-vs-observation files:')
        for var_names, model_vs_obs_file in zip(var_names_, model_vs_obs_files, strict=True):
            print('File', model_vs_obs_file, 'contains the time series variables')
            for var_name in var_names:
                print('  ' + var_name)
        raise ValueError('Model-vs-observation files contain different time series variables.')

    if verbose:
        print()
        print('Plots will be created for the following time series variables:')
        print()
        for var_name in var_names:
            print(var_name, ':', dss[0][var_name].attrs.get('long_name'))
        print()

    time_series_plots = []

    for var_name in var_names:
        # Variable long name

        var_long_name_ = [ds[var_name].attrs['long_name'] for ds in dss]

        if len(set(var_long_name_)) == 1:
            var_long_name = var_long_name_[0]
        else:
            raise Exception(
                'Model-vs-observation files contain different variable long names for the variable',
                var_name,
                ':',
                var_long_name_,
            )

        # Units

        units_ = [ds[var_name].attrs['units'] for ds in dss]

        if len(set(units_)) == 1:
            units = units_[0]
        else:
            raise Exception('Model-vs-observation files contain different units for the variable', var_name, ':', units_)

        # Labels

        if 'Model-observation' in var_long_name:
            labels = [
                ds.model_region_name
                + ', '
                + ds.model
                + ', '
                + 'init hour = '
                + str(model_forecast_init_hour)
                + ' UTC , '
                + 'lead time = '
                + str(model_forecast_lead_hour)
                + ' h'
                + ' vs '
                + ds.observations
                for ds, model_forecast_init_hour, model_forecast_lead_hour in zip(
                    dss, model_forecast_init_hours, model_forecast_lead_hours, strict=True
                )
            ]
        elif 'Model' in var_long_name:
            labels = [
                ds.model_region_name
                + ', '
                + ds.model
                + ', '
                + 'init hour = '
                + str(model_forecast_init_hour)
                + ' UTC , '
                + 'lead time = '
                + str(model_forecast_lead_hour)
                + ' h'
                for ds, model_forecast_init_hour, model_forecast_lead_hour in zip(
                    dss, model_forecast_init_hours, model_forecast_lead_hours, strict=True
                )
            ]
        else:
            labels = [ds.model_region_name + ', ' + ds.observations for ds in dss]

        #
        # x- and y-values
        #

        x_label = 'Date (UTC)'
        y_label = var_long_name.replace('at reporting', '\nat reporting') + ' (' + units + ')'

        time = [ds['time'] for ds in dss]
        data = [ds[var_name] for ds in dss]

        if t_min is None:
            t_min = min([min(t) for t in time])
        if t_max is None:
            t_max = max([max(t) for t in time])

        y_min = None
        y_max = None

        # Other elements

        if x_size is None:
            x_size = 10
        if y_size is None:
            y_size = 4
        if colors is None:
            colors = [None for ds in dss]
        if linethicknesses is None:
            linethicknesses = [1.5 for ds in dss]
        if linestyles is None:
            linestyles = ['-' for ds in dss]

        regions = [ds.model_region_name for ds in dss]

        plot_file_path = plot_dir / Path('_'.join(regions) + '.' + model_variable_name + '.' + var_name + '.png')

        # Create plot

        line_plot_1d(
            x_size,
            y_size,
            x_label,
            y_label,
            title,
            time,
            data,
            labels,
            colors,
            linethicknesses,
            linestyles,
            x_min=t_min,
            x_max=t_max,
            y_min=y_min,
            y_max=y_max,
            legend_font_size=legend_font_size,
            plot_path=plot_file_path,
        )

        time_series_plots.append(plot_file_path)

        if verbose:
            print('Created plot', plot_file_path)

    return time_series_plots


def plot_evaluation_pdf(
    model_vs_obs_files: list[Path],
    plot_dir: Path,
    title: 'str | None' = None,
    x_size: 'float | None' = None,
    y_size: 'float | None' = None,
    colors: 'list[object] | None' = None,
    linethicknesses: 'list[float] | None' = None,
    linestyles: 'list[str] | None' = None,
    legend_font_size: int = 10,
    verbose: bool = False,
) -> list[Path]:
    '''
    Create probability density evaluation plots from one or more Arotake model-vs-observation NetCDF datasets.

    Parameters
    ----------
    model_vs_obs_files : list[pathlib.Path]
        Paths to Arotake-produced model-vs-observations evaluation datasets that contain
        time series variables marked with the attribute ``description == 'time series'``
        as well as global attributes used for labeling (e.g., model name, region, observations).

    plot_dir : pathlib.Path
        Output directory where PNG plots are written.

    title : str | None, optional
        Figure title applied to all plots.

    x_size, y_size : float | None, optional
        Figure size in inches for width and height, respectively.

    colors : list[object] | None, optional
        Per-probability density color specifications forwarded to the plotting backend.

    linethicknesses : list[float] | None, optional
        Per-probability density line widths.

    linestyles : list[str] | None, optional
        Per-probability density line styles.

    legend_font_size : int, optional
        Font size for figure legend.

    verbose : bool, optional
        If True, print progress information.

    Returns
    -------
    list[pathlib.Path]
        File paths to the created plot images, in the same order as variables are processed.

    Notes
    -----
    This function does not alter the datasets and forwards plotting parameters unchanged
    to ``plotting.line_plot_1d``. It assumes that all provided datasets describe the same
    base variable, forecast initialization hour, and lead time metadata used for labeling.
    '''

    #
    # Load datasets
    #

    dss = [xr.load_dataset(model_vs_obs_file, decode_timedelta=True) for model_vs_obs_file in model_vs_obs_files]

    #
    # Base (model) variable
    #

    model_variable_name_ = [ds.attrs['model_variable'] for ds in dss]

    if len(set(model_variable_name_)) == 1:
        model_variable_name = model_variable_name_[0]
    else:
        raise Exception('Datasets contain different model (base) variables:', model_variable_name_)

    #
    # Forecast initialization and lead times
    #

    model_forecast_init_hours = [
        ds['model_forecast_init_time'].values[0].astype('datetime64[h]').astype(int) % 24 for ds in dss
    ]
    model_forecast_lead_hours = [int(round(ds['model_forecast_lead_time'].values[0] / np.timedelta64(1, 'h'))) for ds in dss]

    #
    # Time series variables
    #

    var_names_ = []

    for ds in dss:
        var_names = [var_name for var_name in ds.data_vars if ds[var_name].attrs.get('description') == 'time series']
        var_names_.append(var_names)

    # Test if the time series variables are present in all model-vs-obs files:
    if all(lst == var_names_[0] for lst in var_names_):
        # The datasets have identical time series variables
        var_names = var_names_[0]
    else:
        # The datasets have different time series variables
        print('Time series variables differ between model-vs-observation files:')
        for var_names, model_vs_obs_file in zip(var_names_, model_vs_obs_files, strict=True):
            print('File', model_vs_obs_file, 'contains the time series variables')
            for var_name in var_names:
                print('  ' + var_name)
        raise ValueError('Model-vs-observation files contain different time series variables.')

    if verbose:
        print()
        print('Plots will be created for the following time series variables:')
        print()
        for var_name in var_names:
            print(var_name, ':', dss[0][var_name].attrs.get('long_name'))
        print()

    pdf_plots = []

    for var_name in var_names:
        # Calculate PDFs

        data_grids = []
        pdfs = []
        y_max = 0

        for ds in dss:
            # PDF function over the entire period using Gaussian kernel density estimation - filter out non-finite data (missing values)
            data = ds[var_name].values
            mask = np.isfinite(data)

            # Data grid for PDF calculation
            data_grid = np.mgrid[min(data[mask]) : max(data[mask]) : 100j]
            data_grids.append(data_grid)

            # PDF
            pdf = gaussian_kde(data[mask])

            # PDF on data grid
            density = pdf(data_grid)

            y_max = max(y_max, max(density))

            pdfs.append(density)

        # Variable long name

        var_long_name_ = [ds[var_name].attrs['long_name'] for ds in dss]

        if len(set(var_long_name_)) == 1:
            var_long_name = var_long_name_[0]
        else:
            raise Exception(
                'Model-vs-observation files contain different variable long names for the variable',
                var_name,
                ':',
                var_long_name_,
            )

        # Units

        units_ = [ds[var_name].attrs['units'] for ds in dss]

        if len(set(units_)) == 1:
            units = units_[0]
        else:
            raise Exception('Model-vs-observation files contain different units for the variable', var_name, ':', units_)

        # Labels

        if 'Model-observation' in var_long_name:
            labels = [
                ds.model_region_name
                + ', '
                + ds.model
                + ', '
                + 'init hour = '
                + str(model_forecast_init_hour)
                + ' UTC , '
                + 'lead time = '
                + str(model_forecast_lead_hour)
                + ' h'
                + ' vs '
                + ds.observations
                + ', '
                + str(ds['time'].values[0].astype('datetime64[D]'))
                + '$-$'
                + str(ds['time'].values[-1].astype('datetime64[D]'))
                for ds, model_forecast_init_hour, model_forecast_lead_hour in zip(
                    dss, model_forecast_init_hours, model_forecast_lead_hours, strict=True
                )
            ]
        elif 'Model' in var_long_name:
            labels = [
                ds.model_region_name
                + ', '
                + ds.model
                + ', '
                + 'init hour = '
                + str(model_forecast_init_hour)
                + ' UTC , '
                + 'lead time = '
                + str(model_forecast_lead_hour)
                + ' h'
                + ', '
                + str(ds['time'].values[0].astype('datetime64[D]'))
                + '$-$'
                + str(ds['time'].values[-1].astype('datetime64[D]'))
                for ds, model_forecast_init_hour, model_forecast_lead_hour in zip(
                    dss, model_forecast_init_hours, model_forecast_lead_hours, strict=True
                )
            ]
        else:
            labels = [
                ds.observations
                + ', '
                + str(ds['time'].values[0].astype('datetime64[D]'))
                + '$-$'
                + str(ds['time'].values[-1].astype('datetime64[D]'))
                for ds in dss
            ]

        #
        # x- and y-values
        #

        x_label = var_long_name + ' (' + units + ')'
        y_label = 'Probability density'

        x_min = None
        x_max = None
        y_min = 0
        y_max = 1.15 * y_max

        # Other elements

        if x_size is None:
            x_size = 6.5
        if y_size is None:
            y_size = 5
        if colors is None:
            colors = [None for ds in dss]
        if linethicknesses is None:
            linethicknesses = [1.5 for ds in dss]
        if linestyles is None:
            linestyles = ['-' for ds in dss]

        regions = [ds.model_region_name for ds in dss]

        plot_file_path = plot_dir / Path('_'.join(regions) + '.PDF.' + model_variable_name + '.' + var_name + '.png')

        # Create plot

        line_plot_1d(
            x_size,
            y_size,
            x_label,
            y_label,
            title,
            data_grids,
            pdfs,
            labels,
            colors,
            linethicknesses,
            linestyles,
            x_min=x_min,
            x_max=x_max,
            y_min=y_min,
            y_max=y_max,
            legend_font_size=legend_font_size,
            plot_path=plot_file_path,
        )

        pdf_plots.append(plot_file_path)

        if verbose:
            print('Created plot', plot_file_path)

    return pdf_plots
