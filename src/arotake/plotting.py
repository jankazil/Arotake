'''
Plotting utilities for time-series diagnostics.

This module provides:
- plot_df_timeseries: Generate and save one figure per statistic (column) by overlaying
  the same statistic across multiple pandas.DataFrame inputs that share a common schema.
- pandas_series_plottable: Heuristic check for whether a pandas.Series can be plotted
  directly by Matplotlib without pre-conversion.

Assumptions shared by callers:
- Each DataFrame in a plotting call uses its first column as the time axis with datetime
  values and carries time unit metadata in "series.attrs['units']".
- Each DataFrame provides "attrs['model_variable_long_name']" for titling and a
  'model_forecast_lead_time' column with timedelta-like values.
'''

from datetime import datetime, timedelta
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pandas.api.types as ptypes
from cartopy.mpl.ticker import LatitudeFormatter, LatitudeLocator, LongitudeFormatter, LongitudeLocator

matplotlib.use("Agg")  # Important to avoid runaway memory use upon creating plots repeatedly.


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
        * <column_label> is the statistic name,
        * <suffix> is the underscore-joined per-DataFrame suffix list:
            - Without "legend_strings":
              "<HH:MM:SS><time_units>_f<lead_hours>h"
              (spaces removed from the valid-time string)
            - With "legend_strings[j]":
              "<legend_strings[j] with spaces→_>.<HH:MM:SS><time_units>f<lead_hours>h"
              (spaces removed from the valid-time string; no underscore before 'f').

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

        forecast_valid_time_string = time.iloc[0].strftime("%H:%M:%S") + ' ' + time.attrs['units']

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
        ax.set_ylabel(time_series.attrs['description'] + ' (' + time_series.attrs['units'] + ')')

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

    # Make sure that all lists have the same length:

    assert (
        len(data_obss) == len(start_dates) == len(end_dates) == len(obs_lat_names) == len(obs_lon_names)
    ), 'Arguments are lists of different length. Aborting.'

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

    # Draw locations

    markersizes = list(range(len(data_obss), 0, -1))

    for data_obs, start_date, end_date, obs_lat_name, obs_lon_name, markersize in zip(
        data_obss, start_dates, end_dates, obs_lat_names, obs_lon_names, markersizes, strict=True
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

    # Add legend
    ax.legend()

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
    ax.legend(title=None)

    if x_min is not None and x_max is not None:
        ax.set_xlim(x_min, x_max)

    if y_min is not None and y_max is not None:
        ax.set_ylim(y_min, y_max)

    if plot_path:
        plot_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(plot_path, bbox_inches='tight')

    plt.close()
