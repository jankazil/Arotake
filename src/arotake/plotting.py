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

import matplotlib.pyplot as plt
import pandas as pd
import pandas.api.types as ptypes


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
