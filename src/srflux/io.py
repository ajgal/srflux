"""Reading your own data in: CSV -> clean regular series -> blocks.

Surface renewal is computed one averaging block at a time, so before any detector
runs your data has to be turned into a list of equal-length, regularly-sampled
blocks. Doing that by hand is where most of the mistakes happen, so it lives here
rather than in the examples.

Two files are needed for a full run:

  high-rate scalar   one temperature column at 1 Hz or faster, with a timestamp.
                     This is what the ramps are detected in -- a fine-wire
                     thermocouple, or a radiometric surface temperature.

  reference flux     one row per averaging block (usually 30 min), carrying the
                     measured H used to calibrate alpha, and optionally NETRAD and
                     G if you want the energy-balance route to ET.

The cleaning step that matters is REINDEXING onto a regular time grid. A logger
file with missing rows is shorter than the period it covers, and a block cut by
row count would then span more time than it should and quietly mis-scale the ramp
rate. Reindexing turns those missing rows into NaN, which :func:`~srflux.preprocess.prepare_block`
then either interpolates or rejects.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def infer_fs(index: pd.DatetimeIndex) -> float:
    """Sampling frequency [Hz] from the median timestamp spacing.

    The median, not the mean: a single large gap in the file would drag a mean
    spacing far off and silently halve the inferred rate.

    Examples
    --------
    >>> idx = pd.date_range("2023-05-11", periods=100, freq="1s")
    >>> infer_fs(idx)
    1.0
    """
    if len(index) < 2:
        raise ValueError("need at least two timestamps to infer the sampling rate")
    dt = np.median(np.diff(index.values).astype("timedelta64[ns]").astype(float)) / 1e9
    if dt <= 0:
        raise ValueError("timestamps are not increasing")
    return float(1.0 / dt)


def load_scalar_csv(path, value_col=None, time_col="TIMESTAMP", fs=None,
                    value_range=None) -> pd.Series:
    """Load one high-rate scalar column and put it on a regular time grid.

    Parameters
    ----------
    path : str, Path or file-like
        Anything ``pandas.read_csv`` accepts, including a ``.csv.gz``.
    value_col : str, optional
        Column holding the scalar. Defaults to the first non-timestamp column,
        which is the common case for a one-sensor export.
    time_col : str
        Timestamp column. Parsed with pandas' inference; ISO-like strings such as
        ``2023-05-11 00:00:00`` are safest.
    fs : float, optional
        Sampling frequency [Hz]. Inferred from the timestamps when omitted.
    value_range : tuple, optional
        Values outside this range become NaN here, before blocking. Use it for a
        logger sentinel that would otherwise survive as a real number
        (``-8190``, ``-3.5e8``). ``prepare_block`` applies its own range later.

    Returns
    -------
    pandas.Series
        Indexed by timestamp on a complete regular grid, gaps present as NaN.

    Examples
    --------
    >>> import io
    >>> rows = ["2023-05-11 00:00:0%d,%.1f" % (t, 11.0 + t)
    ...         for t in (0, 1, 2, 4, 5)]          # 00:00:03 is missing
    >>> s = load_scalar_csv(io.StringIO("TIMESTAMP,FWT\\n" + "\\n".join(rows) + "\\n"))
    >>> len(s)                      # the missing row is restored on the grid
    6
    >>> bool(s.isna().sum() == 1)
    True
    """
    df = pd.read_csv(path)
    if time_col not in df.columns:
        raise KeyError(f"no {time_col!r} column; found {list(df.columns)}")
    if value_col is None:
        others = [c for c in df.columns if c != time_col]
        if not others:
            raise ValueError("file has a timestamp column and nothing else")
        value_col = others[0]
    if value_col not in df.columns:
        raise KeyError(f"no {value_col!r} column; found {list(df.columns)}")

    # t = pd.to_datetime(df[time_col], errors="coerce") # original line
    t = pd.to_datetime(df[time_col], format="%Y-%m-%d %H:%M:%S", errors="coerce")
    # errors="coerce" turns junk into NaT/NaN rather than raising, so one bad row
    # in a million-row logger file does not abort the load
    v = pd.to_numeric(df[value_col], errors="coerce")
    s = pd.Series(v.to_numpy(), index=pd.DatetimeIndex(t)).sort_index()
    s = s[s.index.notna()]
    # a logger that restarts mid-file can repeat a timestamp; keep the first
    s = s[~s.index.duplicated(keep="first")]
    if s.empty:
        raise ValueError("no valid rows after parsing timestamps")

    if value_range is not None:
        lo, hi = value_range
        s = s.where((s >= lo) & (s <= hi))

    rate = infer_fs(s.index) if fs is None else float(fs)
    # reindex onto the complete grid: this is what makes a block's sample count
    # and its duration agree even when rows are missing from the file
    grid = pd.date_range(s.index[0], s.index[-1], freq=pd.Timedelta(seconds=1.0 / rate))
    return s.reindex(grid)


def load_flux_csv(path, time_col="TIMESTAMP", columns=None) -> pd.DataFrame:
    """Load the block-averaged reference file (H, and optionally NETRAD and G).

    Parameters
    ----------
    columns : sequence of str, optional
        Restrict to these columns. Missing ones raise rather than being filled,
        so a typo in a column name fails loudly instead of producing all-NaN.

    Returns
    -------
    pandas.DataFrame
        Indexed by timestamp, sorted, numeric.

    Examples
    --------
    >>> import io
    >>> csv = io.StringIO("TIMESTAMP,NETRAD,G,H\\n2023-05-11 12:00:00,600,70,250\\n")
    >>> float(load_flux_csv(csv)["H"].iloc[0])
    250.0
    """
    df = pd.read_csv(path)
    if time_col not in df.columns:
        raise KeyError(f"no {time_col!r} column; found {list(df.columns)}")
    t = pd.to_datetime(df[time_col], errors="coerce")
    # float throughout: an all-integer column would otherwise come back as int64
    # and propagate a dtype that cannot hold the NaN a missing block needs
    df = (df.drop(columns=[time_col])
            .apply(pd.to_numeric, errors="coerce")
            .astype(float))
    df.index = pd.DatetimeIndex(t)
    df = df[df.index.notna()].sort_index()
    if columns is not None:
        missing = [c for c in columns if c not in df.columns]
        if missing:
            raise KeyError(f"missing column(s) {missing}; found {list(df.columns)}")
        df = df[list(columns)]
    return df


def blocks_from_series(series: pd.Series, block_s: float = 1800.0, fs=None,
                       min_valid_frac: float = 0.5):
    """Cut a regular series into averaging blocks.

    Blocks are aligned to whole multiples of ``block_s`` past midnight, which is
    what makes them line up with an eddy-covariance file written on the same
    clock. A partial block at either end is dropped rather than padded.

    Returns
    -------
    (timestamps, blocks)
        ``timestamps`` is a DatetimeIndex of block START times -- matching the
        convention of the reference files this package ships -- and ``blocks`` a
        list of 1-D arrays, each ``block_s * fs`` samples long.

    Examples
    --------
    >>> idx = pd.date_range("2023-05-11 00:00:00", periods=7200, freq="1s")
    >>> s = pd.Series(np.arange(7200.0), index=idx)
    >>> ts, blocks = blocks_from_series(s, block_s=1800)
    >>> len(blocks), len(blocks[0]), str(ts[1])
    (4, 1800, '2023-05-11 00:30:00')
    """
    if not isinstance(series.index, pd.DatetimeIndex):
        raise TypeError("series must be indexed by timestamp")
    rate = infer_fs(series.index) if fs is None else float(fs)
    n = int(round(block_s * rate))
    if n < 2:
        raise ValueError("block_s * fs must be at least 2 samples")

    # snap the start forward to the next block boundary so block edges are on the
    # clock (00:00, 00:30, ...) rather than wherever the file happens to begin
    start = series.index[0].ceil(pd.Timedelta(seconds=block_s))
    s = series.loc[series.index >= start]

    stamps, blocks = [], []
    values = s.to_numpy(dtype=float)
    index = s.index
    for i in range(0, len(values) - n + 1, n):
        chunk = values[i:i + n]
        # drop a block that is mostly missing before it reaches the detector; the
        # same test runs again in prepare_block, but failing here keeps the
        # returned lists aligned with what actually got processed
        if np.isfinite(chunk).sum() < min_valid_frac * n:
            continue
        stamps.append(index[i])
        blocks.append(chunk)
    return pd.DatetimeIndex(stamps), blocks


def align_reference(stamps: pd.DatetimeIndex, reference: pd.DataFrame,
                    column: str = "H") -> np.ndarray:
    """Line a reference column up with block start times.

    Returns NaN where the reference has no row for a block, so the calibration
    simply drops that block instead of silently pairing it with the wrong half
    hour.

    Examples
    --------
    >>> stamps = pd.DatetimeIndex(["2023-05-11 12:00", "2023-05-11 12:30"])
    >>> ref = pd.DataFrame({"H": [250.0]}, index=pd.DatetimeIndex(["2023-05-11 12:00"]))
    >>> align_reference(stamps, ref).tolist()
    [250.0, nan]
    """
    if column not in reference.columns:
        raise KeyError(f"no {column!r} column; found {list(reference.columns)}")
    return reference[column].reindex(stamps).to_numpy(dtype=float)
