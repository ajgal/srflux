"""Tests for CSV loading and blocking."""
import io

import numpy as np
import pandas as pd
import pytest

from srflux.io import (align_reference, blocks_from_series, infer_fs, load_flux_csv,
                       load_scalar_csv)


def scalar_csv(rows, header="TIMESTAMP,FWT"):
    return io.StringIO(header + "\n" + "\n".join(rows) + "\n")


def one_hz(n, start="2023-05-11 00:00:00", skip=()):
    """n seconds of 1 Hz rows, optionally omitting some of them."""
    t0 = pd.Timestamp(start)
    return [f"{t0 + pd.Timedelta(seconds=i)},{20.0 + 0.01 * i}"
            for i in range(n) if i not in skip]


def test_infer_fs_reads_the_median_spacing():
    idx = pd.date_range("2023-05-11", periods=50, freq="100ms")
    assert infer_fs(idx) == pytest.approx(10.0)


def test_infer_fs_ignores_a_single_large_gap():
    # one 60 s hole must not drag the inferred rate off 1 Hz
    idx = pd.DatetimeIndex(list(pd.date_range("2023-05-11", periods=30, freq="1s"))
                           + [pd.Timestamp("2023-05-11 00:01:30")])
    assert infer_fs(idx) == pytest.approx(1.0)


def test_missing_rows_are_restored_as_nan():
    s = load_scalar_csv(scalar_csv(one_hz(20, skip=(5, 6, 7))))
    assert len(s) == 20                       # full grid, not the 17 rows present
    assert int(s.isna().sum()) == 3


def test_duplicate_timestamps_keep_the_first():
    rows = one_hz(5)
    rows.insert(3, rows[2])                   # logger restart repeats a stamp
    s = load_scalar_csv(scalar_csv(rows))
    assert len(s) == 5
    assert not s.index.duplicated().any()


def test_unsorted_rows_are_sorted():
    rows = one_hz(10)
    shuffled = [rows[i] for i in (4, 0, 9, 2, 1, 8, 3, 7, 5, 6)]
    s = load_scalar_csv(scalar_csv(shuffled))
    assert s.index.is_monotonic_increasing
    assert len(s) == 10


def test_unparseable_rows_do_not_abort_the_load():
    rows = one_hz(10)
    rows[4] = "not-a-timestamp,20.5"
    rows[6] = "2023-05-11 00:00:06,not-a-number"
    s = load_scalar_csv(scalar_csv(rows))
    assert len(s) == 10                       # grid still spans the full period
    assert np.isnan(s.iloc[6])                # the bad value became NaN


def test_value_range_masks_logger_sentinels():
    rows = one_hz(10)
    rows[3] = "2023-05-11 00:00:03,-8190"
    s = load_scalar_csv(scalar_csv(rows), value_range=(-40, 75))
    assert np.isnan(s.iloc[3])
    assert int(s.isna().sum()) == 1


def test_value_col_defaults_to_the_first_data_column():
    s = load_scalar_csv(scalar_csv(one_hz(5)))
    assert s.notna().all()


def test_missing_columns_raise_rather_than_silently_emptying():
    with pytest.raises(KeyError):
        load_scalar_csv(scalar_csv(one_hz(5)), value_col="NOPE")
    with pytest.raises(KeyError):
        load_scalar_csv(scalar_csv(one_hz(5)), time_col="NOPE")


def test_flux_csv_is_float_even_when_all_values_are_integers():
    csv = io.StringIO("TIMESTAMP,H\n2023-05-11 12:00:00,250\n")
    assert load_flux_csv(csv)["H"].dtype == float


def test_flux_csv_column_typo_raises():
    csv = io.StringIO("TIMESTAMP,H\n2023-05-11 12:00:00,250\n")
    with pytest.raises(KeyError):
        load_flux_csv(csv, columns=["H", "NETRAD"])


def test_blocks_are_aligned_to_the_clock():
    # start at 00:07 -- the first whole block must begin at 00:30, not 00:07
    idx = pd.date_range("2023-05-11 00:07:00", periods=7200, freq="1s")
    s = pd.Series(np.zeros(7200), index=idx)
    stamps, blocks = blocks_from_series(s, block_s=1800)
    assert str(stamps[0]) == "2023-05-11 00:30:00"
    assert all(len(b) == 1800 for b in blocks)


def test_partial_trailing_block_is_dropped():
    idx = pd.date_range("2023-05-11 00:00:00", periods=2000, freq="1s")
    s = pd.Series(np.zeros(2000), index=idx)
    _, blocks = blocks_from_series(s, block_s=1800)
    assert len(blocks) == 1                   # 200 leftover samples are not padded


def test_mostly_empty_blocks_are_skipped():
    idx = pd.date_range("2023-05-11 00:00:00", periods=3600, freq="1s")
    v = np.zeros(3600)
    v[:1700] = np.nan                         # first block is 94 % missing
    stamps, blocks = blocks_from_series(pd.Series(v, index=idx), block_s=1800)
    assert len(blocks) == 1
    assert str(stamps[0]) == "2023-05-11 00:30:00"


def test_align_reference_pairs_on_timestamp_not_position():
    stamps = pd.DatetimeIndex(["2023-05-11 12:00", "2023-05-11 12:30",
                               "2023-05-11 13:00"])
    ref = pd.DataFrame({"H": [250.0, 300.0]},
                       index=pd.DatetimeIndex(["2023-05-11 12:00", "2023-05-11 13:00"]))
    out = align_reference(stamps, ref)
    assert out[0] == 250.0
    assert np.isnan(out[1])                   # no 12:30 row -> NaN, not 300
    assert out[2] == 300.0


def test_round_trip_on_the_shipped_example_day():
    from pathlib import Path
    data = Path(__file__).resolve().parent.parent / "examples" / "data"
    s = load_scalar_csv(data / "ola_2023-05-11_1hz.csv.gz", value_col="IRT")
    ref = load_flux_csv(data / "ola_2023-05-11_flux30.csv")
    stamps, blocks = blocks_from_series(s, block_s=1800)
    assert len(blocks) == 48
    assert len(blocks[0]) == 1800
    assert np.isfinite(align_reference(stamps, ref, "H")).all()
