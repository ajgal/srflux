#!/usr/bin/env python
"""Run srflux on your own CSV data. Copy this file and edit the SETTINGS block.

    python examples/run_my_data.py

As shipped it runs on the example day in examples/data/, so you can check the
whole chain works before pointing it at anything of your own. Change the paths
and column names in SETTINGS below and it runs on your data instead.

What it does, in order:

    1. read a high-rate temperature CSV and a block-averaged reference CSV
    2. cut the temperature into averaging blocks
    3. detect ramps in each block and turn them into an uncalibrated flux
    4. fit the calibration coefficient alpha against the measured H
    5. apply alpha, get the flux direction, close the energy balance for ET
    6. write a per-block CSV and print a summary

What it writes:

    srflux_output.csv   one row per block: ramp statistics, calibrated H, LE, and
                        the reference H it was compared against
    srflux_output.png   H time series against the reference (only if matplotlib
                        is installed; pip install "srflux[examples]")
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from srflux import (HaarDetector, VanAttaDetector, align_reference, blocks_from_series,
                    daily_et, fit_and_score, latent_heat_residual, load_flux_csv,
                    load_scalar_csv, prepare_block, ramp_flux, sensible_heat,
                    sign_from_skewness)

# ----------------------------------------------------------------------------
# SETTINGS -- this is the part you edit
# ----------------------------------------------------------------------------

HERE = Path(__file__).parent

# Your high-rate temperature file. One timestamp column plus one or more data
# columns, 1 Hz or faster. A .csv.gz is read directly, no need to unzip.
SCALAR_CSV = HERE / "data" / "ola_2023-05-11_1hz.csv.gz"
SCALAR_COL = "IRT"          # which column holds the temperature [deg C]
TIME_COL = "TIMESTAMP"      # timestamp column name, in both files

# Your block-averaged reference file: one row per averaging period, holding the
# measured sensible heat flux used to calibrate, and Rn/G if you want ET.
REFERENCE_CSV = HERE / "data" / "ola_2023-05-11_flux30.csv"
H_COL = "H"                 # measured sensible heat flux [W m-2]
NETRAD_COL = "NETRAD"       # net radiation [W m-2]; set to None to skip ET
G_COL = "G"                 # soil heat flux [W m-2]; set to None to skip ET

# Measurement setup.
FS = 1.0                    # sampling frequency [Hz]; None infers it from the timestamps
BLOCK_S = 1800.0            # averaging period [s], must match the reference file
Z = 5.5                     # length scale [m]: measurement height for an air
                            # temperature, canopy height for a surface temperature

# Which detector. "haar" is the more robust default; "vanatta" fits the classical
# structure-function ramp instead.
DETECTOR = "haar"

# Calibration. Leave as None to fit alpha from this dataset. Set a number to
# apply a coefficient fitted elsewhere -- which is the honest choice when you
# want to SCORE the result, since fitting and scoring on the same blocks always
# flatters the method.
ALPHA = None

# Flux direction. The detectors measure ramp size, not which way the heat goes.
# "skewness" infers it from the temperature itself; "reference" borrows the sign
# of the measured H, which is only for checking the rest of the chain.
SIGN_SOURCE = "skewness"
SIGN_LAG_S = 3.0            # increment lag for the skewness test [s]
SIGN_TAU = 0.0              # flip threshold; fit it with srflux.fit_skewness_sign

OUTPUT_CSV = Path("srflux_output.csv")
OUTPUT_PNG = Path("srflux_output.png")

# ----------------------------------------------------------------------------
# From here down you should not need to change anything.
# ----------------------------------------------------------------------------


def build_detector():
    """Make the detector named in SETTINGS."""
    if DETECTOR == "haar":
        return HaarDetector(fs=FS or 1.0)
    if DETECTOR == "vanatta":
        # period_mode="unit" is the safe choice for a 1 Hz radiometric surface
        # temperature, whose fitted ramp period is too noisy to divide by
        return VanAttaDetector(fs=FS or 1.0, lag="chen", period_mode="unit")
    raise ValueError(f"DETECTOR must be 'haar' or 'vanatta', not {DETECTOR!r}")


def main() -> None:
    # --- 1. load ------------------------------------------------------------
    scalar = load_scalar_csv(SCALAR_CSV, value_col=SCALAR_COL, time_col=TIME_COL, fs=FS)
    reference = load_flux_csv(REFERENCE_CSV, time_col=TIME_COL)
    print(f"loaded {len(scalar)} samples of {SCALAR_COL} "
          f"({scalar.index[0]} to {scalar.index[-1]})")
    print(f"loaded {len(reference)} reference rows: {list(reference.columns)}")

    # --- 2. cut into blocks -------------------------------------------------
    stamps, blocks = blocks_from_series(scalar, block_s=BLOCK_S, fs=FS)
    if not blocks:
        raise SystemExit("no usable blocks -- check BLOCK_S and that the file has data")
    print(f"{len(blocks)} blocks of {len(blocks[0])} samples")

    # --- 3. detect ramps, one block at a time -------------------------------
    detector = build_detector()
    rate = FS or (len(blocks[0]) / BLOCK_S)
    rows = []
    for t, raw in zip(stamps, blocks):
        # prepare_block clips out-of-range values, fills short gaps, and rejects
        # a block that is too incomplete to trust
        qc = prepare_block(raw, fs=rate)
        if not qc.ok:
            rows.append(dict(timestamp=t, count=0, amplitude=np.nan, period=np.nan,
                             F=np.nan, sign=np.nan))
            continue

        res = detector.detect(qc.values)
        if not res.valid:
            rows.append(dict(timestamp=t, count=0, amplitude=np.nan, period=np.nan,
                             F=np.nan, sign=np.nan))
            continue

        # the two detectors use different flux forms: count/block for the front
        # picker, amplitude/period for Van Atta
        if detector.name == "haar":
            F = ramp_flux(res.amplitude, count=res.count, z=Z, block_s=BLOCK_S)
        else:
            F = ramp_flux(res.amplitude, period=res.period, z=Z, block_s=BLOCK_S)

        sign = sign_from_skewness(qc.values, fs=rate, tau=SIGN_TAU, lag_s=SIGN_LAG_S)
        rows.append(dict(timestamp=t, count=res.count, amplitude=res.amplitude,
                         period=res.period, F=float(F), sign=sign))

    out = pd.DataFrame(rows).set_index("timestamp")

    # --- 4. calibrate -------------------------------------------------------
    out["H_reference"] = align_reference(out.index, reference, H_COL)
    if ALPHA is None:
        # fit on unsigned flux against |H|: the sign is applied afterwards, so
        # letting it into the fit would flip the stable blocks twice
        fit = fit_and_score(out["F"].to_numpy(), np.abs(out["H_reference"].to_numpy()))
        alpha = fit.alpha
        print(f"\nfitted alpha = {alpha:.5g} on {fit.n} blocks "
              f"(r = {fit.r:.2f}, NSE = {fit.nse:.2f}, RMSE = {fit.rmse:.1f} W m-2)")
        print("  NOTE: alpha was fitted on these same blocks, so this score is "
              "optimistic.\n        Set ALPHA to a coefficient from other days to "
              "score honestly.")
    else:
        alpha = float(ALPHA)
        print(f"\napplying alpha = {alpha:.5g} from SETTINGS")

    # --- 5. calibrated flux, direction, and ET ------------------------------
    if SIGN_SOURCE == "reference":
        sign = np.sign(out["H_reference"].to_numpy())
    else:
        sign = out["sign"].to_numpy()
    out["H_srflux"] = sensible_heat(out["F"].to_numpy(), alpha, sign=sign)

    have_energy = NETRAD_COL in reference.columns and G_COL in reference.columns
    if have_energy:
        rn = align_reference(out.index, reference, NETRAD_COL)
        g = align_reference(out.index, reference, G_COL)
        out["LE_srflux"] = latent_heat_residual(rn, g, out["H_srflux"].to_numpy())
        et = daily_et(out["LE_srflux"].to_numpy(), block_s=BLOCK_S)
        print(f"daily ET from srflux  = {et:.2f} mm")
        if H_COL in reference.columns:
            le_ref = latent_heat_residual(rn, g, out["H_reference"].to_numpy())
            print(f"daily ET from measured H = {daily_et(le_ref, block_s=BLOCK_S):.2f} mm")
    else:
        print(f"skipping ET: need {NETRAD_COL!r} and {G_COL!r} in the reference file")

    # --- 6. write out -------------------------------------------------------
    out.to_csv(OUTPUT_CSV)
    print(f"\nwrote {OUTPUT_CSV}  ({len(out)} rows)")
    make_plot(out)


def make_plot(out: pd.DataFrame) -> None:
    """Save an H comparison figure, if matplotlib is available."""
    try:
        import matplotlib
        matplotlib.use("Agg")          # write a file without needing a display
        import matplotlib.pyplot as plt
    except ImportError:
        print('no figure written; pip install "srflux[examples]" for matplotlib')
        return

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(out.index, out["H_reference"], "k-", lw=1.2, label="measured H")
    ax.plot(out.index, out["H_srflux"], "o-", ms=3, lw=1.0, label="srflux H")
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_ylabel("sensible heat flux [W m$^{-2}$]")
    ax.set_xlabel("time")
    ax.legend(frameon=False)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=150)
    print(f"wrote {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
