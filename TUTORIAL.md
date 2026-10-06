# srflux tutorial

A start-to-finish guide: install the package, check it works, then run it on your
own data. No prior knowledge of the code is assumed. If you only want the API
reference, read the [README](README.md) instead.

**Contents**

1. [What the model actually does](#1-what-the-model-actually-does)
2. [What data you need](#2-what-data-you-need)
3. [Installing](#3-installing)
4. [Checking the install works](#4-checking-the-install-works)
5. [Running it on your own data](#5-running-it-on-your-own-data)
6. [Reading the output](#6-reading-the-output)
7. [Calibrating properly](#7-calibrating-properly)
8. [Common problems](#8-common-problems)

---

## 1. What the model actually does

Eddy covariance is widely regarded as the reference method for measuring
canopy–atmosphere fluxes, and it is what surface renewal is calibrated and
evaluated against. But it needs a sonic anemometer, which is expensive, and it
places real demands on the site — an adequate upwind fetch, careful maintenance.
Those constraints are why alternatives exist.

Surface renewal is one of them. It estimates the same flux from a single fast
thermometer, using no wind measurement at all.

The premise is that turbulence does not exchange heat with the surface smoothly.
Air parcels near the canopy are repeatedly replaced by parcels from above: a
parcel rests there long enough to exchange heat with the surface, and is then
swept away and replaced as a coherent eddy sweeps in and ejects it. A thermometer
sitting in that air records the cycle directly — a gradual warming, then a sharp
reset, over and over. On a temperature trace it looks like a row of ramps, or
sawtooth teeth.

How much heat is moving depends on two properties of those teeth: how tall they
are (the amplitude) and how much time passes between successive renewals.
Combine those with the heat capacity of air and a length scale and you have a
flux — with no vertical wind measurement anywhere in the calculation.

The idea is older than its use in micrometeorology. It began in 1930s–50s
chemical engineering, as a way of describing transfer across an interface as
parcels being renewed rather than heat diffusing through a static film. It
entered atmospheric science once ramp-shaped structures were observed in
temperature records over plant canopies and tied to the sweep–ejection cycle.

**Surface renewal is not calibration-free.** The number the ramps give you is the
right *shape* but the wrong *size*, because a single-point thermometer does not
resolve the whole volume of air being exchanged, and the idealised ramp is not
quite the real one. So it is scaled by one coefficient, `alpha`, found by
comparing against measured fluxes over a calibration period. You need a reference
flux at least once. What the method buys you is that afterwards, temperature
alone gives you the flux.

This package does that in four steps:

```
temperature  ->  find the ramps  ->  uncalibrated flux  ->  x alpha  ->  H
```

and then, optionally, one more: if you also have net radiation and soil heat
flux, the energy balance turns `H` into evapotranspiration.

**Two ways to find the ramps** are included, and you pick one:

| | how it works | when to use it |
|---|---|---|
| **SR-WL** (Haar) | Locates each ejection front directly and measures the temperature jump across it | The default. Makes no structural assumption about the ramp shape |
| **SR-VA** (Van Atta) | Infers amplitude and period from the block's structure functions, under an assumed ramp model | The established benchmark; use it to compare against published work |

Where the sweep–ejection cycle is well organised the two agree, which is a useful
check: if they disagree badly on your data, suspect the data before the method.

**How well does it work?** Tested against eddy covariance over grass, vineyard
and almond orchard, the wavelet detector matched or exceeded the established
air-temperature benchmark for sensible heat flux, and surface temperature
reproduced it within the eddy-covariance scatter. Agreement is best at the daily
scale — daily ET recovered with slopes of 0.95–1.01 and biases within
±0.10 mm d⁻¹. Individual half-hour estimates are noticeably noisier than that,
so treat the daily total as the reliable product.

---

## 2. What data you need

Two CSV files.

### File 1 — high-rate temperature

The signal the ramps are found in. **1 Hz or faster.** One row per sample.

```
TIMESTAMP,IRT
2023-05-11 00:00:00,10.740
2023-05-11 00:00:01,10.740
2023-05-11 00:00:02,10.698
```

- A timestamp column, and at least one temperature column in °C.
- Either an air temperature (fine-wire thermocouple) or a canopy surface
  temperature (infrared radiometer). Both work; they need different settings,
  covered below.
- Missing rows are fine — the loader puts the series back on a regular grid and
  marks the holes. You do not need to fill them yourself.
- `.csv.gz` is read directly; no need to unzip.

### File 2 — block-averaged reference

Used to calibrate `alpha`, and to check the answer. **One row per averaging
period**, normally 30 minutes, on the same clock as the temperature file.

```
TIMESTAMP,NETRAD,G,H,LE
2023-05-11 00:00:00,-94.9,-34.47,-30.43,14.31
2023-05-11 00:30:00,-94.47,-34.88,-53.25,27.71
```

- `H` — measured sensible heat flux [W m⁻²]. **Required**, this is what
  calibration fits against. Eddy covariance, ideally energy-balance corrected.
- `NETRAD` and `G` — net radiation and soil heat flux [W m⁻²]. Optional; needed
  only if you want evapotranspiration out the far end.
- `LE` is not used by the model; it is in the example file for comparison.

Both files ship as working examples in `examples/data/`, so you can look at the
exact format rather than guessing from this description.

### Once you are calibrated

You only need File 2 for the calibration period. Once you have an `alpha` for
your site and sensor, running on new data needs the temperature file alone —
which is the entire point of the method.

---

## 3. Installing

These commands assume macOS or Linux; Windows differences are noted.

### Step 1 — get the code

```bash
git clone https://github.com/nbambach/srflux.git
cd srflux
```

Everything below is run from inside that `srflux` directory.

### Step 2 — make a virtual environment

A virtual environment keeps this project's packages separate from everything else
on your machine, so installing srflux cannot break another project. Optional but
recommended.

```bash
python3 -m venv srflux_env
```

Activate it — **macOS / Linux**:

```bash
source srflux_env/bin/activate
```

**Windows** (PowerShell):

```powershell
srflux_env\Scripts\Activate.ps1
```

Your prompt should now start with `(srflux_env)`. You need to run that activate
command again in each new terminal session. To leave it, type `deactivate`.

### Step 3 — install srflux

```bash
pip install -e .
```

This installs the package *and* its dependencies (numpy, pandas). The `-e` means
"editable": Python reads the code from this folder, so if you change a source
file the change takes effect immediately with no reinstall. That is what you want
while working from a clone.

Two optional extras:

```bash
pip install -e ".[examples]"   # adds matplotlib and jupyter, for the plots and notebook
pip install -e ".[test]"       # adds pytest, for running the test suite
```

`pytest` is a test runner. You do not need it to use the model — it is there so
you can confirm the code behaves correctly on your machine, which is worth doing
once after installing. If you want it, install the `[test]` extra and run:

```bash
pytest
```

Every test should pass. If any fail, something is wrong with the install and the
results would not be trustworthy.

Requires Python 3.9 or newer. `requirements.txt` lists the dependencies if you
prefer to install those separately first.

---

## 4. Checking the install works

Before pointing anything at your own data, run the synthetic example:

```bash
python examples/quickstart.py
```

It builds a fake day of temperature data with a *known* flux hidden inside it,
runs the whole chain, and reports how closely it recovered the answer. No data
files needed. If this prints a table of alpha values and a daily ET, the install
is good.

---

## 5. Running it on your own data

Use `examples/run_my_data.py`. It is written to be copied and edited — all the
settings are in one clearly marked block at the top, and nothing below that block
needs changing.

As shipped it runs on the example day, so start by confirming that works:

```bash
python examples/run_my_data.py
```

Then copy it and edit for your data:

```bash
cp examples/run_my_data.py my_site.py
```

Open `my_site.py` and work down the `SETTINGS` block:

| setting | what to put |
|---|---|
| `SCALAR_CSV` | path to your temperature file |
| `SCALAR_COL` | the column name inside it, e.g. `"IRT"` |
| `REFERENCE_CSV` | path to your block-averaged file |
| `H_COL`, `NETRAD_COL`, `G_COL` | the column names in that file |
| `FS` | sampling rate in Hz, or `None` to work it out from the timestamps |
| `BLOCK_S` | averaging period in seconds; must match the reference file |
| `Z` | length scale in metres — see below |
| `DETECTOR` | `"haar"` or `"vanatta"` |
| `ALPHA` | `None` to fit it, or a number to apply an existing one |

**Choosing `Z`.** This is the length scale of the air being exchanged.

- Air temperature sensor → use the **measurement height**.
- Surface (radiometric) temperature → use the **canopy height**, because the air
  being renewed is the canopy layer.

`Z` and `alpha` are not separately identifiable — only their product affects the
flux — so choosing `Z` differently is absorbed by the calibration and does not
change your fluxes. It matters when comparing coefficients with someone else, who
must have made the same choice. Which convention is *correct* is a genuinely open
question in the literature: applied work usually scales by measurement height,
the original derivation ties it to canopy height, and there is a case for the
displacement height `z − d` above the roughness sublayer. Reporting the
height-invariant product `alpha × canopy height` sidesteps the argument.

**If you are using a surface temperature with the Van Atta detector**, set
`DETECTOR = "vanatta"` and leave the `period_mode="unit"` that the script already
uses for that case. The reason is in the README; the short version is that the
fitted ramp period is unreliable at 1 Hz on a radiometer, so it is better to drop
it and let `alpha` absorb it.

Then run it:

```bash
python my_site.py
```

---

## 6. Reading the output

The script prints a summary and writes two files.

**`srflux_output.csv`** — one row per averaging block:

| column | meaning |
|---|---|
| `timestamp` | start of the block |
| `count` | ramps detected |
| `amplitude` | ramp amplitude [°C] |
| `period` | mean time between ramps [s] |
| `F` | uncalibrated flux [W m⁻²] |
| `sign` | inferred flux direction, +1 up / −1 down |
| `H_reference` | the measured H it was compared against |
| `H_srflux` | **the calibrated result** [W m⁻²] |
| `LE_srflux` | latent heat by energy-balance residual [W m⁻²] |

**`srflux_output.png`** — the calibrated `H` against the measured `H` over the
day. This is the plot to look at first: if the two curves have the same shape,
the detector is working, and if they differ only by a constant factor, the
calibration is off rather than the physics.

The printed summary reports `alpha` and three skill numbers. **Read NSE, not r.**
`r` only measures whether the shape is right — an estimate that is correct in
shape but twice too large still scores `r = 1.0`. NSE penalises being the wrong
size. An NSE of 0 means you would have done just as well predicting the average.

Blank rows are normal: a block that was too incomplete, or in which no ramps rose
above the detection threshold, is reported as missing rather than guessed.

---

## 7. Calibrating properly

The one mistake that makes results look far better than they are: **fitting
`alpha` on the same blocks you then report skill for.** The script warns when you
do this. It is fine while setting things up, but the resulting NSE is not a real
measure of performance.

Do it properly like this:

1. Take a stretch of days where you have both temperature and measured `H`.
2. Run the script on those days with `ALPHA = None` to fit a coefficient.
3. Note the `alpha` it prints.
4. Set `ALPHA` to that number, and run on *different* days.

Practical points, in rough order of how much trouble they cause:

- **Use 10–21 days, spread out**, not one continuous block. Days close together
  share weather, so a continuous window converges three to five times more
  slowly for the same number of days.
- **Know how far a coefficient travels.** Reusing last year's `alpha` at the same
  site costs little — about 0.4–0.5 mm d⁻¹ of extra daily-ET error, close to the
  floor set by calibrating on the site itself. Borrowing from another site of the
  same crop roughly doubles that. Borrowing across crops does not work: median
  error around 1.2 mm d⁻¹, and the worst site pairs exceed 6 mm d⁻¹. Calibrate on
  your own site if you possibly can.
- **Fit separate coefficients per regime.** Daytime convection and the weak
  night-time flux do not share an `alpha`. `examples/data/make_ola_calibration.py`
  shows the pattern.
- **Check each day before pooling them.** The fit is weighted by `F²`, so a
  single day with a mis-scaled sensor has enormous leverage and can drag the
  pooled coefficient nearly to zero.
- **Watch for sensor drift.** Plot the median ramp amplitude per day. If it
  trends, `alpha` will move without the flux having changed, and your window is
  too long.

---

## 8. Common problems

**`no usable blocks`** — `BLOCK_S` may be longer than the data you have, or the
timestamp column is not being parsed. Print the loaded series and check the first
and last timestamps look right.

**`KeyError: no 'TIMESTAMP' column`** — your file uses a different name for the
time column. Set `TIME_COL` to match. The error lists the columns it did find.

**Most blocks come out blank** — usually a sampling-rate mismatch. If `FS` says
1 Hz but the file is 10 Hz, every block covers a tenth of the intended period.
Set `FS = None` to have it inferred from the timestamps.

**`alpha` looks wrong** — it is not an arbitrary tuning knob, and published
values give you a sanity check:

| your setup | expect roughly |
|---|---|
| air temperature, Van Atta | 0.24–0.34 |
| air temperature, Haar | 0.8–1.8 |
| surface temperature, Haar | 2.8–15.9, larger the further the sensor sits above the canopy |

A value far outside the row that matches your setup usually means a scaling error
upstream — check `Z`, and check the units of your temperature column. The large
values for surface temperature are expected rather than a symptom: the transport
length is the canopy height rather than the sensor height, the surface-temperature
ramp amplitude is damped by the surface-to-aerodynamic resistance ratio, and the
renewal timescale is prolonged by thermal relaxation of the canopy elements.

**The sign is wrong at dawn or dusk** — expected with the default settings. The
flip threshold defaults to zero, and near the `H → 0` transition the raw
skewness sign turns over too early. Fit it with `srflux.fit_skewness_sign` on a
calibration window and set `SIGN_TAU`.

**Everything works but ET looks wrong** — ET is a residual, so it inherits every
error in `Rn`, `G` and `H` at once. Check the energy balance closes on your
reference data before blaming the surface-renewal step.

---

Questions and bug reports: <https://github.com/nbambach/srflux/issues>
