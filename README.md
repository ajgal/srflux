# srflux

[![tests](https://github.com/nbambach/srflux/actions/workflows/tests.yml/badge.svg)](https://github.com/nbambach/srflux/actions/workflows/tests.yml)

**Sensible heat flux from a high-rate temperature series, without vertical wind measurements.**

Eddy covariance is widely regarded as the reference method for canopy–atmosphere flux
measurement (Baldocchi 2014; Foken 2017), but its cost and its site requirements often make
alternatives necessary. Surface renewal (SR) is one: it infers turbulent fluxes from the
high-frequency fluctuations of a scalar alone (Paw U et al. 1995, 2005).

The physical premise is that air parcels near the surface are repeatedly replaced by parcels
from above. A parcel resides near the canopy long enough to exchange heat, water vapour or
CO₂, and turbulence then carries it away. In a scalar trace the cycle appears as a
characteristic ramp — a gradual warming or cooling terminated by a sharp reset — imprinted as
coherent eddies sweep in and eject the modified parcel. Given the amplitude of the scalar
change and the timescale between successive renewals, the corresponding flux can be estimated
without the direct vertical-wind measurement that eddy covariance requires.

The idea predates its use in micrometeorology by decades, originating in mid-twentieth-century
chemical-engineering treatments of interfacial mass transfer as turbulent parcel renewal
rather than diffusion through a static film (Higbie 1935; Danckwerts 1951). Gao et al. (1989)
linked ramp structures over plant canopies to the sweep–ejection cycle; Paw U and Brunet
(1991) and Paw U et al. (1995) formalised SR as a flux-estimation method; Snyder et al. (1996)
extended it to latent heat as a residual of the surface energy balance.

SR is not calibration-free: both detectors here return a flux of the correct form but
uncalibrated magnitude, closed by a coefficient `alpha` fitted against a reference flux. What
the method buys is that once `alpha` is known, the scalar series alone suffices.

This package implements the chain end to end: find the ramps, convert them to a flux,
calibrate, and optionally close the energy balance to get evapotranspiration.

**What to expect.** Evaluated against eddy covariance at three California sites — a reference
grass, vineyard and almond orchard — the wavelet detector matched or exceeded the established
air-temperature benchmark for H, and the surface-temperature retrieval reproduced H within the
eddy-covariance scatter. Agreement is strongest at the daily scale, where daily ET was
recovered with slopes of 0.95–1.01 and biases within ±0.10 mm d⁻¹. Half-hourly estimates are
correspondingly noisier; see [Citation](#citation).

**New here? Read the [tutorial](TUTORIAL.md)** — installation, the data you need, and a
worked run on your own files. The rest of this README is the reference.

| | method | what it measures |
|---|---|---|
| **SR-WL** | Haar wavelet front picker | counts microfronts and measures their amplitude directly |
| **SR-VA** | Van Atta structure-function cubic | fits an idealised ramp to S2/S3/S5 |

Both produce an *uncalibrated* flux, closed by one dimensionless coefficient `alpha`.

```python
from srflux import HaarDetector, VanAttaDetector, prepare_block, ramp_flux, fit_and_score

block = prepare_block(temperature_1hz, fs=1.0)                  # clip, gap-fill, validate
wl = HaarDetector(fs=1.0).detect(block.values)                  # SR-WL
va = VanAttaDetector(fs=1.0, lag="chen").detect(block.values)   # SR-VA

F_wl = ramp_flux(wl.amplitude, count=wl.count, z=2.0)           # uncalibrated, W m-2
F_va = ramp_flux(va.amplitude, period=va.period, z=2.0)

fit = fit_and_score(F_wl_all_blocks, H_reference)               # alpha, r, NSE, RMSE, bias
H_sr = fit.apply(F_wl_all_blocks)
```

## Install

```bash
git clone https://github.com/nbambach/srflux.git
cd srflux

python3 -m venv srflux_env          # optional but recommended
source srflux_env/bin/activate      # Windows: srflux_env\Scripts\Activate.ps1

pip install -e .                    # the package and its dependencies
```

Optional extras: `pip install -e ".[examples]"` adds matplotlib and jupyter for the plots and
the notebook; `pip install -e ".[test]"` adds pytest so you can run `pytest` and confirm the
install behaves correctly. Requires Python ≥ 3.9, numpy and pandas — also listed in
[`requirements.txt`](requirements.txt).

The [tutorial](TUTORIAL.md#3-installing) walks through the same steps with more explanation.

## What data you need

Two CSV files, both shown as working examples in `examples/data/`.

**High-rate temperature**, 1 Hz or faster, one row per sample. Air temperature (fine wire) or
canopy surface temperature (radiometer) both work. Missing rows are handled — the loader
restores the regular time grid rather than silently shortening blocks.

```
TIMESTAMP,IRT
2023-05-11 00:00:00,10.740
2023-05-11 00:00:01,10.740
```

**Block-averaged reference**, one row per averaging period (usually 30 min), on the same
clock. `H` is required and is what `alpha` is fitted against; `NETRAD` and `G` are needed only
for evapotranspiration.

```
TIMESTAMP,NETRAD,G,H,LE
2023-05-11 00:00:00,-94.9,-34.47,-30.43,14.31
```

Loading and blocking are in `srflux.io`, so you do not have to write that yourself:

```python
from srflux import load_scalar_csv, load_flux_csv, blocks_from_series, align_reference

series = load_scalar_csv("my_temps.csv.gz", value_col="IRT")   # regular grid, gaps as NaN
stamps, blocks = blocks_from_series(series, block_s=1800)      # aligned to the clock
H_ref = align_reference(stamps, load_flux_csv("my_flux.csv"), "H")
```

Once `alpha` is known for your site and sensor, the reference file is no longer needed — the
temperature file alone gives you the flux, which is the point of the method.

## The three examples

| file | what it is for | needs your data? |
|---|---|---|
| `examples/quickstart.py` | Confirms the install works. Builds synthetic data with a known flux and checks the chain recovers it. Run this first. | no |
| `examples/run_my_data.py` | **The one to copy.** Reads your CSVs and runs the full chain, writing a per-block CSV and a figure. All settings in one block at the top. | yes (ships working defaults) |
| `examples/ola_one_day.ipynb` | Explains the method visually — ramp anatomy, day and night behaviour, the route to ET — on one real day from an almond orchard, with α loaded from a calibration fitted on other days. | no |

The notebook is for *understanding* the method; `run_my_data.py` is for *using* it. Both run
the same library code.

`run_my_data.py` writes `srflux_output.csv` (one row per block: ramp statistics, calibrated H,
LE, and the reference it was compared against) and `srflux_output.png` (calibrated against
measured H). `quickstart.py` prints its results rather than writing files.

## The two detectors

Both answer the same question — how large are the ramps and how often do they arrive — by
different routes. SR-VA infers the two from the block's structure functions under an assumed
ramp model; SR-WL locates each front directly and measures it, without those structural
assumptions. Where the sweep–ejection cycle is well organised the two agree, which is evidence
that the wavelet is capturing the same physical property. Start with SR-WL.

**SR-WL — Haar (`srflux.detectors.haar`).** A surface-renewal ramp is, to first order, a
step: a gradual approach terminated by an abrupt front. A step is not a narrow-band
oscillation — it is a localised discontinuity, spread across frequencies but concentrated in
time — which is what makes a time–scale method appropriate and a Fourier decomposition not.

At a fixed scale the continuous wavelet transform reduces to a convolution, so the wavelet
acts simultaneously as a **matched filter** (largest output where the local signal resembles
the wavelet) and a **localised differencing operator** (its zero mean annihilates any local
constant or trend). The Haar wavelet is itself a step — the difference of two adjacent block
means — and so is the natural matched template for an ejection front. Four properties follow:

- **Step-matched shape.** The coefficient is the temperature jump across the front, so
  `a = median |c|` is the ramp amplitude directly, in kelvin, with no further scaling.
- **Minimal support.** Haar has the shortest support of any wavelet, giving the sharpest
  temporal localisation and well-separated coefficient peaks — which is what makes reliable
  front counting and a short (15 s) dead-time de-duplication possible.
- **One vanishing moment.** ∫ψ = 0 removes local constants, complementing the 300 s detrend
  and leaving the detector insensitive to slow drift.
- **Real-valued and non-oscillatory.** One sharp peak per front. Smooth or oscillatory
  wavelets — Mexican hat, complex Morlet — spread each front across several alternating-sign
  coefficients, blurring its location.

The implementation does not scan a continuum of scales. Ramps carry a characteristic
ejection-front duration, so a single dyadic scale matched to it is used as a one-scale filter.
This places SR-WL among objective, tunable alternatives to fixed-window conditional sampling
such as VITA (Blackwelder and Kaplan 1976) and WAG (Schols 1984), and to structure-function
ramp models (Van Atta 1977; Chen et al. 1997), following Collineau and Brunet (1993).

The threshold is **sigma-relative** (`k · std(coef)`, default `k = 0.75`), not absolute — a
canopy skin temperature has ramps an order of magnitude smaller than air temperature, and a
fixed threshold does not travel between sites or scalars. The 32 s kernel sits in a broad
flat optimum (20–48 s); neither default is worth tuning per site.

**SR-VA — Van Atta (`srflux.detectors.vanatta`).** Solves `A³ + pA + q = 0` with
`p = 10·S2 − S5/S3`, `q = 10·S3` for the ramp amplitude and `tau = A²r/S2` for the period.
The root is chosen to match `sign(S3)`.

Lag selection dominates the result. A fixed 1 s lag suits a 20 Hz sonic but collapses on a
radiometric surface temperature, whose ramps carry almost no 1 s signal. `lag="chen"` adapts
per block via the first global maximum of `|S3(r)|/r`.

On a low-rate skin temperature also **drop the period**: `period_mode="unit"` sets τ = 1 and
uses the amplitude alone, `F = ρ·cp·z·A`, folding the ramp duration into α. The fitted τ on a
1 Hz IRT can swing by an order of magnitude between blocks (15–324 s observed), so dividing
by it adds noise rather than information. With τ discarded, choose the lag on a calibration
window rather than with the Chen criterion, which targets the period you have just dropped.

## Calibration

The detectors return a flux of the correct form but uncalibrated magnitude, because a point
sensor resolves only part of the renewed volume and the idealised ramp is not the real one.
One coefficient closes the gap, fitted against a reference flux. Everything below is about
fitting it so that it still holds on days you did not fit on.

**Expected magnitudes.** `alpha` is not an arbitrary tuning knob and a value far outside these
ranges usually indicates an error upstream:

| scalar and detector | α |
|---|---|
| air temperature, SR-VA | 0.24–0.34 |
| air temperature, SR-WL | 0.8–1.8 |
| radiometric surface temperature, SR-WL | 2.8–15.9, scaling with `z/h` |

The air-temperature values sit alongside the published range — Fischer et al. (2023) report
0.52–0.69 unstable and 0.15–0.36 stable, and the Shapland et al. (2014) meta-analysis clusters
near 0.3 (0.28–0.38; pooled mean 0.41). SR-WL exceeds SR-VA on the same air-temperature data
chiefly because the two return different ramp periods τ, in part because the Van Atta model
assumes instantaneous termination.

Surface temperature needs an order of magnitude more, for three compounding reasons: the
transport length is canopy height `h` rather than sensor height `z` (`z/h` reaches 6 over
grass); the surface-temperature amplitude is damped by the ratio of surface to aerodynamic
resistance; and τ is prolonged by thermal relaxation of stagnant canopy elements. That
height dependence is a predictable consequence of the surface energy balance rather than an
artefact, so α carries physical information about the canopy.

```
SR-WL:  F = rho·cp·z · N·A / block        SR-VA:  F = rho·cp·z · A / tau
        H = alpha · F
```

`alpha` is fitted through the origin, `alpha = Σ(F·H)/Σ(F²)`, against a reference flux —
ideally energy-balance-closure-corrected eddy covariance — **per regime**, so that a
coefficient fitted on daytime convection is not applied to the weak nocturnal flux and the
ramp direction never enters the fit. `examples/data/make_ola_calibration.py` shows the
pattern: fit on a window of days, publish only the coefficients, apply them elsewhere.

Five things that matter in practice:

- **`alpha` and the length scale are not separately identifiable** — only their product enters
  the flux. For an air temperature the scale is the measurement height; for a surface
  temperature the renewed volume is the canopy layer, so canopy height is the physical choice.
  The ambiguity is longstanding: applied work generally scales by `z`, while the original
  column-budget derivation ties it to `h` (Paw U et al. 1995), and Castellví (2004) argues for
  `z − d` above the roughness sublayer. Two studies that choose differently cannot compare
  alphas — compare the height-invariant `alpha·h`, which collapsed a six-site spread from 4.6×
  to 1.7× and made vineyards and almonds overlap.
- **Transferability is directional.** Year to year at one site is close to free (≈0.4–0.5 mm
  d⁻¹ added to the daily-ET error, near the own-calibration floor). Between sites of the same
  crop roughly doubles that. Across crops it is not viable — median 1.2 mm d⁻¹, worst pairs
  above 6 mm d⁻¹.
- **Budget 10–21 days** spread across conditions, not one continuous campaign: a held-out
  test needed that much before the daily-ET error settled, and a continuous window converges
  three to five times more slowly because days inside it see the same weather.
- **Screen the calibration set per day before pooling.** A through-origin fit weights by F²,
  so one day on which the sensor is mis-scaled carries enormous leverage and can drag the
  pooled coefficient to nearly zero.
- **Bound the window by sensor drift, not the calendar.** Check the median ramp amplitude per
  day; a drifting sensor moves α with no change in the flux.
- **Score with NSE** (`srflux.scores`), not r. r is blind to scale, so an estimate with the
  right shape but the wrong magnitude still scores well.

## Flux direction

A detector measures how big the ramp is, not which way the heat is going. That has to come
from somewhere else, and there are two options depending on what you have.

`sign_from_skewness` is the temperature-only option: a warm ramp rises gradually and
collapses sharply, so the increment skewness is negative. **Fit its two parameters together
on a calibration window** with `fit_skewness_sign`:

- the **increment lag** must match the timescale of the fronts. The 1 s default suits a fine
  wire; a slower skin temperature needs several seconds, and a skin channel that appears to
  carry no direction information usually just needs a longer lag.
- the **flip threshold** is not zero — near the daytime H → 0 transition the raw sign flips
  early, and a small positive τ corrects it.

`sign_from_stability` uses `−sign(zeta)` where an eddy-covariance system is present.

## Layout

```
src/srflux/
  io.py                CSV loading, regular-grid cleaning, blocking
  preprocess.py        block QC, detrending, increment skewness
  detectors/haar.py    SR-WL front picker
  detectors/vanatta.py SR-VA cubic + Chen lag selection
  flux.py              uncalibrated flux, through-origin alpha, skill scores, LE and ET
  sign.py              flux-direction conventions
  synthetic.py         ramp generator for tests and sanity checks
tests/                 tests on signals whose answer is known by construction
examples/              synthetic quickstart, editable user script, real-data notebook
TUTORIAL.md            start-to-finish guide for new users
```

## Citation

The method implemented here, and the evaluation the numbers in this README come from:

> Bambach, N.E., McElrone, A.J., Paw U, K.T., Kustas, W.P., Knipper, K.R., Parry, C.,
> Hipps, L.E., Prueger, J.H., Alfieri, J.G., McKee, L.G., Gal, A., Nocco, M.A., Castro, S.J.,
> Tolentino, P.M. *A novel wavelet-based surface renewal analysis for sensible heat flux
> density estimations from radiometric surface temperature.* Submitted to *Agricultural and
> Forest Meteorology*, 2026.

To cite the software itself, see [CITATION.cff](CITATION.cff) or GitHub's "Cite this
repository" button.

## References

**Surface renewal**

- Higbie (1935) *Trans. Am. Inst. Chem. Eng.* **31**, 365–389; Danckwerts (1951) *Ind. Eng.
  Chem.* **43**, 1460–1467 — parcel-renewal origins in interfacial mass transfer.
- Gao, Shaw et al. (1989), in *Boundary Layer Studies and Applications*, Springer, 349–377 —
  organised structure in flow within and above a canopy.
- Paw U & Brunet (1991); Paw U, Qiu, Su, Watanabe & Brunet (1995) *Agric. For. Meteorol.*
  **74**, 119–137 — SR formalised for scalar fluxes.
- Snyder, Spano & Paw U (1996) *Boundary-Layer Meteorol.* **77**, 249–266 — sensible and
  latent heat flux density.
- Paw U, Snyder, Spano & Su (2005), in *Micrometeorology in Agricultural Systems*, ASA — SR
  in practice.
- Fischer et al. (2023) *Agric. For. Meteorol.* — merging flux-variance with SR; recent
  calibration ranges.

**Ramp detection**

- Van Atta (1977) *Arch. Mech.* **29**, 161–171 — structure-function ramp model.
- Chen, Novak, Black & Lee (1997) *Boundary-Layer Meteorol.* **84**, 99–124 — ramp model with
  finite microfront time; the lag criterion.
- Haar (1910) — the wavelet.
- Collineau & Brunet (1993) *Boundary-Layer Meteorol.* **65**, 357–379; Brunet & Collineau
  (1994) — wavelet detection of coherent motions and ramp fronts.
- Blackwelder & Kaplan (1976) *J. Fluid Mech.* **76**, 89–112 — VITA; Schols (1984)
  *Boundary-Layer Meteorol.* **29**, 39–58 — WAG: the fixed-window conditional sampling that
  wavelet detection replaces.

**Calibration and length scale**

- Spano, Snyder, Duce et al. (1997) *Agric. For. Meteorol.* **86**, 259–271 — the effective
  factor varies with height over one canopy.
- Castellví (2004) *Water Resour. Res.* **40**, W05201 — combining SR with similarity theory;
  referencing to `z − d`.
- Shapland, McElrone, Snyder & Paw U (2012b) *Boundary-Layer Meteorol.* **145**, 5–25 —
  two-scale scalar ramps; Shapland, Snyder, McElrone et al. (2014) *Agric. For. Meteorol.*
  **189**, 36–47 — frequency-response compensation and α convergence.

**Eddy covariance**

- Baldocchi (2014); Aubinet et al. (2012); Foken (2017) — the reference method and its
  requirements.

## License

MIT — see [LICENSE](LICENSE).
