# Reproducibility

**Manuscript:** Computational assessment of oral contraceptive forgiveness: an in silico quantitative systems pharmacology and PK-PD study of missed-dose tolerance for ovulation suppression

## Environment

- R version 4.5.2 (2025-10-31 ucrt), x86_64-w64-mingw32/x64, Windows 11 x64.
- Production integrator: a self-contained fixed-step fourth-order Runge-Kutta solver written in base R (`R/hpo.R`). The pharmacokinetic layer is analytic.
- Cross-validation engine: the compiled-C solver of the `rxode2` R package, with the state-triggered ovulation event applied by segment-and-restart.
- Plotting and reporting packages: `ggplot2`, `patchwork`, `scales`, `officer`, `sensitivity`, `boot` (see `R/README.md`).
- `deSolve` is deliberately not used in production: on this build and platform it was unstable under the repeated-call Monte-Carlo workload. This is stated in the manuscript.

## Running the pipeline

```
Rscript run_protocol.R          # full run, writes results_R/
```

On Windows, run the heavy analyses from a script file rather than with `Rscript -e`.

## Determinism

All randomness is seeded deterministically from a stable string hash of the scenario identifiers, not from the process-salted built-in hash, so repeated runs reproduce every number exactly. The one exception is the Morris elementary-effects design, which the `sensitivity` package generates without a fixed seed, so Table S10 reproduces as a ranking rather than as identical design points. Pack positions share one common-random-number virtual population: virtual subject *i* occupies row *i* at every position, which is what makes the paired bootstrap in the H2 analysis valid.

## Where each reported number comes from

| Manuscript item | Source file under `results_R/tables/` |
|---|---|
| Table 2, maximum tolerated missed pills with bootstrap CI | `max_tolerated_missed.csv` |
| Table 3, recovery actions | `recovery.csv` |
| Table 4, delayed-intake validation | `validation_delayed_intake.csv` |
| Table 5, guideline concordance | `guideline_concordance.csv` |
| Pack-position contrasts (H2) quoted in the Results and abstract | `h2_paired.csv` |
| Supplementary Table S1, solver cross-validation | `rxode2_crossvalidation.csv` |
| Supplementary Table S3, continuous forgiveness index (all pack positions) | `max_tolerated_missed.csv`, `forgiveness_index.csv` |
| Supplementary Table S4, covariate grid (CYP3A4 status, body weight) | `escape_probability_grid.csv`, rows with `analysis` in {`cyp`, `weight`, `main`} |
| Supplementary Table S5, Pearl-Index translation (all pack positions) | `pearl_index.csv` |
| Supplementary Table S6, pharmacokinetic validation | `validation_pk.csv` |
| Supplementary Table S7, perfect-use ovulation | `validation_pd.csv` |
| Supplementary Table S8, untreated-cycle validation | `cycle_validation.csv` |
| Supplementary Table S9, local one-at-a-time sensitivity | `sensitivity_local.csv` |
| Supplementary Table S10, Morris global sensitivity | `sensitivity_morris.csv` |
| Supplementary Table S11, escape-threshold robustness | `threshold_robustness.csv` |
| Levonorgestrel missed-at-start (Elomaa) check | `validation_missed_pill.csv` |
| Full per-cell escape grid | `escape_probability_grid.csv` |

Table 1 and Supplementary Table S2 list model inputs rather than outputs: the per-agent parameters in `R/drugs.R` and the hypothalamic-pituitary-ovarian parameters in `P_default()` in `R/hpo.R`.

## Figures

The TIFF files in `Figures/` and `Supplementary_Figures/` are lossless conversions of the master PNGs in `final files/`, which are themselves produced by `R/reporting.R` into `results_R/figures_hi/`. Each conversion was verified pixel-identical to its source.

## Sample sizes by analysis

Set in `run_protocol.R`. Each table and figure now states its own N.

| Analysis | Constant | N |
|---|---|---|
| Main missed-pill grid, perfect-use and missed-at-start validation | `N_MAIN` | 600 |
| CYP3A4 and body-weight covariates | `N_COV` | 250 |
| Lateness, regimen, first-pack and scattered-adherence analyses | `N_SEC` | 250 |
| Recovery-action comparison | `N_REC` | 250 |
| Local one-at-a-time sensitivity | `N_SENS` | 300 |
| Morris elementary effects, per model evaluation | `n_inner` in `morris_global()` | 60 |

`sample_population()` seeds on (agent, body weight, CYP3A4 status, tag), so the three
pack positions of one agent share a virtual population and support a paired bootstrap,
whereas a covariate arm is an independent draw from its reference and does not.

## Verified end-to-end, from a cleared state

The whole pipeline was re-run from cleared checkpoints and every output table diffed
against the shipped set. **All 17 tables reproduce exactly.** Three files initially
appeared to differ and did not: two differed only in row order and number formatting
introduced by targeted rerun scripts, and `h2_paired.csv` as shipped predated a code
change that added three columns (the shared columns were identical). Compared on a
composite key, `escape_probability_grid.csv` matched on all 255 rows with zero value
mismatches.

A full clean run takes **about 100 minutes** and needs at least two invocations: the
per-invocation wall-clock budget (`OCP_BUDGET_S`, default 1300 s) makes the script exit
gracefully with checkpoints saved, and relaunching resumes. The "3582 s" in the original
run log is cumulative across resumptions, not a single-shot runtime.

### Morris screening: seeded and widened after that check

The one table that did **not** reproduce was `sensitivity_morris.csv`. The design is
randomly generated, and at the r = 6 trajectories originally used the estimator was far
too noisy: SE(mu*) = sigma/sqrt(r) reached 13 percentage points, so levonorgestrel's
leading driver moved from mu* 28.7 to 8.3 between runs. The design is now generated from
a fixed seed and r raised to 40, which brings SE(mu*) to 0.1-3.8; Table S10 reports
SE(mu*) alongside mu* so the precision is explicit. Every number in the package is now
exactly reproducible.

## Model corrections applied during development

Two coding-level corrections were made to `R/` after the first full run, and only the
affected outputs were recomputed:

1. **Ethinylestradiol clearance** now uses its own CYP3A4-dependent fraction
   (`ee_fm_cyp3a4 = 0.5`), which was defined in `R/drugs.R` but never read; previously
   the progestin's fraction was applied to both compounds (`R/pk.R`). Because
   `cyp3a4_clearance_factor(f, "normal")` is exactly 1 and the shared clearance random
   effect is drawn once and reused, every analysis outside the CYP3A4 arm is
   bit-identical; this was asserted on five sampled main-grid cells before the rerun and
   passed. Only the 64 `analysis == "cyp"` rows of `escape_probability_grid.csv` changed.
   The effect is confined to the levonorgestrel combined pill, whose progestin and
   ethinylestradiol fractions differ (0.10 vs 0.50); induced perfect-use escape rose from
   6.4% to 28.0%. Figure 5 was regenerated and Supplementary Table S4 updated.

2. **Fertile window** for the cervical-barrier arm: a six-day window ending at ovulation
   was trialled in place of the fixed 14-day window running from the challenge dose. It
   was **not adopted**. For a suppressed combined-pill user the follicle peaks during the
   hormone-free interval, when the progestin concentration is zero, so a physiologically
   anchored window measures no barrier at all: levonorgestrel's barrier competence fell
   from 0.935 to 0.034 and its perfect-use Pearl Index rose to 4.13 against a published
   0.3, unreachable by any margin value, while drospirenone required a margin of 54.6
   against 2.6. The layer's agreement with published Pearl Indices therefore depends on
   averaging across active-pill days rather than over a fertile window. The 14-day window
   is retained and this result is reported as a limitation.

### Levonorgestrel's inducible fraction corrected

A third correction was made after an expert pharmacology review. `fm_cyp3a4` for
levonorgestrel (both the combined pill and the traditional progestin-only pill) was 0.10,
which contradicts the published PBPK fraction of 0.33-0.47 and the SmPC statement that
CYP3A4 is the main route. Under the four-fold induction factor of Eq. (19) the model
predicts a rifampicin AUC ratio of 1/(1 + 3*fm), so 0.10 gives 0.77 against an observed
0.43. The value is now **0.44**, which reproduces the observed ratio (0.431) and sits
inside the published range. Because the induction factor is exactly 1 at normal CYP3A4
status for any fm, only the 64 `analysis == "cyp"` rows changed; the regression guard on
three main-grid cells passed bit-identically before the rerun. Levonorgestrel's induced
perfect-use escape moved from 28.0% to 29.2%. Figure 5, Table 1 and Supplementary
Table S4 were regenerated.

The `fm_cyp3a4` values are lumped **inducible** fractions anchored to observed interaction
data, not measured metabolic fractions, and Table 1 now says so: drospirenone's 0.50
exceeds its in vitro CYP fraction (4-7%) but is needed to approach the observed
ketoconazole and rifampicin effects. A consequence is that the two combined pills no
longer span a wide metabolic range, which the Methods now state.

## What is calibrated and what is predicted

- Calibrated: the ovarian-dynamics parameters, to untreated-cycle behaviour; the progestin suppression margins, to published perfect-use ovulation rates; and the per-agent cervical-barrier ceilings, to published perfect-use Pearl Indices.
- Not calibrated, and therefore genuine out-of-sample checks: the pharmacokinetic disposition, the levonorgestrel missed-at-start escape curve, the two scheduled-delay ovulation-inhibition trials, and the concordance with FSRH and CDC missed-pill guidance.
- No missed-pill datum enters any calibration. This is what makes the forgiveness ranking a prediction rather than a restatement of the fit.
