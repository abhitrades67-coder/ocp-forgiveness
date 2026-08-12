# Validation datasets

Reference data the OCP-forgiveness model is validated against. No public **individual-level**
dataset of *ovulation under missed pills* exists for these agents, so model qualification is a
**predictive check against published summary statistics** (the appropriate standard for a
literature-based QSP model), supplemented by a **quantitative check of the untreated
hypothalamic-pituitary-ovarian (HPO) oscillator** against published natural-cycle hormone
reference values. The model's hormones are in arbitrary units, so the HPO check is on
**units-independent features** - cycle length, luteal-phase duration, and phase fold-amplitudes
(peri-ovulatory:follicular and luteal:follicular for E2, LH, P4) - i.e. whether the model
reproduces the *shape and amplitude* of the natural cycle.

## Files

| File | Contents | Consumed by |
|------|----------|-------------|
| `natural_cycle_hormones_Anckaert2021.csv` | Serum E2/LH/P4 across the natural cycle (medians + 5th-95th percentiles, by phase and sub-phase) in 85 healthy women | `R/validate_cycle.R` → `results_R/tables/cycle_validation.csv`, Figure S2 |
| `guideline_concordance_rules.csv` | Major OB/GYN missed/late-dose rules (FSRH, CDC) and the model's representation/result for each | Manuscript Guideline-Concordance table (Discussion) |
| `validation_targets_summary.csv` | Published PK, perfect-use ovulation, and Pearl-Index targets used as calibration/qualification anchors | `R/validate.R`, manuscript Validation section |

## Primary natural-cycle source

> Anckaert E, Jank A, Petzold J, Rohsmann F, Paris R, Renggli M, Schönfeld K, Schiettecatte J,
> Kriner M. **Extensive monitoring of the natural menstrual cycle using the serum biomarkers
> estradiol, luteinizing hormone and progesterone.** *Practical Laboratory Medicine* 2021;25:e00211.
> doi:[10.1016/j.plabm.2021.e00211](https://doi.org/10.1016/j.plabm.2021.e00211) · PMID 33869706.

N = 85 healthy women (18-37 y), ~3 serum samples/week over one cycle, Roche Elecsys (cobas e801).
Cycles normalised to 29 days with ovulation on day 15. Open access (CC BY-NC-ND); reference values
transcribed from the article tables. FSH was not reported in this dataset and is therefore validated
only qualitatively. Bibliographic metadata retrieved via PubMed.

## Guideline-concordance sources (retrieved June 2026)

- **FSRH** Clinical Guideline: *Combined Hormonal Contraception*, October 2023.
- **FSRH** Clinical Guideline: *Progestogen-only Pills*, August 2022 (amended April 2026).
- **CDC** *U.S. Selected Practice Recommendations for Contraceptive Use*, 2024 - Recommended Actions
  After Late or Missed Combined Oral Contraceptives.

## PK/PD target provenance

See `../../PROVENANCE_AND_VALIDATION.md` (sources S1-S6, PubMed-verified) for the PK, perfect-use
ovulation, and missed-pill targets summarised in `validation_targets_summary.csv`.
