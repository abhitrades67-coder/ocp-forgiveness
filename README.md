# Oral contraceptive forgiveness: model and analysis pipeline

R implementation of an integrated population-pharmacokinetic and quantitative systems
pharmacology (QSP) model of the hypothalamic-pituitary-ovarian axis, used to quantify how
many missed or delayed active pills different oral contraceptives tolerate before
ovulation suppression is lost.

This repository contains everything needed to reproduce every number reported in the
accompanying paper.

> **Paper:** *Computational assessment of oral contraceptive forgiveness: an in silico
> quantitative systems pharmacology and PK-PD study of missed-dose tolerance for ovulation
> suppression.* [journal, year, DOI to be added on acceptance]

## What the model does

Five progestin formulations - the levonorgestrel and drospirenone combined pills, the
desogestrel (etonogestrel) pill, the traditional levonorgestrel-only pill and the
drospirenone-only pill - are simulated across three pack positions and 0 to 7 consecutive
missed active pills in a virtual population.

A one-compartment pharmacokinetic layer with first-order absorption drives sigmoid-Emax
suppression of the estradiol-triggered luteinising-hormone surge by progestin, and of
follicle-stimulating hormone by ethinylestradiol, within a six-state ODE model of the axis
in which ovulation is a discrete event fired by a surge-permission signal.

## Reproducing the results

```
Rscript run_protocol.R      # writes results_R/tables/
```

**A clean run takes about 100 minutes and needs at least two invocations.** The script
carries a per-invocation wall-clock budget (`OCP_BUDGET_S`, default 1300 s) and exits
gracefully with checkpoints saved; relaunching resumes where it stopped. This is expected
behaviour, not a hang. Set `OCP_BUDGET_S` higher to run longer per invocation.

On Windows, run the heavy analyses from a script file rather than with `Rscript -e`.

Requires R 4.5.2 with `parallel`, `ggplot2`, `patchwork`, `scales`, `officer`,
`sensitivity`, `boot` and `rxode2`. Exact versions are in `results_R/sessionInfo.txt`.

## Determinism

All Monte-Carlo sampling is seeded deterministically from a stable string hash of the
scenario identifiers, so repeated runs reproduce every number exactly. The Morris
elementary-effects design is also seeded (`seed = 20240601` in `R/sensitivity.R`).

This has been verified end to end: the whole pipeline was re-run from cleared checkpoints
and every output table diffed against the committed set. See `REPRODUCIBILITY.md`.

Pack positions share one common-random-number virtual population - virtual subject *i*
occupies row *i* at every position - which is what makes the paired bootstrap valid.
Covariate and recovery arms are seeded separately and are **not** nested subsets of the
reference; those comparisons are unpaired.

## Layout

| Path | Contents |
|---|---|
| `R/hpo.R` | Six-state HPO axis, RK4 integrator, ovulation event rule, default parameters |
| `R/pk.R` | Analytic one-compartment PK with first-order absorption |
| `R/drugs.R` | Per-agent parameters and the CYP3A4 clearance relationship |
| `R/simulate.R` | Scenario construction (pack structure, missed and delayed doses) |
| `R/forgiveness.R` | Virtual population, escape probability, bootstrap endpoints |
| `R/sensitivity.R` | Local one-at-a-time and Morris global screening |
| `R/recovery.R` | Optimal recovery-action comparison |
| `R/translate.R` | Cervical-mucus arm and the illustrative Pearl-Index mapping |
| `R/validate*.R` | Untreated-cycle, PK, PD and rxode2 cross-validation |
| `run_protocol.R` | Driver; writes every table under `results_R/tables/` |
| `data/validation/` | Published summary statistics used as validation targets |
| `results_R/tables/` | Every result table reported in the paper |

`R/hpo_rxode2.R` and `R/validate_rxode2.R` re-implement the model on the compiled-C solver
of `rxode2` as an independent cross-check of the production integrator.

## Scope and limitations

This model is calibrated to **rank** agents and pack positions, not to predict failure
rates. Absolute escape probabilities are illustrative. In particular:

- Ovulation suppression is a surrogate for contraceptive failure.
- Against the pivotal scheduled-delay trials the model over-predicts absolute post-delay
  escape ten- to nineteen-fold, and its delayed-intake arm delays seven consecutive pills
  where those trials delayed three or four non-consecutively.
- The `fm_cyp3a4` values are lumped **inducible** fractions anchored to observed
  interaction data, not measured metabolic fractions.
- The Pearl-Index layer is a calibrated illustration, not an independent translation.

The paper's Discussion sets these out in full.

## Not included

The manuscript-generation scripts (`build_submission.R`, `make_manuscript.R`) are omitted:
they embed prose from an earlier draft of the manuscript and are not needed to
reproduce any result.

## Licence

MIT - see `LICENSE`.

## Citation

See `CITATION.cff`. Please cite both the software (via its Zenodo DOI) and the paper.
