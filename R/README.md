# OCP-Forgiveness - R implementation (journal pipeline)

Reproducible in-silico QSP/PK-PD assessment of oral-contraceptive forgiveness.
This R pipeline implements the model and the full analysis.

## Requirements
- R ≥ 4.5 (developed on 4.5.2, Windows).
- Packages: `parallel` (base), `boot`, `ggplot2`, `patchwork`, `scales`,
  `officer`, `flextable`, `jsonlite`, `sensitivity` (Morris global SA).
  `deSolve` is **not** required - the production HPO ODE uses a self-contained
  pure-R RK4 integrator (deSolve's compiled solver was unstable under the
  repeated-call Monte-Carlo workload on the build used).
  `rxode2` (optional) provides a protocol-named **alternative ODE engine**
  (`R/hpo_rxode2.R`); its compiled-C solver runs cleanly here and reproduces the
  pure-R engine (cross-validated by `R/validate_rxode2.R`). Requires Rtools.

Install any missing packages:
```r
install.packages(c("boot","ggplot2","patchwork","scales","officer",
                   "flextable","jsonlite","sensitivity"))
```

## Run
From the project root (`OCP forgiveness/`):
```sh
Rscript run_protocol.R            # full run  (N = 600/cell)
NQUICK=1 Rscript run_protocol.R   # fast smoke test (small N)
Rscript R/reporting.R             # figures (results_R/figures/*.png) + REPORT.md
Rscript R/make_manuscript.R       # results_R/MANUSCRIPT.docx
Rscript R/validate_rxode2.R       # optional: rxode2-vs-RK4 engine cross-validation
```
`R/validate_rxode2.R` writes `results_R/tables/rxode2_crossvalidation.csv` and
`results_R/figures/figS1_rxode2_crossvalidation.png`; the manuscript picks these up
as Table S1 / Fig S1 if present.
`run_protocol.R` writes all CSV tables to `results_R/tables/`, a serialized
`results_R/bundle.rds` (consumed by the reporting and manuscript scripts), and a
live progress log to `results_R/run_progress.txt`.

Re-lock the calibration targets (perfect-use ovulation rates, follicular recovery):
```sh
NCAL=400 Rscript R/calibrate.R    # writes results_R/calibration.json
```

## Reproducibility
- All randomness is seeded deterministically (`stable_seed()` djb2 → base-R
  Mersenne-Twister), independent of N, so populations are nested and every number
  reproduces bit-for-bit on a given platform.
- **Always run via script files, not `Rscript -e '...'`** on Git-Bash/Windows:
  multiline `-e` strings are mangled at the shell→exe boundary and crash R.
- Pinned environment captured in `results_R/sessionInfo.txt`.

## File map
| File | Role |
|---|---|
| `R/pk.R` | one-compartment analytic PK + dosing engine |
| `R/drugs.R` | agent definitions, EC50 anchoring, cervical-mucus arm |
| `R/hpo.R` | HPO-axis ODE + pure-R RK4 integrator with ovulation events (production engine) |
| `R/hpo_rxode2.R` | same HPO model on rxode2's compiled-C solver (protocol-named alternative; segment-and-restart events) |
| `R/validate_rxode2.R` | cross-validates rxode2 vs pure-R RK4 (Table S1 / Fig S1) |
| `R/simulate.R` | single-subject PK-PD scenario |
| `R/forgiveness.R` | virtual population, MC escape, bootstrap max-tolerated, MC index |
| `R/validate.R` | PK + PD predictive checks (calibration vs independent) |
| `R/sensitivity.R` | local OAT tornado + Morris global screening |
| `R/recovery.R` | take-now / double / skip recovery-action analysis |
| `R/translate.R` | cervical-mucus arm + OR→Pearl-Index with r≈0.52 uncertainty |
| `R/calibrate.R` | locks parameters to documented targets |
| `run_protocol.R` | orchestrator |
| `R/reporting.R` | ggplot2 figures + REPORT.md |
| `R/make_manuscript.R` | officer/flextable → MANUSCRIPT.docx |
