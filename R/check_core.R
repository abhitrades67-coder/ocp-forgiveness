# Core qualification checks for the R port (cf. src/_tune.py).
suppressMessages({
  here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
  if (length(here) == 0 || here == "") here <- "R"
  source(file.path(here, "pk.R")); source(file.path(here, "drugs.R"))
  source(file.path(here, "hpo.R")); source(file.path(here, "simulate.R"))
})

DRUGS <- build_drugs()

# --- untreated cycle ---
t_end <- 120
tg <- seq(0, t_end, by = 0.05)
prof <- list(t = tg, c_prog = numeric(length(tg)), c_ee = numeric(length(tg)))
u <- hpo_simulate(prof, t_end, ec50_prog = 1e9, emax_prog = 0, hill_prog = 2,
                  ec50_ee = NA, emax_ee = 0, hill_ee = 2)
cat("UNTREATED: ovulations (days):", paste(round(u$ov_times, 1), collapse = ", "), "\n")
cat("  intervals (d):", paste(round(diff(u$ov_times), 1), collapse = ", "), "\n")
cat("  E2 range:", round(min(u$Y[iE2, ]), 1), "-", round(max(u$Y[iE2, ]), 1),
    "| Follicle max:", round(max(u$Y[iF, ]), 1),
    "| P4 max:", round(max(u$Y[iP4, ]), 1), "\n\n")

# --- perfect-use suppression per agent ---
for (k in c("LNG", "DRSP", "ETG", "LNGP", "DRSPP")) {
  drug <- DRUGS[[k]]
  if (drug$kind == "POP") sched <- build_pack_schedule(28, 0, n_packs = 5)
  else sched <- build_pack_schedule(21, 7, n_packs = 5)
  taken <- apply_adherence(sched)
  prof <- make_profile(drug, taken, t_end)
  ec_p <- ec50_prog(drug); ec_e <- if (drug$has_ee) ec50_ee(drug) else NA_real_
  s <- hpo_simulate(prof, t_end, ec_p, drug$emax_prog, drug$hill_prog,
                    ec_e, drug$ee_emax, drug$ee_hill)
  cat(sprintf("PERFECT %-4s: ovulations = %d | follicle max = %.1f | E2 max = %.1f | EC50_prog = %.4f\n",
              k, length(s$ov_times), max(s$Y[iF, ]), max(s$Y[iE2, ]), ec_p))
}
