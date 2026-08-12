# =============================================================================
# calibrate.R -- Lock model parameters to documented qualification targets.
# (COC perfect-use ovulation ~0-1%; follicular
# recovery so a normal 7-day HFI yields no ovulation and escape emerges only as
# the pill-free interval extends; Rice 1999 / Elomaa 1998 targets.)
#
# This is a MEASUREMENT + check harness: it reports the perfect-use ovulation
# rate per agent and the LNG missed-at-start (extended-HFI) escape curve, against
# the targets, so margin_prog / emax_prog / HPO params in drugs.R / hpo.R can be
# locked. Writes results_R/calibration.json.
# =============================================================================
suppressMessages({
  source("R/pk.R"); source("R/drugs.R"); source("R/hpo.R")
  source("R/simulate.R"); source("R/forgiveness.R"); source("R/parallel_util.R")
  library(jsonlite)
})

NCAL <- as.integer(Sys.getenv("NCAL", "300"))
DRUGS <- build_drugs()
cl <- make_model_cluster()
on.exit(stopCluster(cl), add = TRUE)

# ---- Targets (documented) ----
TARGETS <- list(
  LNG  = c(0, 2.0),   # perfect-use COC ovulation ~1-2% (~98-99% suppression)
  DRSP = c(0, 2.0),
  ETG  = c(0, 3),     # Rice 1999: desogestrel inhibition >97%
  LNGP = c(40, 60),   # Rice 1999: LNG-POP partial inhibition (ovulation arm)
  DRSPP= c(0, 3)      # drospirenone-only reliably inhibits ovulation
)

cat(sprintf("== Perfect-use ovulation rate (N=%d) vs target ==\n", NCAL))
pu <- list()
for (k in c("LNG", "DRSP", "ETG", "LNGP", "DRSPP")) {
  drug <- DRUGS[[k]]
  pos <- if (drug$kind == "POP") "mid" else "start"
  subs <- sample_population(NCAL, drug = drug, tag = "calib")
  ev <- escape_prob_par(cl, k, 0, pos, subs)
  rate <- 100 * mean(ev)
  ci <- 100 * wilson_ci(sum(ev), length(ev))
  tg <- TARGETS[[k]]
  ok <- (rate >= tg[1] - 0.5 && rate <= tg[2] + 0.5)
  cat(sprintf("  %-4s perfect-use ovulation = %5.1f%% (95%% CI %.1f-%.1f) | target %g-%g | %s | margin=%.2f emax=%.2f\n",
              k, rate, ci[1], ci[2], tg[1], tg[2], if (ok) "OK" else "ADJUST",
              drug$margin_prog, drug$emax_prog))
  pu[[k]] <- list(rate = rate, ci = ci, target = tg, margin = drug$margin_prog,
                  emax = drug$emax_prog)
}

# ---- Follicular-recovery / extended-HFI check (LNG, missed at pack start) ----
# Target (Elomaa 1998): omitting <=2-3 pills at start (extended HFI) -> follicular
# growth but escape stays low; escape rises as the pill-free interval extends.
cat("\n== LNG missed-at-start escape curve (extended HFI; Elomaa 1998) ==\n")
drug <- DRUGS$LNG
subs <- sample_population(NCAL, drug = drug, tag = "calib")
ns <- 0:7
lng_curve <- sapply(ns, function(n) 100 * mean(escape_prob_par(cl, "LNG", n, "start", subs)))
for (i in seq_along(ns))
  cat(sprintf("  missed=%d at start: escape = %4.1f%%\n", ns[i], lng_curve[i]))

# ---- Time-to-preovulatory follicle after suppression release (recovery rate) ----
# A typical subject: stop perfect-use LNG, measure days from last dose until the
# follicle reaches preovulatory size (proxy for the follicular-recovery timescale).
cat("\n== Follicular-recovery timescale (typical subject, LNG stop) ==\n")
sched <- build_pack_schedule(21, 7, n_packs = 3); taken <- apply_adherence(sched)
last_dose <- max(taken); t_end <- last_dose + 40
prof <- make_profile(drug, taken, t_end)
s <- hpo_simulate(prof, t_end, ec50_prog(drug), drug$emax_prog, drug$hill_prog,
                  ec50_ee(drug), drug$ee_emax, drug$ee_hill)
fov <- P_default()$F_ov
after <- s$t >= last_dose
t_pre <- s$t[after][which(s$Y[iF, after] >= fov)[1]]
cat(sprintf("  days from last dose to preovulatory follicle (F>=%.0f): %s\n",
            fov, if (length(t_pre) && !is.na(t_pre)) sprintf("%.1f", t_pre - last_dose) else "none within 40 d"))
cat(sprintf("  first ovulation after stop: %s d\n",
            if (length(s$ov_times)) sprintf("%.1f", s$ov_times[1] - last_dose) else "none"))

writeLines(toJSON(list(N = NCAL, perfect_use = pu,
                       lng_missed_start = setNames(as.list(round(lng_curve, 1)), paste0("n", ns)),
                       targets = TARGETS), pretty = TRUE, auto_unbox = TRUE),
           "results_R/calibration.json")
cat("\nWrote results_R/calibration.json\n")
