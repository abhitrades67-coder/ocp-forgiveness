# =============================================================================
# translate.R -- Clinical translation: cervical-mucus arm + ovulation-rate ->
# Pearl-Index mapping with uncertainty propagation.
# (Protocol Sections 4.4, 5.7.)
#
# Two components:
#
#  (1) CERVICAL-MUCUS ARM. POP (and, secondarily, COC) protection is not
#      ovulation-only: progestin thickens cervical mucus, impeding sperm. We model
#      a barrier competence B(t) = mucus_emax * hill(C_prog, EC50_mucus, n_mucus),
#      summarised over the fertile window. Because EC50_mucus is set so the
#      barrier wanes within ~1 day of a missed/late dose, this reproduces the
#      known POP timing windows (LNG-POP ~3 h, desogestrel ~12 h) and gives
#      traditional LNG-POP a DEFINED forgiveness on the failure endpoint even
#      though its ovulation arm is ~50% (resolving the prior 'undefined' artifact).
#
#  (2) OR -> PEARL-INDEX MAPPING. Ovulation escape only risks pregnancy
#      if sperm reach the egg. Per-cycle pregnancy risk =
#         P(ovulation escape) * fecundability * (1 - barrier competence),
#      annualised over 13 cycles. The published OR-PI correlation is only
#      r ~ 0.52, so OR explains ~27% of PI variance; we propagate this by
#      Monte-Carlo sampling fecundability AND a residual reflecting (1 - r^2),
#      yielding a CREDIBLE INTERVAL on the approximate Pearl Index rather than a
#      false point estimate. The cervical-barrier ceilings (drugs.R, mucus_emax)
#      are calibrated so the PERFECT-USE Pearl Indices reproduce published values
#      (COC ~0.3, desogestrel POP ~0.3, LNG-POP ~1.1 per 100 woman-years); the
#      missed-pill escalation and absolute typical-use values remain illustrative.
# =============================================================================

# Mean cervical-barrier competence over [t0, t1] for a (typical-subject) profile.
mucus_competence <- function(drug, profile, t0, t1) {
  ec <- ec50_mucus(drug)
  if (is.na(ec) || drug$mucus_emax <= 0) return(0.0)
  m <- profile$t >= t0 & profile$t <= t1
  if (!any(m)) return(0.0)
  B <- drug$mucus_emax * hill(profile$c_prog[m], ec, drug$mucus_hill)
  mean(B)
}

# Barrier competence for a (drug, position, n_missed) cell, typical subject.
cell_mucus <- function(drug, n_missed, position, scen_kwargs = list()) {
  args <- c(list(drug = drug, n_missed = n_missed, position = position,
                 return_trace = TRUE), scen_kwargs)
  r <- do.call(run_scenario, args)
  mucus_competence(drug, r$profile, r$t_challenge, r$t_challenge + 14)
}

OR_PI <- list(
  fecundability = 0.20,   # per-cycle conception prob given an unprotected ovulation
  fec_cv = 0.35,          # between-woman/illustrative spread
  cycles_per_year = 13,
  r_or_pi = 0.52,         # published OR-PI correlation (meta-analysis)
  resid_sd_log = NA       # set from r in pearl_index()
)

# Map an ovulation-escape probability (+ barrier competence) to an approximate
# annual Pearl Index with a propagated credible interval.
pearl_index <- function(escape_prob, barrier = 0.0, B = 4000, seed = 7) {
  set.seed(seed)
  # residual multiplicative noise reflecting that OR explains only r^2 of PI var
  r <- OR_PI$r_or_pi
  resid_sd <- sqrt(max(1e-6, 1 - r^2)) * 0.6      # illustrative scale on log-PI
  fec <- OR_PI$fecundability * exp(rnorm(B, 0, OR_PI$fec_cv))
  per_cycle <- pmin(0.95, escape_prob * fec * (1 - barrier))
  pi_annual <- 100 * (1 - (1 - per_cycle)^OR_PI$cycles_per_year)
  pi_annual <- pi_annual * exp(rnorm(B, 0, resid_sd))
  c(pi = stats::median(pi_annual),
    lo = unname(quantile(pi_annual, 0.05)),
    hi = unname(quantile(pi_annual, 0.95)))
}

# Build a Pearl-Index translation table from the main escape grid + mucus arm.
# grid_df: data.frame with columns drug, position, n_missed, escape_prob (main).
pearl_table <- function(DRUGS, grid_df) {
  rows <- list()
  for (i in seq_len(nrow(grid_df))) {
    g <- grid_df[i, ]
    drug <- DRUGS[[g$drug]]
    reg <- if (drug$kind == "POP") "continuous" else "21/7"
    barrier <- cell_mucus(drug, g$n_missed, g$position,
                          scen_kwargs = list(regimen = reg))
    pe <- pearl_index(g$escape_prob, barrier)
    rows[[i]] <- data.frame(drug = g$drug, position = g$position, n_missed = g$n_missed,
                            escape_prob = round(g$escape_prob, 3),
                            barrier_competence = round(barrier, 3),
                            pearl_index = round(pe["pi"], 2),
                            pearl_lo = round(pe["lo"], 2), pearl_hi = round(pe["hi"], 2),
                            stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
