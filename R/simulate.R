# =============================================================================
# simulate.R -- Integrated PK-PD simulation for a single virtual subject
# (Protocol Section 3, Figure 1). R port of src/simulate.py.
#
# A scenario = a regimen (pack structure) + a missed/delayed-dose pattern applied
# in a 'challenge' pack after a perfect-use lead-in that establishes steady-state
# suppression. Escape = an ovulation from the start of the challenge window
# onward (within one follow-up cycle).
# =============================================================================

REGIMENS <- list("21/7" = c(21, 7), "24/4" = c(24, 4),
                 "26/2" = c(26, 2), "continuous" = c(28, 0))

regimen_for <- function(drug) if (drug$kind == "POP") "continuous" else "21/7"

run_scenario <- function(drug, n_missed = 0, position = "mid", challenge_pack = 2,
                         n_packs = 4, regimen = NULL, delay_h = NULL, delay_idx = NULL,
                         cl_factor = 1.0, vd_factor = 1.0, ee_cl_factor = NULL,
                         pd = NULL, scattered = FALSE,
                         follow_days = 28.0, return_trace = FALSE) {
  if (is.null(regimen)) regimen <- regimen_for(drug)
  rg <- REGIMENS[[regimen]]; active <- rg[1]; hfi <- rg[2]
  sched <- build_pack_schedule(active, hfi, n_packs = n_packs)

  # 1-based index of the first pill of the challenge window (R is 1-based;
  # base = pills before the challenge pack, off = within-pack offset).
  base <- challenge_pack * active
  off <- switch(position, start = 0, mid = 7, late = 14, 0)
  first_chal <- base + off + 1L
  if (scattered) {
    miss <- first_chal + 2 * (seq_len(n_missed) - 1L)
  } else {
    miss <- if (n_missed > 0) first_chal + (seq_len(n_missed) - 1L) else integer(0)
  }
  miss <- miss[miss >= 1 & miss <= length(sched)]

  delays <- NULL
  if (!is.null(delay_h)) {
    if (is.null(delay_idx)) {
      k <- if (n_missed > 0) n_missed else 7
      delay_idx <- first_chal + (seq_len(k) - 1L)
    }
    delay_idx <- delay_idx[delay_idx >= 1 & delay_idx <= length(sched)]
    delays <- setNames(rep(delay_h, length(delay_idx)), as.character(delay_idx))
  }

  taken <- apply_adherence(sched, missed_idx = miss, delay_h = delays)

  idx_chal <- min(first_chal, length(sched))
  t_challenge <- sched[idx_chal]
  t_end <- as.numeric(t_challenge + follow_days + 2)

  profile <- make_profile(drug, taken, t_end, cl_factor = cl_factor, vd_factor = vd_factor,
                          ee_cl_factor = ee_cl_factor)

  p <- P_default()
  ec_scale <- 1.0
  if (!is.null(pd)) {
    if (!is.null(pd$ec_scale)) ec_scale <- pd$ec_scale
    for (nm in setdiff(names(pd), "ec_scale")) p[[nm]] <- pd[[nm]]
  }
  ec_p <- ec50_prog(drug) * ec_scale
  ec_e <- if (isTRUE(drug$has_ee)) ec50_ee(drug) * ec_scale else NA_real_

  sim <- hpo_simulate(profile, t_end, ec_p, drug$emax_prog, drug$hill_prog,
                      ec_e, drug$ee_emax, drug$ee_hill, p = p, need_Y = return_trace)
  ov <- sim$ov_times
  escape_ov <- ov[ov >= (t_challenge - 0.5) & ov <= (t_challenge + follow_days)]
  res <- list(escaped = length(escape_ov) > 0, ov_times = ov, escape_ov = escape_ov,
              t_challenge = t_challenge, ec_prog = ec_p, ec_ee = ec_e)
  if (return_trace) { res$t <- sim$t; res$Y <- sim$Y; res$profile <- profile }
  res
}
