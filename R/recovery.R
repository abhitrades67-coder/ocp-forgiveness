# =============================================================================
# recovery.R -- Optimal recovery-action analysis (Protocol Section 5.2;
# Counterman-Lawley framing).
#
# A subject misses a single active pill, discovered `discover_h` later. We compare
# three recovery actions and report which minimises ovulation-escape probability:
#   take_now : take the missed pill immediately (it is `discover_h` late), resume
#   double   : skip until the next scheduled dose, then take two pills together
#   skip     : omit the missed pill entirely, resume on schedule
# Implemented exactly with the linear dosing engine: two doses at the same instant
# superpose to a double dose.
# =============================================================================

# Build the taken-dose schedule for one recovery action around a missed pill at
# 1-based index `mk` in `sched`.
recovery_taken <- function(sched, mk, action, discover_h = 12) {
  n <- length(sched)
  if (action == "skip") {
    return(apply_adherence(sched, missed_idx = mk))
  } else if (action == "take_now") {
    delay <- setNames(discover_h, as.character(mk))   # take the missed pill late
    return(apply_adherence(sched, delay_h = delay))
  } else if (action == "double") {
    taken <- apply_adherence(sched, missed_idx = mk)  # missed at its slot...
    if (mk + 1 <= n) taken <- sort(c(taken, sched[mk + 1]))  # ...taken with next dose
    return(taken)
  }
  stop("unknown action")
}

# Escape probability for a recovery action over a population (parallel).
recovery_escape <- function(cl, drug_key, position, action, subs, discover_h = 12) {
  res <- clusterApply(cl, subs, function(s, dk, pos, act, dh) {
    drug <- DRUGS[[dk]]
    reg <- if (drug$kind == "POP") "continuous" else "21/7"
    rg <- REGIMENS[[reg]]; active <- rg[1]; hfi <- rg[2]
    sched <- build_pack_schedule(active, hfi, n_packs = 4)
    base <- 2 * active
    off <- switch(pos, start = 0, mid = 7, late = 14, 0)
    mk <- base + off + 1L
    taken <- recovery_taken(sched, mk, act, dh)
    t_chal <- sched[base + off + 1L]; t_end <- t_chal + 30
    prof <- make_profile(drug, taken, t_end, cl_factor = s$cl_factor, vd_factor = s$vd_factor,
                         ee_cl_factor = s$ee_cl_factor)
    p <- P_default(); if (!is.null(s$pd$ec_scale)) {}
    ec_scale <- if (!is.null(s$pd$ec_scale)) s$pd$ec_scale else 1
    for (nm in setdiff(names(s$pd), "ec_scale")) p[[nm]] <- s$pd[[nm]]
    ec_p <- ec50_prog(drug) * ec_scale
    ec_e <- if (isTRUE(drug$has_ee)) ec50_ee(drug) * ec_scale else NA_real_
    sim <- hpo_simulate(prof, t_end, ec_p, drug$emax_prog, drug$hill_prog,
                        ec_e, drug$ee_emax, drug$ee_hill, p = p)
    ov <- sim$ov_times
    as.integer(any(ov >= (t_chal - 0.5) & ov <= (t_chal + 28)))
  }, dk = drug_key, pos = position, act = action, dh = discover_h)
  mean(unlist(res))
}

recovery_table <- function(cl, DRUGS, n_sub = 200, discover_h = 12) {
  rows <- list()
  combos <- list(c("LNG", "start"), c("LNG", "mid"), c("DRSP", "start"),
                 c("ETG", "mid"), c("LNGP", "mid"), c("DRSPP", "mid"))
  for (cb in combos) {
    dk <- cb[1]; pos <- cb[2]
    subs <- sample_population(n_sub, drug = DRUGS[[dk]], tag = "recovery")
    e <- sapply(c("take_now", "double", "skip"),
                function(a) 100 * recovery_escape(cl, dk, pos, a, subs, discover_h))
    best <- names(which.min(e))
    rows[[length(rows) + 1]] <- data.frame(
      drug = dk, position = pos, discover_h = discover_h,
      take_now_pct = round(e["take_now"], 1), double_pct = round(e["double"], 1),
      skip_pct = round(e["skip"], 1), recommended = best, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
