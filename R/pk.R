# =============================================================================
# pk.R -- Population PK layer and dosing engine (Protocol Section 4.2, Appendix A.1)
#
# One-compartment, first-order-absorption model. Plasma concentration is the
# analytic superposition of single-dose responses (each taken pill is an
# independent first-order absorption/elimination bolus). Exact for the linear
# one-compartment model and far faster than integrating the PK ODE alongside the
# PD system.
#
# Concentrations are in ng/mL when dose is in micrograms and V/F in litres.
# R port of src/pk.py (verified equivalent), used by the journal R pipeline.
# =============================================================================

# Nominal once-daily dosing times (days) for a packs-with-HFI regimen.
build_pack_schedule <- function(active_pills, hfi_days, n_packs,
                                regimen_start_day = 0, dose_hour = 8) {
  cycle_len <- active_pills + hfi_days
  times <- numeric(0)
  for (p in seq_len(n_packs) - 1L) {
    pack_start <- regimen_start_day + p * cycle_len
    times <- c(times, pack_start + (seq_len(active_pills) - 1L) + dose_hour / 24)
  }
  times
}

# Turn a nominal schedule into an actual intake schedule.
#   missed_idx : 1-based indices of OMITTED pills
#   delay_h    : named numeric (names = 1-based indices) of hours late
apply_adherence <- function(scheduled_times, missed_idx = NULL, delay_h = NULL) {
  keep <- setdiff(seq_along(scheduled_times), as.integer(missed_idx))
  taken <- scheduled_times[keep]
  if (!is.null(delay_h) && length(delay_h)) {
    idx <- as.integer(names(delay_h))
    for (j in seq_along(delay_h)) {
      pos <- which(keep == idx[j])
      if (length(pos)) taken[pos] <- taken[pos] + delay_h[j] / 24
    }
  }
  sort(taken)
}

# Plasma concentration on t_eval (days) from a vector of taken-dose times.
# Linear one-compartment, first-order absorption; superposition of doses.
conc_timecourse <- function(dose_times_days, t_eval_days, dose_amt, F, Vd, ke, Ka) {
  if (length(dose_times_days) == 0) return(numeric(length(t_eval_days)))
  ke_d <- ke * 24
  Ka_d <- Ka * 24
  if (abs(Ka_d - ke_d) < 1e-9) ke_d <- ke_d * (1 - 1e-6)   # flip-flop singularity
  coef <- (F * dose_amt * Ka_d) / (Vd * (Ka_d - ke_d))
  C <- numeric(length(t_eval_days))
  for (td in dose_times_days) {
    dt <- t_eval_days - td
    m <- dt >= 0
    if (any(m)) C[m] <- C[m] + coef * (exp(-ke_d * dt[m]) - exp(-Ka_d * dt[m]))
  }
  C
}

# Concentration profile (progestin and, for COCs, EE) on a fine grid, with
# linear-interpolation accessors the ODE can query cheaply at any time.
make_profile <- function(drug, dose_times_days, t_end_days, dt_h = 0.25,
                         cl_factor = 1.0, vd_factor = 1.0, ee_cl_factor = NULL) {
  t_grid <- seq(0, t_end_days + dt_h / 24, by = dt_h / 24)
  ke <- drug$ke_h * cl_factor / vd_factor
  Vd <- drug$Vd_L * vd_factor
  c_prog <- conc_timecourse(dose_times_days, t_grid, drug$dose_ug, drug$F, Vd, ke, drug$Ka_h)
  c_ee <- NULL
  if (isTRUE(drug$has_ee)) {
    # Ethinylestradiol carries its own CYP3A4-dependent fraction, so under
    # induction/inhibition its clearance moves independently of the progestin's.
    ee_cf <- if (is.null(ee_cl_factor)) cl_factor else ee_cl_factor
    ee_ke <- drug$ee_ke_h * ee_cf / vd_factor
    ee_Vd <- drug$ee_Vd_L * vd_factor
    ee_F <- if (!is.null(drug$ee_F)) drug$ee_F else drug$F   # EE-specific oral bioavailability
    c_ee <- conc_timecourse(dose_times_days, t_grid, drug$ee_dose_ug, ee_F, ee_Vd, ee_ke, drug$ee_Ka_h)
  }
  list(t = t_grid, c_prog = c_prog,
       c_ee = if (is.null(c_ee)) numeric(length(t_grid)) else c_ee)
}

prof_prog <- function(profile, t) approx(profile$t, profile$c_prog, t, rule = 2)$y
prof_ee   <- function(profile, t) approx(profile$t, profile$c_ee,   t, rule = 2)$y
