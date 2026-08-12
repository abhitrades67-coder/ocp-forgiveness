# =============================================================================
# hpo.R -- QSP HPO-axis (pharmacodynamic) core (Protocol Section 4.3, Appendix A.3)
#
# Reduced ODE model of the hypothalamic-pituitary-ovarian axis: FSH, LH, a
# dominant-follicle pathway, estradiol (E2), corpus luteum and progesterone (P4),
# with endocrine feedback so the untreated model self-sustains ~monthly ovulatory
# cycles. Contraceptive effects enter as concentration-driven sigmoid-Emax
# suppression: progestin blocks the E2-triggered LH surge (ovulation block) and
# tonic LH; ethinylestradiol (COC) suppresses FSH (follicular recruitment).
#
# Ovulation is a discrete EVENT: it fires when 'surge permission'
# (1 - sLH) * hill(E2) crosses surge_crit upward; if a preovulatory follicle is
# present (F >= F_ov) it is recorded as an ovulation and the follicle luteinises
# (corpus luteum forms, P4 rises, luteal suppression), else it luteinises without
# being counted. A refractory window prevents the post-ovulation downward
# crossing from re-triggering.
#
# NOTE ON THE LH STATE: the ovulation event is
# governed by the surge-PERMISSION signal Perm(t) = (1 - sLH)*hill(E2) (Eq. 15),
# i.e. a preovulatory follicle is released when estradiol crosses the surge
# set-point and the progestin is not blocking. The LH ODE (Eq. 10) integrates the
# corresponding gonadotropin surge as a DIAGNOSTIC readout of that same permission
# signal; it is not itself the trigger, so LH inaccuracies (it carries a fast,
# bounded surge source) cannot bias the ovulation decision. This is the standard
# reduction for a suppression-focused model and is stated as such in the paper.
#
# Event handling uses a rootfun-style upward-
# crossing test + a physiological luteinisation reset with a refractory flag,
# replacing the prior code's non-physical 'F*=0.5, E2*=0.5 to clear the crossing'
# device.
#
# Parameters: see P_default(). The ovarian-dynamics parameters kgrow / F_ov /
# surge_crit / E2_surge are set so the UNTREATED oscillator reproduces (i) a
# ~28-day cycle length and luteal-phase duration and (ii) a physiological
# follicular-recovery timescale, checked against natural-cycle reference values
# (Anckaert 2021; R/validate_cycle.R). They are NOT fitted to any missed-pill
# escape data, so the Elomaa 1998 missed-at-start escape curve remains an
# out-of-sample (independent) check. The progestin ovulation-suppression
# potencies (margin_prog) are calibrated separately, to published perfect-use
# ovulation rates (R/calibrate.R).
# =============================================================================

# NOTE: deSolve is intentionally NOT used -- its compiled lsoda/lsodar path
# segfaults under repeated calls on this R build. The integrator below is pure R.

# state indices (1-based)
iFSH <- 1; iLH <- 2; iF <- 3; iE2 <- 4; iLUT <- 5; iP4 <- 6
STATE_NAMES <- c("FSH", "LH", "Follicle", "E2", "CorpusLuteum", "P4")

P_default <- function() list(
  # FSH (basal = kFSH_syn/kFSH_el = 2.0); n_fb = endocrine-feedback Hill power
  kFSH_syn = 12.0, kFSH_el = 6.0, E2_fb = 120.0, P4_fb = 8.0, n_fb = 2.0,
  # Follicle (FSH-driven growth, slow atresia, logistic cap)
  kgrow = 0.85, Fmax = 30.0, k_atr = 0.025, FSH_eff = 0.40,
  # E2 (~4 * follicle at quasi-steady state)
  kE2 = 8.0, kE2_el = 2.0,
  # LH (diagnostic surge trace): tonic + E2-triggered, progestin-suppressible
  kLH_syn = 4.0, kLH_el = 6.0, P4_gn = 6.0, ksurge = 240.0,
  # ovulation event ("surge permission" = (1-sLH)*hill(E2) crosses surge_crit)
  E2_surge = 46.0, n_surge = 8.0, surge_crit = 0.50, F_ov = 11.0,
  # luteal phase
  LUT0 = 10.0, klut_reg = 0.07, kP4 = 3.0, kP4_el = 2.0,
  # ovulation-event refractory (days): blocks re-trigger on the post-surge
  # downward crossing; << inter-ovulation interval so real ovulations are kept
  refractory_d = 6.0
)

hill <- function(x, K, n) { x <- pmax(x, 0); x^n / (K^n + x^n) }

supp_prog <- function(Cpg, ec50, emax, hillc) emax * hill(Cpg, ec50, hillc)
supp_ee   <- function(Cee, ec50, emax, hillc) if (is.na(ec50)) 0.0 else emax * hill(Cee, ec50, hillc)

initial_state <- function() c(FSH = 1.5, LH = 2.0, Follicle = 0.5, E2 = 5.0,
                              CorpusLuteum = 0.0, P4 = 0.0)

# Surge-permission signal that the ovulation event roots on.
surge_permission <- function(t, y, profile, ec50_prog, emax_prog, hill_prog, p) {
  E2 <- max(y[iE2], 0)
  sLH <- supp_prog(prof_prog(profile, t), ec50_prog, emax_prog, hill_prog)
  (1 - sLH) * hill(E2, p$E2_surge, p$n_surge)
}

# Integrate the HPO axis under a concentration profile, with ovulation events.
# Returns list(t, Y [6 x n on the requested grid], ov_times).
#
# Implementation: a self-contained pure-R RK4
# integrator with the ovulation event detected as an upward crossing of the
# surge-permission threshold and applied in-loop (physiological luteinisation +
# refractory window). This replaces deSolve, whose compiled lsoda/lsodar path
# segfaults under repeated calls on this R build -- unacceptable for a 10^5-sim
# sweep. The HPO axis is non-stiff in the states that drive ovulation (rates
# <= 12/day); LH carries only a fast diagnostic surge term that is a bounded
# source, well inside RK4 stability at the chosen step. Drug forcing (sLH, sFSH)
# is pre-evaluated on the integration half-grid so the inner loop is allocation-
# light. need_Y=FALSE skips trajectory assembly for the escape-only sweep.
hpo_simulate <- function(profile, t_end_days, ec50_prog, emax_prog, hill_prog,
                         ec50_ee, emax_ee, hill_ee, p = P_default(),
                         y0 = NULL, dt_out_h = 3.0, h_int = 0.02, need_Y = TRUE) {
  if (is.null(y0)) y0 <- initial_state()
  nfb <- p$n_fb
  E2_fb <- p$E2_fb; P4_fb <- p$P4_fb; P4_gn <- p$P4_gn
  kFSH_syn <- p$kFSH_syn; kFSH_el <- p$kFSH_el; FSH_eff <- p$FSH_eff
  kgrow <- p$kgrow; Fmax <- p$Fmax; k_atr <- p$k_atr
  kE2 <- p$kE2; kE2_el <- p$kE2_el
  kLH_syn <- p$kLH_syn; kLH_el <- p$kLH_el; ksurge <- p$ksurge
  E2_surge <- p$E2_surge; n_surge <- p$n_surge; surge_crit <- p$surge_crit
  F_ov <- p$F_ov; LUT0 <- p$LUT0; klut_reg <- p$klut_reg; kP4 <- p$kP4; kP4_el <- p$kP4_el

  nstep <- max(1L, as.integer(ceiling(t_end_days / h_int)))
  h <- t_end_days / nstep
  # forcing on the half-grid (2*nstep+1 nodes): sLH (progestin) and sFSH (EE)
  tf <- seq(0, t_end_days, by = h / 2)
  cprog <- approx(profile$t, profile$c_prog, tf, rule = 2)$y
  sLHf <- emax_prog * (function(x){x<-pmax(x,0); x^hill_prog/(ec50_prog^hill_prog+x^hill_prog)})(cprog)
  if (!is.na(ec50_ee)) {
    cee <- approx(profile$t, profile$c_ee, tf, rule = 2)$y
    sFSHf <- emax_ee * (function(x){x<-pmax(x,0); x^hill_ee/(ec50_ee^hill_ee+x^hill_ee)})(cee)
  } else sFSHf <- numeric(length(tf))

  # derivative given state vector y (length 6) and forcing (sLH, sFSH) at a node
  deriv <- function(y, sLH, sFSH) {
    FSH <- if (y[1] > 0) y[1] else 0; LH <- y[2]
    Fol <- if (y[3] > 0) y[3] else 0; E2 <- if (y[4] > 0) y[4] else 0
    LUT <- if (y[5] > 0) y[5] else 0; P4 <- if (y[6] > 0) y[6] else 0
    fb_fsh  <- 1 / (1 + (E2 / E2_fb)^nfb + (P4 / P4_fb)^nfb)
    fb_gnrh <- 1 / (1 + (P4 / P4_gn)^nfb)
    FSH_drive <- if (FSH - FSH_eff > 0) FSH - FSH_eff else 0
    surge <- ksurge * (1 - sLH) * (E2^n_surge / (E2_surge^n_surge + E2^n_surge))
    c(kFSH_syn * (1 - sFSH) * fb_fsh - kFSH_el * FSH,
      kLH_syn * (1 - sLH) * fb_gnrh + surge - kLH_el * LH,
      kgrow * FSH_drive * (1 - Fol / Fmax) - k_atr * Fol,
      kE2 * Fol - kE2_el * E2,
      -klut_reg * LUT,
      kP4 * LUT - kP4_el * P4)
  }

  y <- as.numeric(y0)
  ov <- numeric(0); last_event <- -Inf
  if (need_Y) { Yint <- matrix(0, nstep + 1L, 6); Yint[1, ] <- y }
  sp_prev <- (1 - sLHf[1]) * (max(y[4], 0)^n_surge / (E2_surge^n_surge + max(y[4], 0)^n_surge)) - surge_crit

  for (i in seq_len(nstep)) {
    a <- 2L * i - 1L; b <- 2L * i; cc <- 2L * i + 1L     # half-grid indices
    k1 <- deriv(y,               sLHf[a],  sFSHf[a])
    k2 <- deriv(y + (h/2) * k1,  sLHf[b],  sFSHf[b])
    k3 <- deriv(y + (h/2) * k2,  sLHf[b],  sFSHf[b])
    k4 <- deriv(y + h * k3,      sLHf[cc], sFSHf[cc])
    y <- y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    y[y < 0] <- 0
    ti <- i * h
    E2i <- y[4]
    sp <- (1 - sLHf[cc]) * (E2i^n_surge / (E2_surge^n_surge + E2i^n_surge)) - surge_crit
    if (sp_prev <= 0 && sp > 0 && (ti - last_event) >= p$refractory_d) {
      last_event <- ti
      if (y[3] >= F_ov) {                 # genuine ovulation
        ov <- c(ov, ti)
        y[3] <- y[3] * 0.10; y[5] <- LUT0; y[4] <- y[4] * 0.5
      } else {                            # marginal: luteinise, not counted
        y[3] <- y[3] * 0.5; y[5] <- 0.5 * LUT0; y[4] <- y[4] * 0.5
      }
      sp <- (1 - sLHf[cc]) * (y[4]^n_surge / (E2_surge^n_surge + y[4]^n_surge)) - surge_crit
    }
    sp_prev <- sp
    if (need_Y) Yint[i + 1L, ] <- y
  }

  times <- seq(0, t_end_days, by = dt_out_h / 24)
  if (need_Y) {
    tint <- seq(0, t_end_days, length.out = nstep + 1L)
    Y <- vapply(1:6, function(j) approx(tint, Yint[, j], times, rule = 2)$y, numeric(length(times)))
    Y <- t(Y); rownames(Y) <- STATE_NAMES
  } else {
    Y <- NULL
  }
  list(t = times, Y = Y, ov_times = ov)
}
