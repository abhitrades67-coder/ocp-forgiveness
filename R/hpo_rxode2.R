# =============================================================================
# hpo_rxode2.R -- rxode2 implementation of the HPO-axis ODE engine.
#
# WHY THIS EXISTS. The study protocol named R (mrgsolve / rxode2) as an intended
# toolchain. The production engine in R/hpo.R is a self-contained pure-R RK4
# integrator (deSolve's compiled lsoda path segfaulted under the repeated-call
# Monte-Carlo workload on this build). This module provides the SAME HPO model on
# rxode2's compiled-C ODE solver, as a named-toolchain alternative, and is
# cross-validated against the pure-R engine in R/validate_rxode2.R.
#
# rxode2 compiles the right-hand side to C ONCE and then solves via compiled code
# (it does NOT call back into R for the RHS), which is a different and more stable
# path than deSolve-with-an-R-RHS; it runs cleanly here.
#
# DESIGN. rxode2 integrates ODEs but does not natively apply a STATE-triggered
# discrete reset mid-integration. Ovulation here fires when 'surge permission'
# (1-sLH)*hill(E2) crosses surge_crit upward with a preovulatory follicle, then
# luteinises the follicle. We therefore use SEGMENT-AND-RESTART: integrate with
# rxode2 between events, detect the first upward crossing (respecting the
# refractory window), apply the physiological reset, and restart integration from
# the reset state -- exactly the crossing/refractory/reset rule of R/hpo.R, with
# rxode2's compiled solver replacing the inner RK4 loop. Drug forcing (sLH, sFSH)
# is supplied as time-varying covariates (linear interpolation), matching the
# pure-R forcing.
#
# The function hpo_simulate_rx() has the SAME signature and return shape
# (list(t, Y[6 x n], ov_times)) as hpo_simulate(), so it is drop-in for checks.
# Requires R/hpo.R (P_default, initial_state, hill, STATE_NAMES) to be sourced.
# =============================================================================

suppressMessages(library(rxode2))

# Compile the HPO RHS once; cache the compiled model object.
.hpo_rx_cache <- new.env(parent = emptyenv())
.hpo_rx_model <- function() {
  if (!is.null(.hpo_rx_cache$mod)) return(.hpo_rx_cache$mod)
  mod <- rxode2({
    # state clamps via boolean*value (C evaluates (x>0) as 0/1), matching the
    # max(.,0) guards in the pure-R deriv; LH is used unclamped as in hpo.R.
    FSHp = FSH * (FSH > 0)
    Folp = Fol * (Fol > 0)
    E2p  = E2  * (E2  > 0)
    LUTp = LUT * (LUT > 0)
    P4p  = P4  * (P4  > 0)
    fb_fsh    = 1 / (1 + (E2p / E2_fb)^nfb + (P4p / P4_fb)^nfb)
    fb_gnrh   = 1 / (1 + (P4p / P4_gn)^nfb)
    FSH_drive = (FSHp - FSH_eff) * ((FSHp - FSH_eff) > 0)
    surge     = ksurge * (1 - sLH) * (E2p^n_surge / (E2_surge^n_surge + E2p^n_surge))
    d/dt(FSH) = kFSH_syn * (1 - sFSH) * fb_fsh - kFSH_el * FSHp
    d/dt(LH)  = kLH_syn  * (1 - sLH)  * fb_gnrh + surge - kLH_el * LH
    d/dt(Fol) = kgrow * FSH_drive * (1 - Folp / Fmax) - k_atr * Folp
    d/dt(E2)  = kE2 * Folp - kE2_el * E2p
    d/dt(LUT) = -klut_reg * LUTp
    d/dt(P4)  = kP4 * LUTp - kP4_el * P4p
  })
  .hpo_rx_cache$mod <- mod
  mod
}

# rxode2 twin of hpo_simulate(). Same arguments (need_Y accepted for signature
# compatibility; trajectory is always assembled here).
hpo_simulate_rx <- function(profile, t_end_days, ec50_prog, emax_prog, hill_prog,
                            ec50_ee, emax_ee, hill_ee, p = P_default(),
                            y0 = NULL, dt_out_h = 3.0, h_int = 0.02, need_Y = TRUE) {
  mod <- .hpo_rx_model()
  if (is.null(y0)) y0 <- initial_state()

  # fine grid (crossing-detection resolution = pure-R h_int)
  nstep <- max(1L, as.integer(ceiling(t_end_days / h_int)))
  tf <- seq(0, t_end_days, length.out = nstep + 1L)
  cprog <- approx(profile$t, profile$c_prog, tf, rule = 2)$y
  sLHf <- emax_prog * hill(cprog, ec50_prog, hill_prog)
  if (!is.na(ec50_ee)) {
    cee <- approx(profile$t, profile$c_ee, tf, rule = 2)$y
    sFSHf <- emax_ee * hill(cee, ec50_ee, hill_ee)
  } else sFSHf <- numeric(length(tf))

  pars <- c(E2_fb = p$E2_fb, P4_fb = p$P4_fb, P4_gn = p$P4_gn, nfb = p$n_fb,
            kFSH_syn = p$kFSH_syn, kFSH_el = p$kFSH_el, FSH_eff = p$FSH_eff,
            kgrow = p$kgrow, Fmax = p$Fmax, k_atr = p$k_atr,
            kE2 = p$kE2, kE2_el = p$kE2_el, kLH_syn = p$kLH_syn, kLH_el = p$kLH_el,
            ksurge = p$ksurge, E2_surge = p$E2_surge, n_surge = p$n_surge,
            klut_reg = p$klut_reg, kP4 = p$kP4, kP4_el = p$kP4_el)

  E2s <- p$E2_surge; ns <- p$n_surge; sc <- p$surge_crit
  scol <- c("FSH", "LH", "Fol", "E2", "LUT", "P4")

  Yfull <- matrix(NA_real_, nstep + 1L, 6); Yfull[1, ] <- as.numeric(y0)
  y <- as.numeric(y0); ov <- numeric(0); last_event <- -Inf
  cursor <- 1L

  repeat {
    idx <- cursor:(nstep + 1L)
    if (length(idx) < 2L) break
    tt <- tf[idx]; t0 <- tt[1]
    dat <- data.frame(id = 1L, time = tt - t0, sLH = sLHf[idx], sFSH = sFSHf[idx])
    inits <- c(FSH = y[1], LH = y[2], Fol = y[3], E2 = y[4], LUT = y[5], P4 = y[6])
    sol <- rxSolve(mod, params = pars, events = dat, inits = inits,
                   covsInterpolation = "linear", returnType = "matrix",
                   cores = 1L, atol = 1e-8, rtol = 1e-8)
    Yseg <- sol[, scol, drop = FALSE]
    Yfull[idx, ] <- Yseg

    E2seg <- pmax(Yseg[, 4], 0)
    sp <- (1 - sLHf[idx]) * (E2seg^ns / (E2s^ns + E2seg^ns)) - sc
    cross <- which(sp[-length(sp)] <= 0 & sp[-1] > 0)   # crossing between k and k+1

    fired <- FALSE
    for (k in cross) {
      tc <- tt[k + 1L]
      if ((tc - last_event) >= p$refractory_d) {
        yk <- as.numeric(Yseg[k + 1L, ]); last_event <- tc
        if (yk[3] >= p$F_ov) {                       # genuine ovulation
          ov <- c(ov, tc); yk[3] <- yk[3] * 0.10; yk[5] <- p$LUT0; yk[4] <- yk[4] * 0.5
        } else {                                     # marginal: luteinise, uncounted
          yk[3] <- yk[3] * 0.5; yk[5] <- 0.5 * p$LUT0; yk[4] <- yk[4] * 0.5
        }
        y <- yk; cursor <- idx[k + 1L]; Yfull[cursor, ] <- y
        fired <- TRUE; break
      }
    }
    if (!fired) break
  }

  times <- seq(0, t_end_days, by = dt_out_h / 24)
  Y <- vapply(1:6, function(j) approx(tf, Yfull[, j], times, rule = 2)$y, numeric(length(times)))
  Y <- t(Y); rownames(Y) <- STATE_NAMES
  list(t = times, Y = Y, ov_times = ov)
}
