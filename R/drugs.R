# =============================================================================
# drugs.R -- Drug definitions (Protocol Section 4.1) for the journal R pipeline.
#
# Representative agents span the PK / metabolic range so conclusions generalise:
#   LNG  (COC)  long t1/2, accumulates, low CYP3A4 dependence  -> high-forgiveness anchor
#   DRSP (COC)  substantial CYP3A4 metabolism                  -> DDI-sensitive anchor
#   EE   (COC)  estrogen component, suppresses FSH             -> standard estrogen
#   ETG  (POP)  desogestrel/etonogestrel, ovulation-inhibiting -> POP forgiveness contrast
#   LNGP (POP)  traditional LNG-POP, narrow window             -> low-forgiveness extreme
#   DRSPP(POP)  drospirenone-only (Slynd), long t1/2, 24h window-> most-forgiving POP
#
# PK: one-compartment, first-order absorption. Disposition (t1/2, F, V/F) sits
# within published ranges (DRSP/EE Blode 2000, Reif 2013; LNG/EE Carol 1992) and
# Ka reproduces t_max 1.5-2.0 h. EC50 for ovulation suppression is anchored to
# each agent's perfect-use steady-state trough via a 'forgiveness margin'
# (trough / EC50) -- the cushion-above-threshold of Protocol Figure 3.
#
# CALIBRATION NOTES:
#  * margin_prog anchored so that perfect-use COC ovulation is
#    ~0-1% (vs 3-5%), matching real perfect-use suppression; values fixed by the
#    calibration harness R/calibrate.R against documented targets.
#  * Cervical-mucus arm, calibrated: a concentration-driven
#    cervical-barrier competence so POP protection is no longer ovulation-arm-only.
#    Used by the OR->Pearl-Index translation (R/translate.R). The barrier ceilings
#    (mucus_emax) are now set so perfect-use Pearl Indices reproduce published
#    values (COC ~0.3, desogestrel POP ~0.3, LNG-POP ~1.1 per 100 woman-years).
#    The barrier represents the COMPOSITE per-escape-cycle protection -- progestin-
#    thickened cervical mucus PLUS the reduced fertility of a hormonally-induced
#    breakthrough ovulation -- which is why its perfect-use value is high (0.93-0.99);
#    it wanes as doses are missed and mucus hostility recovers within the known POP
#    windows (LNG-POP ~3 h, desogestrel ~12 h grace).
#
# FURTHER CALIBRATION NOTES:
#  * ee_F: ethinylestradiol now carries its own oral bioavailability (~0.43)
#    instead of reusing the progestin's F. Outcome-neutral -- the EE EC50 is
#    anchored to the EE steady-state trough (both trough and profile scale with
#    ee_F, so the concentration/EC50 ratio, and hence FSH suppression, are
#    invariant) -- but it removes a latent parameter error (a reader would have
#    seen EE F = 0.95 for the LNG COC).
#  * grace_h REMOVED: it was defined and documented as if mechanistic but was
#    never consumed by the simulation. The missed-dose "grace windows" are wholly
#    EMERGENT from washout of the trough below EC50; they are not an input. The
#    nominal per-agent windows (LNG-POP ~3 h, desogestrel ~12 h, drospirenone-only
#    ~24 h) are reported in the text as the clinical comparators the emergent
#    behaviour is checked against, not as model parameters.
# =============================================================================

# Analytic steady-state trough, one-compartment first-order absorption, once-daily.
ss_trough <- function(dose_ug, F, Vd_L, ke, Ka, tau_h = 24) {
  if (abs(Ka - ke) < 1e-9) ke <- ke * (1 - 1e-6)
  amp <- (F * dose_ug * Ka) / (Vd_L * (Ka - ke))
  amp * (exp(-ke * tau_h) / (1 - exp(-ke * tau_h)) -
         exp(-Ka * tau_h) / (1 - exp(-Ka * tau_h)))
}

make_drug <- function(name, label, kind, dose_ug, F, Vd_L, t_half_h, Ka_h,
                      fm_cyp3a4, margin_prog, emax_prog, hill_prog = 2.0,
                      has_ee = FALSE, ee_dose_ug = 0, ee_F = 0.43, ee_Vd_L = 350,
                      ee_t_half_h = 18, ee_Ka_h = 2.5, ee_fm_cyp3a4 = 0.5,
                      ee_margin = 3.0, ee_emax = 0.85, ee_hill = 2.0,
                      # cervical-mucus arm (illustrative; see header)
                      mucus_margin = NA_real_, mucus_emax = 0.0,
                      mucus_hill = 3.0) {
  d <- list(name = name, label = label, kind = kind,
            dose_ug = dose_ug, F = F, Vd_L = Vd_L, t_half_h = t_half_h, Ka_h = Ka_h,
            fm_cyp3a4 = fm_cyp3a4, margin_prog = margin_prog, emax_prog = emax_prog,
            hill_prog = hill_prog,
            has_ee = has_ee, ee_dose_ug = ee_dose_ug, ee_F = ee_F, ee_Vd_L = ee_Vd_L,
            ee_t_half_h = ee_t_half_h, ee_Ka_h = ee_Ka_h, ee_fm_cyp3a4 = ee_fm_cyp3a4,
            ee_margin = ee_margin, ee_emax = ee_emax, ee_hill = ee_hill,
            mucus_margin = mucus_margin, mucus_emax = mucus_emax, mucus_hill = mucus_hill)
  d$ke_h <- log(2) / d$t_half_h
  d$CL_Lh <- d$ke_h * d$Vd_L
  d$ee_ke_h <- log(2) / d$ee_t_half_h
  d
}

ec50_prog <- function(drug) {
  tr <- ss_trough(drug$dose_ug, drug$F, drug$Vd_L, drug$ke_h, drug$Ka_h)
  tr / drug$margin_prog
}
ec50_ee <- function(drug) {
  if (!isTRUE(drug$has_ee)) return(NA_real_)
  tr <- ss_trough(drug$ee_dose_ug, drug$ee_F, drug$ee_Vd_L, drug$ee_ke_h, drug$ee_Ka_h)
  tr / drug$ee_margin
}
# Cervical-mucus EC50 (NA if the agent has no modelled mucus arm).
ec50_mucus <- function(drug) {
  if (is.na(drug$mucus_margin) || drug$mucus_emax <= 0) return(NA_real_)
  tr <- ss_trough(drug$dose_ug, drug$F, drug$Vd_L, drug$ke_h, drug$Ka_h)
  tr / drug$mucus_margin
}

# Multiplicative factor on clearance for CYP3A4 status (Section 4.7).
# Induction (rifampicin / enzyme-inducing antiepileptics) ~4x the CYP3A4-mediated
# fraction; strong inhibition ~0.3x.
cyp3a4_clearance_factor <- function(drug_fm, status) {
  fold <- c(normal = 1.0, induced = 4.0, inhibited = 0.3)[[status]]
  drug_fm * fold + (1 - drug_fm)
}

# ---------------------------------------------------------------------------
# Representative agents. margin_prog / emax_prog / mucus_* are LOCKED by
# R/calibrate.R to documented qualification targets (see results_R/calibration.json).
# Values below are the calibrated set; re-running calibrate.R reproduces them.
# ---------------------------------------------------------------------------
build_drugs <- function() {
  list(
    LNG = make_drug(
      name = "LNG", label = "Levonorgestrel COC (0.15 mg)", kind = "COC",
      dose_ug = 150, F = 0.95, Vd_L = 126, t_half_h = 26, Ka_h = 2.5,
      # fm_cyp3a4: levonorgestrel is substantially CYP3A4-cleared (published
      # PBPK fraction 0.33-0.47; the SmPC names CYP3A4 as the main route).
      # 0.44 reproduces the observed rifampicin AUC ratio of 0.43, since
      # AUC_ratio = 1/(1 + 3*fm) under the 4-fold induction factor of Eq. (19).
      fm_cyp3a4 = 0.44, margin_prog = 9.5, emax_prog = 0.98,
      has_ee = TRUE, ee_dose_ug = 30, ee_F = 0.43, ee_t_half_h = 18, ee_fm_cyp3a4 = 0.5,
      ee_margin = 4.0, ee_emax = 0.90,
      mucus_margin = 2.5, mucus_emax = 0.98, mucus_hill = 3.0),
    DRSP = make_drug(
      name = "DRSP", label = "Drospirenone COC (3 mg)", kind = "COC",
      dose_ug = 3000, F = 0.76, Vd_L = 280, t_half_h = 31, Ka_h = 2.5,
      fm_cyp3a4 = 0.50, margin_prog = 7.5, emax_prog = 0.98,
      has_ee = TRUE, ee_dose_ug = 30, ee_F = 0.43, ee_t_half_h = 18, ee_fm_cyp3a4 = 0.5,
      ee_margin = 4.0, ee_emax = 0.90,
      mucus_margin = 2.6, mucus_emax = 0.986, mucus_hill = 3.0),
    ETG = make_drug(
      # Desogestrel/etonogestrel: the reliably ovulation-inhibiting POP. A
      # comfortable perfect-use margin keeps even high-clearance subjects
      # suppressed at steady state; the short missed-pill window emerges via
      # washout once doses are skipped.
      name = "ETG", label = "Desogestrel POP (etonogestrel)", kind = "POP",
      dose_ug = 75, F = 0.76, Vd_L = 150, t_half_h = 25, Ka_h = 2.5,
      fm_cyp3a4 = 0.30, margin_prog = 3.2, emax_prog = 0.97,
      has_ee = FALSE,
      mucus_margin = 2.0, mucus_emax = 0.90, mucus_hill = 3.0),
    LNGP = make_drug(
      # Traditional LNG-POP suppresses ovulation in only ~half of cycles; its
      # contraceptive effect is substantially cervical-mucus-mediated. Low margin
      # + sub-maximal Emax place the perfect-use trough near the surge-block
      # threshold, so a sizeable fraction ovulate even with perfect use. The
      # cervical-mucus arm (mucus_emax) now carries its residual forgiveness on
      # the failure endpoint, so 'forgiveness' is not undefined.
      # Because ~48% of cycles ovulate, the mucus barrier must block nearly all
      # breakthrough cycles to reproduce the observed perfect-use Pearl Index
      # (~1.1/100 WY) -- hence the near-complete mucus_emax.
      name = "LNGP", label = "Levonorgestrel POP (0.03 mg)", kind = "POP",
      dose_ug = 30, F = 0.95, Vd_L = 126, t_half_h = 26, Ka_h = 2.5,
      fm_cyp3a4 = 0.44, margin_prog = 1.30, emax_prog = 0.70,   # same drug as LNG COC
      has_ee = FALSE,
      mucus_margin = 3.2, mucus_emax = 0.992, mucus_hill = 4.5),
    DRSPP = make_drug(
      # Drospirenone-only pill (Slynd, 4 mg): estrogen-free, reliably ovulation-
      # inhibiting like desogestrel but with the LONGEST missed-pill window of the
      # POPs (24 h) owing to drospirenone's long t1/2 (~31 h). Clinically a 24/4
      # regimen; modelled continuous for comparability with the other POPs (the
      # challenge is mid-pack missed active pills). CYP3A4-metabolised (fm 0.5), so
      # enzyme induction erodes its margin -- shown in the covariate analysis. A
      # comfortable margin + the long half-life make the 24 h grace window emerge
      # via slow washout once doses are skipped.
      name = "DRSPP", label = "Drospirenone POP (4 mg)", kind = "POP",
      dose_ug = 4000, F = 0.76, Vd_L = 280, t_half_h = 31, Ka_h = 2.5,
      fm_cyp3a4 = 0.50, margin_prog = 3.6, emax_prog = 0.97,
      has_ee = FALSE,
      mucus_margin = 2.2, mucus_emax = 0.96, mucus_hill = 3.0)
  )
}
