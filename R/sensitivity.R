# =============================================================================
# sensitivity.R -- Sensitivity & identifiability (Protocol Section 5.5).
# Local and global sensitivity analysis.
#
#  * LOCAL one-at-a-time +/-20% tornado (as before), now N>=200 (less noisy).
#  * GLOBAL screening via the Morris elementary-effects method (Fix: the prior
#    build did only local OAT; the protocol asked for "Sobol or Morris"). Morris
#    ranks drivers across the joint parameter space and flags interactions/
#    non-linearity (sigma vs mu*).
#
# SCENARIO PLACEMENT:
# The perturbation scenarios are now evaluated AT each agent's max-tolerated
# BOUNDARY (LNG: 2 missed at start, its escape~5% crossing; ETG: 1 missed mid,
# just past its 0-pill window), NOT deep in the saturated escape regime. In the
# an earlier placement probed LNG at 4 missed pills -- twice its tolerance -- where the
# progestin has fully washed out, so every PK/PD drug parameter had exactly zero
# elementary effect and only ovarian-dynamics parameters moved the (already
# saturated) endpoint. That was an artdefact of scenario placement. At the
# tolerance boundary the suppression margin and PK persistence re-enter alongside
# the ovarian-dynamics parameters, giving an informative ranking of what actually
# sets the endpoint. The qualitative agent ranking remains preserved under all
# perturbations; the absolute endpoint remains sensitive to the least
# data-constrained ovarian-dynamics parameters, which we state explicitly.
# =============================================================================

DRUG_PARAMS <- c("t_half_h", "Vd_L", "Ka_h", "margin_prog", "emax_prog", "hill_prog")
HPO_PARAMS  <- c("kgrow", "surge_crit", "F_ov", "E2_surge")
PARAM_LABEL <- c(t_half_h = "progestin t1/2", Vd_L = "progestin V/F", Ka_h = "absorption Ka",
                 margin_prog = "forgiveness margin", emax_prog = "Emax (surge block)",
                 hill_prog = "Hill coefficient", kgrow = "follicle growth rate",
                 surge_crit = "surge threshold", F_ov = "preovulatory size",
                 E2_surge = "E2 surge set-point")
# Scenarios placed at each agent's max-tolerated boundary (see header).
SENS_SCEN <- list(
  LNG = list(label = "COC: LNG, 2 missed at start (tolerance boundary)", n_missed = 2, position = "start"),
  ETG = list(label = "POP: desogestrel, 1 missed mid (tolerance boundary)", n_missed = 1, position = "mid"))
DELTA <- 0.20

perturb_drug <- function(drug, param, factor) {
  d <- drug
  d[[param]] <- d[[param]] * factor
  if (param == "emax_prog") d$emax_prog <- min(d$emax_prog, 0.999)
  d$ke_h <- log(2) / d$t_half_h; d$CL_Lh <- d$ke_h * d$Vd_L
  d$ee_ke_h <- log(2) / d$ee_t_half_h
  d
}

# Escape prob over population with optional drug + HPO perturbation, parallel.
sens_escape <- function(cl, drug, subs, scen, drug_param = NULL, hpo_over = NULL,
                        kgrow_factor = 1.0) {
  res <- clusterApply(cl, subs, function(s, drug, scen, hpo_over, kgrow_factor) {
    pdv <- s$pd
    if (kgrow_factor != 1.0) pdv$kgrow <- (if (is.null(pdv$kgrow)) P_default()$kgrow else pdv$kgrow) * kgrow_factor
    if (!is.null(hpo_over)) for (nm in names(hpo_over)) pdv[[nm]] <- hpo_over[[nm]]
    as.integer(isTRUE(run_scenario(drug, n_missed = scen$n_missed, position = scen$position,
                                   cl_factor = s$cl_factor, vd_factor = s$vd_factor,
                                   ee_cl_factor = s$ee_cl_factor, pd = pdv)$escaped))
  }, drug = drug, scen = scen, hpo_over = hpo_over, kgrow_factor = kgrow_factor)
  mean(unlist(res))
}

# ---- LOCAL one-at-a-time tornado ----
tornado_for <- function(cl, DRUGS, drug_key, n_sub = 300) {
  drug <- DRUGS[[drug_key]]; scen <- SENS_SCEN[[drug_key]]
  # Use the SAME nested main-grid population (tag="main") so the tornado is the
  # sensitivity OF the reported endpoint, and its base escape lands at the ~5%
  # decision boundary rather than at a fluke rate from an independent draw.
  subs <- sample_population(n_sub, drug = drug, tag = "main")
  base <- sens_escape(cl, drug, subs, scen)
  rows <- list()
  for (param in c(DRUG_PARAMS, HPO_PARAMS)) {
    if (param %in% DRUG_PARAMS) {
      e_lo <- sens_escape(cl, perturb_drug(drug, param, 1 - DELTA), subs, scen)
      e_hi <- sens_escape(cl, perturb_drug(drug, param, 1 + DELTA), subs, scen)
    } else if (param == "kgrow") {
      e_lo <- sens_escape(cl, drug, subs, scen, kgrow_factor = 1 - DELTA)
      e_hi <- sens_escape(cl, drug, subs, scen, kgrow_factor = 1 + DELTA)
    } else {
      p0 <- P_default()[[param]]
      e_lo <- sens_escape(cl, drug, subs, scen, hpo_over = setNames(list(p0 * (1 - DELTA)), param))
      e_hi <- sens_escape(cl, drug, subs, scen, hpo_over = setNames(list(p0 * (1 + DELTA)), param))
    }
    rows[[length(rows) + 1]] <- data.frame(
      drug = drug_key, scenario = scen$label, param = param, label = PARAM_LABEL[[param]],
      escape_base = round(100 * base, 1), escape_minus20 = round(100 * e_lo, 1),
      escape_plus20 = round(100 * e_hi, 1), swing_pp = round(100 * abs(e_hi - e_lo), 1),
      stringsAsFactors = FALSE)
  }
  df <- do.call(rbind, rows)
  df[order(-df$swing_pp), ]
}

# ---- GLOBAL Morris elementary-effects screening ----
# NOTE: the Morris design is randomly generated. It is seeded here so the
# screening is exactly reproducible; r is set high enough that mu* is estimated
# with a useful standard error (sigma/sqrt(r)), which is reported alongside it.
morris_global <- function(cl, DRUGS, drug_key, r = 40, n_inner = 60, seed = 20240601) {
  if (!requireNamespace("sensitivity", quietly = TRUE)) return(NULL)
  set.seed(seed)
  drug <- DRUGS[[drug_key]]; scen <- SENS_SCEN[[drug_key]]
  subs <- sample_population(n_inner, drug = drug, tag = "main")  # same nested population as the grid
  params <- c(DRUG_PARAMS, HPO_PARAMS)
  base_vals <- c(sapply(DRUG_PARAMS, function(p) drug[[p]]),
                 sapply(HPO_PARAMS, function(p) P_default()[[p]]))
  lower <- base_vals * (1 - DELTA); upper <- base_vals * (1 + DELTA)
  des <- sensitivity::morris(model = NULL, factors = params, r = r,
                             design = list(type = "oat", levels = 6, grid.jump = 3),
                             binf = lower, bsup = upper)
  X <- des$X
  y <- numeric(nrow(X))
  for (i in seq_len(nrow(X))) {
    d2 <- drug
    for (p in DRUG_PARAMS) d2[[p]] <- X[i, p]
    d2$emax_prog <- min(d2$emax_prog, 0.999)           # fractional-suppression ceiling
    d2$ke_h <- log(2) / d2$t_half_h; d2$ee_ke_h <- log(2) / d2$ee_t_half_h
    hpo_over <- as.list(X[i, HPO_PARAMS]); names(hpo_over) <- HPO_PARAMS
    kgf <- hpo_over$kgrow / P_default()$kgrow; hpo_over$kgrow <- NULL
    y[i] <- sens_escape(cl, d2, subs, scen, hpo_over = hpo_over, kgrow_factor = kgf)
  }
  des <- sensitivity::tell(des, y)
  mu_star <- apply(des$ee, 2, function(e) mean(abs(e)))
  sigma <- apply(des$ee, 2, sd)
  data.frame(drug = drug_key, scenario = scen$label, param = params,
             label = PARAM_LABEL[params],
             mu_star = round(100 * mu_star, 2), sigma = round(100 * sigma, 2),
             se_mu_star = round(100 * sigma / sqrt(r), 2), r = r, n_inner = n_inner,
             stringsAsFactors = FALSE)[order(-mu_star), ]
}
