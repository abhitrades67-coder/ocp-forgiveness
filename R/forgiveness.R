# =============================================================================
# forgiveness.R -- Forgiveness quantification + virtual population
# (Protocol Sections 4.5-4.7, 5.2).
#
#  * Virtual population: log-normal IIV on PK clearance/volume and PD axis
#    sensitivity/potency, plus CYP3A4 and body-weight covariates.
#  * Monte-Carlo ovulation-escape probability per (drug, position, n_missed) cell.
#  * Primary metric: max consecutive missed active pills tolerated before escape
#    probability exceeds 5%.
#  * Continuous forgiveness metric.
#
# DESIGN NOTES:
#  * N_SUBJECTS = 600: the integer max-tolerated endpoint does not
#    longer flips on +/-1 subject (1 subject = 0.17 pp, not 1.67 pp).
#  * Bootstrap CI on the max-tolerated integer: the headline endpoint
#    now carries uncertainty, so claimed drug/position differences are testable.
#  * Single shared population per (drug, weight, cyp), seeded independently of N
#    : the perfect-use ovulation rate is ONE number, reused by the
#    main grid and the validation module (no more 0% / 3.3% / 5.0% disagreement).
#  * MC-derived continuous forgiveness index: the missed-pill count at
#    which the escape curve crosses 5% (linear interpolation), replacing the
#    bimodal trough-vs-threshold index that collapsed to +/-24 h.
# =============================================================================

ESCAPE_THRESHOLD <- 0.05
N_SUBJECTS <- 600
IIV <- list(cl_cv = 0.30, vd_cv = 0.25, kgrow_cv = 0.15, ec_cv = 0.15)

# Deterministic 31-bit string hash (djb2) -> stable RNG seed, reproducible across
# platforms (base R Mersenne-Twister). Independent of n so populations are nested.
stable_seed <- function(...) {
  key <- paste(..., sep = "|")
  h <- 5381
  for (cc in utf8ToInt(key)) h <- (h * 33 + cc) %% 2147483647
  as.integer(h)
}

sample_population <- function(n, weight_kg = 70.0, cyp_status = "normal", drug = NULL,
                              tag = "pop") {
  seed <- stable_seed(if (is.null(drug)) "nodrug" else drug$name,
                      weight_kg, cyp_status, tag)
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed)
  on.exit({ if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv) }, add = TRUE)

  cyp_factor <- if (!is.null(drug)) cyp3a4_clearance_factor(drug$fm_cyp3a4, cyp_status) else 1.0
  # Ethinylestradiol has its own CYP3A4-dependent fraction (ee_fm_cyp3a4); under
  # normal status both factors are exactly 1, so non-CYP analyses are unchanged.
  ee_cyp_factor <- if (!is.null(drug) && !is.null(drug$ee_fm_cyp3a4) &&
                       !is.na(drug$ee_fm_cyp3a4) && isTRUE(drug$has_ee))
    cyp3a4_clearance_factor(drug$ee_fm_cyp3a4, cyp_status) else cyp_factor
  wt_v  <- weight_kg / 70.0
  wt_cl <- (weight_kg / 70.0)^0.75
  kgrow0 <- P_default()$kgrow
  lapply(seq_len(n), function(i) {
    eta_cl <- exp(rnorm(1, 0, IIV$cl_cv))   # one shared CL random effect
    cl <- eta_cl * cyp_factor * wt_cl
    ee_cl <- eta_cl * ee_cyp_factor * wt_cl
    vd <- exp(rnorm(1, 0, IIV$vd_cv)) * wt_v
    kgrow <- kgrow0 * exp(rnorm(1, 0, IIV$kgrow_cv))
    ec_scale <- exp(rnorm(1, 0, IIV$ec_cv))
    list(cl_factor = cl, ee_cl_factor = ee_cl, vd_factor = vd,
         pd = list(kgrow = kgrow, ec_scale = ec_scale))
  })
}

# Escape indicator vector (0/1) over a population -- keep per-subject so we can
# bootstrap. scen_kwargs passed through to run_scenario.
escape_vector <- function(drug, n_missed, position, subjects, scen_kwargs = list()) {
  vapply(subjects, function(s) {
    args <- c(list(drug = drug, n_missed = n_missed, position = position,
                   cl_factor = s$cl_factor, vd_factor = s$vd_factor,
                   ee_cl_factor = s$ee_cl_factor, pd = s$pd),
              scen_kwargs)
    as.integer(isTRUE(do.call(run_scenario, args)$escaped))
  }, integer(1))
}

wilson_ci <- function(k, n, z = 1.96) {
  if (n == 0) return(c(NA, NA))
  p <- k / n; denom <- 1 + z * z / n
  centre <- (p + z * z / (2 * n)) / denom
  half <- (z / denom) * sqrt(p * (1 - p) / n + z * z / (4 * n * n))
  c(max(0, centre - half), min(1, centre + half))
}

# Max consecutive missed pills with escape prob <= threshold (first breach stops).
max_tolerated_from <- function(ns, probs, threshold = ESCAPE_THRESHOLD) {
  o <- order(ns); ns <- ns[o]; probs <- probs[o]
  tol <- -1L
  for (i in seq_along(ns)) {
    if (probs[i] <= threshold) tol <- as.integer(ns[i]) else break
  }
  tol
}

# Continuous forgiveness index: missed-pill count at which the escape
# curve first crosses `threshold`, by linear interpolation between grid points.
# Returns the interpolated count (>= 0); if perfect use already exceeds the
# threshold returns 0 with attribute defined=FALSE (ovulation not reliably
# suppressed -- the mucus arm carries residual protection, see translate.R).
escape_threshold_crossing <- function(ns, probs, threshold = ESCAPE_THRESHOLD) {
  o <- order(ns); ns <- ns[o]; probs <- probs[o]
  if (probs[1] > threshold) return(structure(0, defined = FALSE))
  cross <- NA_real_
  for (i in seq_len(length(ns) - 1L)) {
    if (probs[i] <= threshold && probs[i + 1] > threshold) {
      frac <- (threshold - probs[i]) / (probs[i + 1] - probs[i])
      cross <- ns[i] + frac * (ns[i + 1] - ns[i]); break
    }
  }
  if (is.na(cross)) cross <- max(ns)   # never crosses within the grid
  structure(cross, defined = TRUE)
}

# Bootstrap CI for the max-tolerated integer endpoint. Resamples the
# per-subject escape indicator matrix (rows = subjects, cols = n_missed grid),
# preserving the common-random-numbers structure across n_missed.
bootstrap_max_tolerated <- function(esc_matrix, ns, threshold = ESCAPE_THRESHOLD,
                                    B = 2000, seed = 1) {
  set.seed(seed)
  N <- nrow(esc_matrix)
  pt <- max_tolerated_from(ns, colMeans(esc_matrix), threshold)
  cr <- as.numeric(escape_threshold_crossing(ns, colMeans(esc_matrix), threshold))
  boots_mt <- integer(B); boots_cr <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(N, N, replace = TRUE)
    pr <- colMeans(esc_matrix[idx, , drop = FALSE])
    boots_mt[b] <- max_tolerated_from(ns, pr, threshold)
    boots_cr[b] <- as.numeric(escape_threshold_crossing(ns, pr, threshold))
  }
  list(point = pt, lo = as.integer(quantile(boots_mt, 0.025, type = 1)),
       hi = as.integer(quantile(boots_mt, 0.975, type = 1)),
       cross = cr, cross_lo = unname(quantile(boots_cr, 0.025)),
       cross_hi = unname(quantile(boots_cr, 0.975)))
}

# Paired (common-random-number) bootstrap test of the pack-position effect (H2).
# `mats` is a named list
# position -> N x length(ns) 0/1 escape matrix, ALL sampled on the SAME shared
# population -- sample_population() seeds on (drug, weight, cyp, tag) and IGNORES
# position, so row i is the identical virtual subject across positions. Resampling
# subject indices ONCE per bootstrap and applying that same resample to every
# position yields a PAIRED difference in the continuous forgiveness index
# (interpolated missed-pills-to-threshold) that cancels the shared subject-level
# variance. This is materially more powerful than comparing the marginal per-
# position bootstrap CIs (which overlap even when the paired difference is
# consistently signed), so it can resolve a sub-pill position effect the integer
# endpoint and the marginal CIs cannot. Returns per-position point crossings and,
# for every position pair, the paired difference with its 95% CI and the bootstrap
# probability that the difference is >= 0.
paired_position_bootstrap <- function(mats, ns, threshold = ESCAPE_THRESHOLD,
                                      B = 2000, seed = 1) {
  positions <- names(mats)
  N <- nrow(mats[[1]])
  stopifnot(all(vapply(mats, nrow, integer(1)) == N))   # must be the same population
  cross_of <- function(mat, idx)
    as.numeric(escape_threshold_crossing(ns, colMeans(mat[idx, , drop = FALSE]), threshold))
  pt <- vapply(positions, function(p) cross_of(mats[[p]], seq_len(N)), numeric(1))
  set.seed(seed)
  bcross <- matrix(NA_real_, B, length(positions), dimnames = list(NULL, positions))
  for (b in seq_len(B)) {
    idx <- sample.int(N, N, replace = TRUE)
    for (p in positions) bcross[b, p] <- cross_of(mats[[p]], idx)
  }
  pairs <- utils::combn(positions, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pr) {
    d <- bcross[, pr[1]] - bcross[, pr[2]]
    lo <- unname(quantile(d, 0.025)); hi <- unname(quantile(d, 0.975))
    data.frame(contrast = paste(pr[1], "-", pr[2]),
               diff = round(pt[pr[1]] - pt[pr[2]], 3),
               ci_lo = round(lo, 3), ci_hi = round(hi, 3),
               p_diff_ge0 = round(mean(d >= 0), 3),
               excludes_zero = (lo > 0) || (hi < 0), stringsAsFactors = FALSE)
  })
  list(point = pt, contrasts = do.call(rbind, rows))
}
