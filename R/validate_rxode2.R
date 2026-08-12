# =============================================================================
# validate_rxode2.R -- Cross-validate the rxode2 HPO engine against the
# production pure-R RK4 engine. Run:  Rscript R/validate_rxode2.R
#
# Confirms the model ports faithfully to rxode2 (a protocol-named toolchain):
# identical ovulation counts and near-identical trajectories on scenarios that
# exercise (i) autonomous cycling + events, (ii) strong COC suppression with the
# EE covariate, (iii) event firing under POP partial washout, and (iv) an
# extended-HFI missed-pill escape. Writes results_R/tables/rxode2_crossvalidation.csv
# and results_R/figures/figS1_rxode2_crossvalidation.png.
# =============================================================================
suppressMessages({
  library(rxode2); library(ggplot2)
  source("R/pk.R"); source("R/drugs.R"); source("R/hpo.R")
  source("R/simulate.R"); source("R/hpo_rxode2.R")
})
DRUGS <- build_drugs()
dir.create("results_R/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results_R/tables",  recursive = TRUE, showWarnings = FALSE)

# agreement metrics between two engine outputs on a shared output grid
agree <- function(a, b, state_row) {
  x <- a$Y[state_row, ]; y <- b$Y[state_row, ]
  n <- min(length(x), length(y)); x <- x[1:n]; y <- y[1:n]
  c(rmse = sqrt(mean((x - y)^2)), maxabs = max(abs(x - y)),
    scale = diff(range(x)))
}
ovmatch <- function(a, b, tol = 0.6) {
  oa <- a$ov_times; ob <- b$ov_times
  matched <- if (length(oa) == length(ob) && length(oa) > 0) all(abs(oa - ob) <= tol) else length(oa) == length(ob)
  list(n_rk4 = length(oa), n_rx = length(ob), max_dt = if (length(oa) == length(ob) && length(oa) > 0) max(abs(oa - ob)) else NA,
       ov_consistent = matched)
}

rows <- list(); traces <- list()
record <- function(tag, rk4, rx) {
  fo <- agree(rk4, rx, iF); eo <- agree(rk4, rx, iE2)
  om <- ovmatch(rk4, rx)
  rows[[tag]] <<- data.frame(
    scenario = tag, ov_rk4 = om$n_rk4, ov_rx = om$n_rx,
    ov_max_dt_d = round(om$max_dt, 3), ov_consistent = om$ov_consistent,
    follicle_rmse = round(fo["rmse"], 4), follicle_maxabs = round(fo["maxabs"], 4),
    follicle_scale = round(fo["scale"], 2),
    e2_rmse = round(eo["rmse"], 4), e2_maxabs = round(eo["maxabs"], 4),
    stringsAsFactors = FALSE, row.names = NULL)
  traces[[tag]] <<- data.frame(t = rk4$t, rk4 = rk4$Y[iF, ],
                               rx = approx(rx$t, rx$Y[iF, ], rk4$t, rule = 2)$y,
                               scenario = tag)
}

cat("Compiling rxode2 model (first use) ...\n"); invisible(.hpo_rx_model())

# ---- A. Untreated 90 d: autonomous cycling + ovulation events, no drug ----
tg <- seq(0, 90, by = 0.05)
profU <- list(t = tg, c_prog = numeric(length(tg)), c_ee = numeric(length(tg)))
A_rk4 <- hpo_simulate   (profU, 90, 1e9, 0, 2, NA, 0, 2)
A_rx  <- hpo_simulate_rx(profU, 90, 1e9, 0, 2, NA, 0, 2)
record("A_untreated", A_rk4, A_rx)
cat(sprintf("A untreated: RK4 ov=%s | rx ov=%s\n",
            paste(round(A_rk4$ov_times,1), collapse=","), paste(round(A_rx$ov_times,1), collapse=",")))

# ---- B. Perfect-use LNG-COC 112 d: strong suppression + EE covariate ----
schedL <- apply_adherence(build_pack_schedule(21, 7, n_packs = 4))
profL <- make_profile(DRUGS$LNG, schedL, 112)
ecpL <- ec50_prog(DRUGS$LNG); eceL <- ec50_ee(DRUGS$LNG)
B_rk4 <- hpo_simulate   (profL, 112, ecpL, DRUGS$LNG$emax_prog, DRUGS$LNG$hill_prog, eceL, DRUGS$LNG$ee_emax, DRUGS$LNG$ee_hill)
B_rx  <- hpo_simulate_rx(profL, 112, ecpL, DRUGS$LNG$emax_prog, DRUGS$LNG$hill_prog, eceL, DRUGS$LNG$ee_emax, DRUGS$LNG$ee_hill)
record("B_perfect_LNG", B_rk4, B_rx)
cat(sprintf("B perfect LNG: RK4 ov=%d | rx ov=%d\n", length(B_rk4$ov_times), length(B_rx$ov_times)))

# ---- C. Perfect-use LNGP-POP 112 d: event firing under partial suppression ----
schedP <- apply_adherence(build_pack_schedule(28, 0, n_packs = 4))
profP <- make_profile(DRUGS$LNGP, schedP, 112)
ecpP <- ec50_prog(DRUGS$LNGP)
C_rk4 <- hpo_simulate   (profP, 112, ecpP, DRUGS$LNGP$emax_prog, DRUGS$LNGP$hill_prog, NA, 0, 2)
C_rx  <- hpo_simulate_rx(profP, 112, ecpP, DRUGS$LNGP$emax_prog, DRUGS$LNGP$hill_prog, NA, 0, 2)
record("C_perfect_LNGP", C_rk4, C_rx)
cat(sprintf("C perfect LNGP: RK4 ov=%d | rx ov=%d\n", length(C_rk4$ov_times), length(C_rx$ov_times)))

# ---- D. LNG missed 6 at start (extended-HFI escape) via the full run_scenario path ----
D_rk4 <- run_scenario(DRUGS$LNG, n_missed = 6, position = "start", challenge_pack = 2,
                      follow_days = 34, return_trace = TRUE)
hpo_simulate_orig <- hpo_simulate
hpo_simulate <<- hpo_simulate_rx            # rebind so run_scenario uses the rx engine
D_rx <- run_scenario(DRUGS$LNG, n_missed = 6, position = "start", challenge_pack = 2,
                     follow_days = 34, return_trace = TRUE)
hpo_simulate <<- hpo_simulate_orig          # restore
record("D_LNG_miss6_start", D_rk4, D_rx)
cat(sprintf("D LNG miss6 start: RK4 escaped=%s ov=%d | rx escaped=%s ov=%d\n",
            D_rk4$escaped, length(D_rk4$ov_times), D_rx$escaped, length(D_rx$ov_times)))

# ---- assemble ----
tbl <- do.call(rbind, rows)
write.csv(tbl, "results_R/tables/rxode2_crossvalidation.csv", row.names = FALSE)
cat("\n==== rxode2 vs pure-R RK4 cross-validation ====\n"); print(tbl, row.names = FALSE)

td <- do.call(rbind, traces)
pS1 <- ggplot(td, aes(t)) +
  geom_line(aes(y = rk4, colour = "pure-R RK4"), linewidth = 0.7) +
  geom_line(aes(y = rx,  colour = "rxode2"), linetype = 2, linewidth = 0.7) +
  facet_wrap(~scenario, scales = "free", ncol = 2) +
  scale_colour_manual(values = c("pure-R RK4" = "#1f77b4", "rxode2" = "#d62728"), name = NULL) +
  labs(x = "time (days)", y = "follicle (a.u.)",
       title = "Engine cross-validation: HPO follicle trajectory, pure-R RK4 vs rxode2") +
  theme_bw(base_size = 11) + theme(legend.position = "top")
.figdir <- Sys.getenv("OCP_FIG_DIR", "results_R/figures"); dir.create(.figdir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(.figdir, "figS1_rxode2_crossvalidation.png"), pS1, width = 9, height = 6,
       dpi = as.numeric(Sys.getenv("OCP_FIG_DPI", "130")))
cat("\nWrote results_R/tables/rxode2_crossvalidation.csv and figures/figS1_rxode2_crossvalidation.png\n")
cat(if (all(tbl$ov_consistent)) "ALL_OV_CONSISTENT\n" else "OV_MISMATCH\n")
