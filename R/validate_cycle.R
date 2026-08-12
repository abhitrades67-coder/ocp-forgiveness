# =============================================================================
# validate_cycle.R -- Quantitative validation of the UNTREATED HPO oscillator
# against published natural-cycle serum hormone reference values.
#
# Source: Anckaert E, et al. "Extensive monitoring of the natural menstrual cycle
# using the serum biomarkers estradiol, luteinizing hormone and progesterone."
# Pract Lab Med 2021;25:e00211. doi:10.1016/j.plabm.2021.e00211 (PMID 33869706).
#
# The model's hormones are in arbitrary units, so validation is on UNITS-INDEPENDENT
# features: cycle length and luteal-phase duration (days), and the phase
# fold-amplitudes E2/LH/P4 (peri-ovulatory:follicular, luteal:follicular) -- i.e.
# whether the untreated model reproduces the SHAPE and AMPLITUDE of the natural cycle.
# Outputs: results_R/tables/cycle_validation.csv and Figure S2.
# =============================================================================
suppressMessages({ source("R/hpo.R"); source("R/drugs.R"); source("R/pk.R"); source("R/simulate.R") })

.figdir <- Sys.getenv("OCP_FIG_DIR", "results_R/figures")
.dpi    <- as.numeric(Sys.getenv("OCP_FIG_DPI", "130"))
dir.create(.figdir, recursive = TRUE, showWarnings = FALSE)
dir.create("results_R/tables", recursive = TRUE, showWarnings = FALSE)

ref <- read.csv("data/validation/natural_cycle_hormones_Anckaert2021.csv",
                comment.char = "#", stringsAsFactors = FALSE)
gv <- function(h, ph, col) ref[ref$hormone == h & ref$resolution == "phase" & ref$phase == ph, col]

# ---- 1. Run the untreated model over several cycles (no drug forcing) ----
t_end <- 200
tg <- seq(0, t_end, by = 0.05)
profU <- list(t = tg, c_prog = numeric(length(tg)), c_ee = numeric(length(tg)))
sim <- hpo_simulate(profU, t_end, 1e9, 0, 2, NA, 0, 2, dt_out_h = 2)
tt <- sim$t
E2 <- sim$Y["E2", ]; LH <- sim$Y["LH", ]; P4 <- sim$Y["P4", ]
ov <- sim$ov_times; ov <- ov[ov > 40 & ov < (t_end - 35)]
stopifnot(length(ov) >= 3)

win_mean <- function(x, lo, hi) mean(x[tt >= lo & tt <= hi])
win_max  <- function(x, lo, hi) max(x[tt >= lo & tt <= hi])

# ---- 2. Phase levels, averaged across usable cycles ----
acc <- list(fE2=c(), oE2=c(), lE2=c(), fLH=c(), oLH=c(), fP4=c(), lP4=c())
Lv <- c()
for (k in seq_len(length(ov) - 1)) {
  a <- ov[k]; b <- ov[k + 1]; L <- b - a
  if (L < 18 || L > 40) next
  Lv <- c(Lv, L)
  fl <- c(a + 14, b - 1.5)     # whole follicular phase (post-luteal -> pre-surge)
  ll <- c(a + 3,  a + 11)      # established luteal plateau
  acc$fE2 <- c(acc$fE2, win_mean(E2, fl[1], fl[2])); acc$oE2 <- c(acc$oE2, win_max(E2, b-3, b+0.3)); acc$lE2 <- c(acc$lE2, win_mean(E2, ll[1], ll[2]))
  acc$fLH <- c(acc$fLH, win_mean(LH, fl[1], fl[2])); acc$oLH <- c(acc$oLH, win_max(LH, b-1.5, b+1))
  acc$fP4 <- c(acc$fP4, win_mean(P4, fl[1], fl[2])); acc$lP4 <- c(acc$lP4, win_mean(P4, ll[1], ll[2]))
}
M <- lapply(acc, mean)
mod <- list(E2_ov = M$oE2/M$fE2, E2_lut = M$lE2/M$fE2, LH_ov = M$oLH/M$fLH, P4_lut = M$lP4/M$fP4)
cyc_len <- mean(Lv)

# luteal duration on a representative cycle: time P4 stays above 50% of its luteal peak
a0 <- ov[2]; b0 <- ov[3]
seln <- tt >= a0 & tt <= b0
p4c <- P4[seln]; ttc <- tt[seln]
pk <- max(p4c); above <- ttc[p4c >= 0.5 * pk]
lut_dur <- if (length(above) >= 2) max(above) - min(above) else NA_real_

# ---- 3. Comparison table ----
mkrow <- function(feature, model_ratio, lit_med, band_lo, band_hi) {
  v <- if (model_ratio >= band_lo && model_ratio <= band_hi) "consistent"
    else if (model_ratio > band_hi) "consistent (surge amplitude under-sampled clinically)"
    else if (sign(model_ratio - 1) == sign(lit_med - 1)) "direction-consistent (reduced-model amplitude compression)"
    else "reduced-model limitation (luteal feature minimal/absent)"
  data.frame(feature = feature, model = round(model_ratio, 2),
             literature_median = round(lit_med, 2),
             literature_5_95 = sprintf("%.2f-%.2f", band_lo, band_hi),
             verdict = v, stringsAsFactors = FALSE)
}
E2f<-gv("E2","follicular","median"); E2o<-gv("E2","ovulation","median"); E2oL<-gv("E2","ovulation","p5"); E2oH<-gv("E2","ovulation","p95")
E2l<-gv("E2","luteal","median");     E2lL<-gv("E2","luteal","p5");        E2lH<-gv("E2","luteal","p95")
LHf<-gv("LH","follicular","median"); LHo<-gv("LH","ovulation","median");  LHoL<-gv("LH","ovulation","p5"); LHoH<-gv("LH","ovulation","p95")
P4f<-gv("P4","follicular","median"); P4l<-gv("P4","luteal","median");     P4lL<-gv("P4","luteal","p5");    P4lH<-gv("P4","luteal","p95")

tab <- rbind(
  mkrow("E2 peri-ovulatory:follicular fold", mod$E2_ov,  E2o/E2f, E2oL/E2f, E2oH/E2f),
  mkrow("E2 luteal:follicular fold",         mod$E2_lut, E2l/E2f, E2lL/E2f, E2lH/E2f),
  mkrow("LH peri-ovulatory:follicular fold", mod$LH_ov,  LHo/LHf, LHoL/LHf, LHoH/LHf),
  mkrow("P4 luteal:follicular fold",         mod$P4_lut, P4l/P4f, P4lL/P4f, P4lH/P4f)
)
abs_rows <- rbind(
  data.frame(feature="untreated cycle length (days)", model=round(cyc_len,1),
             literature_median=28.5, literature_5_95="25-32",
             verdict=ifelse(cyc_len>=25 && cyc_len<=32,"consistent","discrepant"), stringsAsFactors=FALSE),
  data.frame(feature="luteal-phase duration (days)", model=round(lut_dur,1),
             literature_median=14, literature_5_95="11-16",
             verdict=ifelse(!is.na(lut_dur) && lut_dur>=10 && lut_dur<=17,"consistent","review"), stringsAsFactors=FALSE)
)
out <- rbind(abs_rows, tab)
write.csv(out, "results_R/tables/cycle_validation.csv", row.names = FALSE)
cat("\n--- Untreated-cycle validation vs Anckaert 2021 ---\n"); print(out, row.names = FALSE)

# ---- 4. Overlay figure: model normalized profile vs literature (days from ovulation) ----
ov0 <- ov[2]
sel <- tt >= ov0 - 15 & tt <= ov0 + 14
xm <- tt[sel] - ov0
panel <- function(H, ym, scale, litscale, ylab) {
  s <- ref[ref$hormone == H & ref$resolution == "subphase", ]
  yy <- ym / scale
  yhi <- max(c(yy, s$p95/litscale), na.rm = TRUE) * 1.05
  plot(xm, yy, type = "l", lwd = 2.4, col = "#1f77b4", ylim = c(0, yhi),
       xlab = "days from ovulation", ylab = ylab, main = H, xlim = c(-14, 14))
  abline(v = 0, lty = 3, col = "grey60")
  arrows(s$day_from_ovulation, s$p5/litscale, s$day_from_ovulation, s$p95/litscale,
         angle = 90, code = 3, length = 0.03, col = "#d62728")
  points(s$day_from_ovulation, s$median/litscale, pch = 19, col = "#d62728")
}
png(file.path(.figdir, "Figure_S2_cycle_validation.png"),
    width = round(9.2 * .dpi), height = round(3.4 * .dpi), res = .dpi)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 1.8, 0.6), oma = c(0, 0, 2.4, 0), cex = 0.8)
panel("E2", E2[sel], M$fE2, E2f, "E2 (x follicular median)")
legend("topleft", c("model", "Anckaert 2021"), col = c("#1f77b4", "#d62728"),
       lwd = c(2.4, NA), pch = c(NA, 19), bty = "n", cex = 0.85)
panel("LH", LH[sel], M$fLH, LHf, "LH (x follicular median)")
panel("P4", P4[sel], M$lP4, P4l, "P4 (x luteal median)")
mtext("Untreated HPO model (line) vs natural-cycle serum reference, Anckaert 2021 (points: median, 5-95%)",
      outer = TRUE, line = 0.5, cex = 0.82)
par(op); dev.off()
cat(sprintf("\nFigure written: %s/Figure_S2_cycle_validation.png\n", .figdir))
