# =============================================================================
# run_protocol.R -- Journal R pipeline orchestrator (in silico execution of
# OCP_Forgiveness_Protocol.docx).
#
# Aim 1  Model development & qualification
# Aim 2  Forgiveness quantification (N=600, bootstrap-CI max-tolerated, MC index)
# Aim 3  Determinants + clinical translation (CYP3A4, weight, regimen, lateness,
#        recovery action, OR->Pearl-Index with mucus arm), validation, sensitivity
#        (local + Morris global).
#
# CHECKPOINTING (2nd-round): every expensive Monte-Carlo unit is memoised to
# results_R/_ckpt/*.rds. On restart, completed units are loaded and skipped, so
# the run RESUMES and completes across multiple invocations (needed where a hard
# background wall-clock cap kills a single ~70-min run). Just relaunch
# `Rscript run_protocol.R` until it prints DONE; the checkpoint dir is cleared on
# successful completion.
#
# Outputs: results_R/tables/*.csv, results_R/bundle.rds (for figures + manuscript).
# Usage:   Rscript run_protocol.R            (full: N=600)
#          NQUICK=1 Rscript run_protocol.R   (fast smoke test: small N)
# =============================================================================
suppressMessages({
  source("R/pk.R"); source("R/drugs.R"); source("R/hpo.R"); source("R/simulate.R")
  source("R/forgiveness.R"); source("R/parallel_util.R"); source("R/validate.R")
  source("R/sensitivity.R"); source("R/recovery.R"); source("R/translate.R")
})
dir.create("results_R/tables", recursive = TRUE, showWarnings = FALSE)
PROG <- "results_R/run_progress.txt"
say <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)),
                           file = PROG, append = TRUE) }

QUICK <- nzchar(Sys.getenv("NQUICK"))
N_MAIN <- if (QUICK) 40 else 600
N_COV  <- if (QUICK) 40 else 250
N_SEC  <- if (QUICK) 40 else 250
N_SENS <- if (QUICK) 30 else 300
N_REC  <- if (QUICK) 30 else 250
NS <- 0:7
COC <- c("LNG", "DRSP"); POP <- c("ETG", "LNGP", "DRSPP"); ALL <- c(COC, POP)
DRUGS <- build_drugs()
t_start <- Sys.time()

# ---- checkpoint infrastructure ----
CKPT <- "results_R/_ckpt"; dir.create(CKPT, recursive = TRUE, showWarnings = FALSE)
RT_FILE <- file.path(CKPT, "_runtime.rds")
# Per-invocation wall-clock budget: exit gracefully (exit 0, all checkpoints
# saved) before the environment's hard background kill (~30 min observed), so the
# run resumes cleanly on relaunch. Longest single unit is ~5 min, so 1300s leaves
# comfortable headroom. QUICK completes in one invocation (huge budget).
.bud <- Sys.getenv("OCP_BUDGET_S")
BUDGET_S <- if (nzchar(.bud) && !is.na(suppressWarnings(as.numeric(.bud)))) as.numeric(.bud) else if (QUICK) 1e6 else 1300
say("=== invocation start (N_MAIN=", N_MAIN, ", budget=", round(BUDGET_S), "s) ===")
.cl <- NULL
get_cl <- function() {
  if (is.null(.cl)) { .cl <<- make_model_cluster(); say("cluster up: ", length(.cl), " workers") }
  .cl
}
graceful_exit <- function(after) {
  if (!is.null(.cl)) try(stopCluster(.cl), silent = TRUE)
  prev <- if (file.exists(RT_FILE)) readRDS(RT_FILE) else 0
  saveRDS(prev + as.numeric(Sys.time() - t_start, units = "secs"), RT_FILE)
  say("[ckpt] budget reached after ", after, " -- graceful exit; relaunch to resume.")
  cat(sprintf("BUDGET reached after %s; checkpoint saved, relaunch to resume.\n", after))
  quit(save = "no", status = 0)
}
cache <- function(name, fn) {
  f <- file.path(CKPT, paste0(name, ".rds"))
  if (file.exists(f)) {
    val <- tryCatch(readRDS(f), error = function(e) NULL)   # tolerate a torn write
    if (!is.null(val)) { say("[ckpt] load ", name); return(val) }
    say("[ckpt] corrupt ", name, " -- recomputing")
  }
  say("[ckpt] compute ", name); val <- fn()
  tmp <- paste0(f, ".tmp"); saveRDS(val, tmp); file.rename(tmp, f)   # atomic
  say("[ckpt] saved ", name)
  if (isTRUE(as.numeric(Sys.time() - t_start, units = "secs") > BUDGET_S)) graceful_exit(name)
  val
}

# --- helper: a (drug,position) series over a shared population -> matrix + rows ---
run_series <- function(drug_key, position, N, weight = 70, cyp = "normal",
                       analysis = "main", regimen = NULL, scattered = FALSE,
                       challenge_pack = 2, ns = NS) {
  drug <- DRUGS[[drug_key]]
  subs <- sample_population(N, weight_kg = weight, cyp_status = cyp, drug = drug, tag = "main")
  sk <- list(scattered = scattered, challenge_pack = challenge_pack)
  if (!is.null(regimen)) sk$regimen <- regimen
  mat <- sapply(ns, function(n) escape_prob_par(get_cl(), drug_key, n, position, subs, sk))
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = N)
  probs <- colMeans(mat)
  rows <- lapply(seq_along(ns), function(i) {
    ci <- wilson_ci(sum(mat[, i]), N)
    data.frame(analysis = analysis, drug = drug_key, position = position, n_missed = ns[i],
               cyp = cyp, weight = weight, regimen = ifelse(is.null(regimen), regimen_for(drug), regimen),
               delay = NA_real_, scattered = scattered, challenge_pack = challenge_pack, N = N,
               n_esc = sum(mat[, i]), escape_prob = probs[i],
               escape_lo = ci[1], escape_hi = ci[2], stringsAsFactors = FALSE)
  })
  list(grid = do.call(rbind, rows), mat = mat, ns = ns, drug = drug_key, position = position)
}

grid_list <- list(); mt_list <- list(); series_store <- list()

# ---- Aim 2 main grid (each cell checkpointed) ----
say("Aim 2 main grid ...")
main_cells <- c(lapply(COC, function(dk) c(dk, "start")), lapply(COC, function(dk) c(dk, "mid")),
                lapply(COC, function(dk) c(dk, "late")), lapply(POP, function(dk) c(dk, "mid")))
main_cells <- main_cells[order(sapply(main_cells, function(x) match(x[1], ALL)))]
for (cell in main_cells) {
  dk <- cell[1]; pos <- cell[2]
  s <- cache(paste0("main_", dk, "_", pos), function() run_series(dk, pos, N_MAIN, analysis = "main"))
  grid_list[[length(grid_list) + 1]] <- s$grid
  bt <- bootstrap_max_tolerated(s$mat, s$ns)
  mt_list[[length(mt_list) + 1]] <- data.frame(
    drug = dk, label = DRUGS[[dk]]$label, kind = DRUGS[[dk]]$kind, position = pos,
    max_tolerated = bt$point, mt_lo = bt$lo, mt_hi = bt$hi,
    fi_cross = round(bt$cross, 2), fi_lo = round(bt$cross_lo, 2),
    fi_hi = round(bt$cross_hi, 2), stringsAsFactors = FALSE)
  series_store[[paste(dk, pos)]] <- s$mat
  say("  main ", dk, " ", pos, " done")
}

# ---- H2 paired pack-position test: a common-random-
# number paired bootstrap of the continuous forgiveness index across start/mid/
# late, on the identical shared population (no extra Monte-Carlo -- reuses the
# main-grid matrices in series_store). Resolves the sub-pill position effect the
# integer endpoint and the overlapping marginal CIs cannot. ----
say("H2 paired position test ...")
h2_paired <- do.call(rbind, lapply(COC, function(dk) {
  mats <- list(start = series_store[[paste(dk, "start")]],
               mid   = series_store[[paste(dk, "mid")]],
               late  = series_store[[paste(dk, "late")]])
  if (any(vapply(mats, is.null, logical(1)))) return(NULL)
  pb <- paired_position_bootstrap(mats, NS)
  df <- pb$contrasts
  df$drug <- dk
  df$fi_start <- round(pb$point[["start"]], 3)
  df$fi_mid <- round(pb$point[["mid"]], 3)
  df$fi_late <- round(pb$point[["late"]], 3)
  df
}))
write.csv(h2_paired, "results_R/tables/h2_paired.csv", row.names = FALSE)
say("H2 paired written")

# ---- Aim 3 covariates: CYP3A4 + weight (each cell checkpointed) ----
say("Aim 3 covariates ...")
for (dk in c("LNG", "DRSP", "ETG", "DRSPP")) {
  pos <- if (DRUGS[[dk]]$kind == "POP") "mid" else "start"
  for (cyp in c("induced", "inhibited"))
    grid_list[[length(grid_list) + 1]] <-
      cache(paste0("cyp_", dk, "_", cyp), function() run_series(dk, pos, N_COV, cyp = cyp, analysis = "cyp")$grid)
  for (wt in c(90, 110))
    grid_list[[length(grid_list) + 1]] <-
      cache(paste0("wt_", dk, "_", wt), function() run_series(dk, pos, N_COV, weight = wt, analysis = "weight")$grid)
  say("  covariate ", dk, " done")
}

# ---- secondary: regimen, firstpack, scattered ----
say("secondary analyses ...")
for (reg in c("21/7", "24/4", "26/2")) {
  regkey <- gsub("/", "-", reg)
  grid_list[[length(grid_list) + 1]] <-
    cache(paste0("reg_", regkey), function() run_series("LNG", "start", N_SEC, regimen = reg, analysis = "regimen")$grid)
}
grid_list[[length(grid_list) + 1]] <-
  cache("firstpack", function() run_series("LNG", "start", N_SEC, challenge_pack = 0, analysis = "firstpack")$grid)
grid_list[[length(grid_list) + 1]] <-
  cache("scattered", function() run_series("LNG", "mid", N_SEC, scattered = TRUE, analysis = "scattered")$grid)

# ---- POP late-dose analysis (lateness in hours; pills taken but delayed) ----
say("POP late-dose ...")
for (dk in POP) {
  df_late <- cache(paste0("late_", dk), function() {
    drug <- DRUGS[[dk]]; subs <- sample_population(N_SEC, drug = drug, tag = "main")
    rr <- list()
    for (dh in c(3, 6, 12, 24, 36)) {
      vec <- escape_prob_par(get_cl(), dk, 0, "mid", subs, list(delay_h = dh, challenge_pack = 2))
      ci <- wilson_ci(sum(vec), N_SEC)
      rr[[length(rr) + 1]] <- data.frame(
        analysis = "late", drug = dk, position = "mid", n_missed = 0, cyp = "normal",
        weight = 70, regimen = "continuous", delay = dh, scattered = FALSE,
        challenge_pack = 2, N = N_SEC, n_esc = sum(vec), escape_prob = mean(vec),
        escape_lo = ci[1], escape_hi = ci[2], stringsAsFactors = FALSE)
    }
    do.call(rbind, rr)
  })
  grid_list[[length(grid_list) + 1]] <- df_late
  say("  late ", dk, " done")
}

escape_grid <- do.call(rbind, grid_list)
mt <- do.call(rbind, mt_list)
write.csv(escape_grid, "results_R/tables/escape_probability_grid.csv", row.names = FALSE)
write.csv(mt, "results_R/tables/max_tolerated_missed.csv", row.names = FALSE)
say("grid + max-tolerated written")

# ---- guideline concordance (CI-aware) ----
gc <- do.call(rbind, lapply(seq_len(nrow(mt)), function(i) {
  r <- mt[i, ]; guide <- if (r$kind == "POP") 1 else 2; model <- r$max_tolerated
  verdict <- if (model < 0) "ovulation not reliably suppressed at perfect use (ovulation surrogate undefined; mucus arm carries residual cover)"
    else if (guide >= r$mt_lo && guide <= r$mt_hi) "consistent with the generic rule (within bootstrap CI)"
    else if (model > guide) "ovulation surrogate tolerates more than the generic rule"
    else "ovulation surrogate tolerates fewer than the generic rule"
  data.frame(drug = r$drug, position = r$position, guideline_max = guide,
             model_max = model, model_ci = sprintf("%d-%d", r$mt_lo, r$mt_hi),
             verdict = verdict, stringsAsFactors = FALSE)
}))
write.csv(gc, "results_R/tables/guideline_concordance.csv", row.names = FALSE)

# ---- continuous forgiveness index (MC-derived) from main series ----
fi <- do.call(rbind, lapply(names(series_store), function(key) {
  parts <- strsplit(key, " ")[[1]]; dk <- parts[1]; pos <- parts[2]
  if (pos != (if (DRUGS[[dk]]$kind == "POP") "mid" else "start")) return(NULL)
  mat <- series_store[[key]]; probs <- colMeans(mat)
  cr <- escape_threshold_crossing(NS, probs)
  data.frame(drug = dk, label = DRUGS[[dk]]$label,
             missed_to_5pct = round(as.numeric(cr), 2),
             defined = attr(cr, "defined"), perfect_use_escape_pct = round(100 * probs[1], 1),
             stringsAsFactors = FALSE)
}))
write.csv(fi, "results_R/tables/forgiveness_index.csv", row.names = FALSE)

# ---- escape-threshold robustness (2nd-round review): max-tolerated recomputed
# at 1% / 5% / 10% escape, to show the agent RANKING is stable even though the
# integer endpoint depends on the (arbitrary) 5% cutoff. Uses the already-
# simulated main series (no extra Monte-Carlo). ----
THRset <- c(0.01, 0.05, 0.10)
thr_tbl <- do.call(rbind, lapply(names(series_store), function(key) {
  parts <- strsplit(key, " ")[[1]]; dk <- parts[1]; pos <- parts[2]
  if (pos != (if (DRUGS[[dk]]$kind == "POP") "mid" else "start")) return(NULL)
  probs <- colMeans(series_store[[key]])
  mt_at <- sapply(THRset, function(th) max_tolerated_from(NS, probs, th))
  ci_at <- sapply(THRset, function(th) as.numeric(escape_threshold_crossing(NS, probs, th)))
  data.frame(drug = dk, position = pos,
             mt_1pct = mt_at[1], mt_5pct = mt_at[2], mt_10pct = mt_at[3],
             cross_1pct = round(ci_at[1], 2), cross_5pct = round(ci_at[2], 2),
             cross_10pct = round(ci_at[3], 2), stringsAsFactors = FALSE)
}))
thr_tbl$drug <- factor(thr_tbl$drug, levels = c("LNG","DRSP","ETG","LNGP","DRSPP"))
thr_tbl <- thr_tbl[order(thr_tbl$drug), ]; thr_tbl$drug <- as.character(thr_tbl$drug)
write.csv(thr_tbl, "results_R/tables/threshold_robustness.csv", row.names = FALSE)

# ---- Validation ----
say("validation ...")
pk_tbl <- pk_validation_table(DRUGS)
pu_tbl <- cache("pu_tbl", function() perfect_use_ovulation(get_cl(), DRUGS, N_MAIN))
miss_tbl <- cache("miss_tbl", function() missed_pill_validation(get_cl(), DRUGS, N_MAIN))
write.csv(pk_tbl, "results_R/tables/validation_pk.csv", row.names = FALSE)
write.csv(pu_tbl, "results_R/tables/validation_pd.csv", row.names = FALSE)
write.csv(miss_tbl, "results_R/tables/validation_missed_pill.csv", row.names = FALSE)

# Delayed-intake INDEPENDENT validation: the POP late-dose
# arm vs the pivotal scheduled-delay trials (Korver 2005 desogestrel 12 h;
# Duijkers 2016 drospirenone-only 24 h). Reuses the already-computed late grid.
di_tbl <- delayed_intake_validation(escape_grid[escape_grid$analysis == "late", ])
write.csv(di_tbl, "results_R/tables/validation_delayed_intake.csv", row.names = FALSE)
say("delayed-intake validation written (ordering_24h_ok=", isTRUE(attr(di_tbl, "ordering_24h_ok")), ")")

# ---- Sensitivity: local tornado + Morris global (checkpointed per drug) ----
say("sensitivity local ...")
sens_local <- rbind(cache("sens_LNG", function() tornado_for(get_cl(), DRUGS, "LNG", N_SENS)),
                    cache("sens_ETG", function() tornado_for(get_cl(), DRUGS, "ETG", N_SENS)))
write.csv(sens_local, "results_R/tables/sensitivity_local.csv", row.names = FALSE)
say("sensitivity Morris global ...")
morris_tbl <- tryCatch(rbind(
  cache("morris_LNG", function() morris_global(get_cl(), DRUGS, "LNG", r = if (QUICK) 4 else 40, n_inner = if (QUICK) 20 else 60)),
  cache("morris_ETG", function() morris_global(get_cl(), DRUGS, "ETG", r = if (QUICK) 4 else 40, n_inner = if (QUICK) 20 else 60))),
  error = function(e) { say("morris error: ", conditionMessage(e)); NULL })
if (!is.null(morris_tbl)) write.csv(morris_tbl, "results_R/tables/sensitivity_morris.csv", row.names = FALSE)

# ---- Recovery action ----
say("recovery action ...")
rec_tbl <- cache("rec_tbl", function() recovery_table(get_cl(), DRUGS, n_sub = N_REC))
write.csv(rec_tbl, "results_R/tables/recovery.csv", row.names = FALSE)

# ---- OR -> Pearl Index translation (main grid) ----
say("OR->Pearl-Index translation ...")
pearl_tbl <- cache("pearl_tbl", function() {
  maingrid <- escape_grid[escape_grid$analysis == "main",
                          c("drug", "position", "n_missed", "escape_prob")]
  pearl_table(DRUGS, maingrid)
})
write.csv(pearl_tbl, "results_R/tables/pearl_index.csv", row.names = FALSE)

if (!is.null(.cl)) stopCluster(.cl)

# ---- Qualification traces (single typical subject; for figures) ----
say("qualification traces ...")
qual <- cache("qual", function() {
  q <- list()
  tg <- seq(0, 90, by = 0.05)
  profU <- list(t = tg, c_prog = numeric(length(tg)), c_ee = numeric(length(tg)))
  q$untreated <- hpo_simulate(profU, 90, 1e9, 0, 2, NA, 0, 2)
  schedL <- apply_adherence(build_pack_schedule(21, 7, 4))
  profL <- make_profile(DRUGS$LNG, schedL, 112)
  q$perfect_lng <- hpo_simulate(profL, 112, ec50_prog(DRUGS$LNG), DRUGS$LNG$emax_prog,
                                DRUGS$LNG$hill_prog, ec50_ee(DRUGS$LNG), DRUGS$LNG$ee_emax, DRUGS$LNG$ee_hill)
  q$escape_example <- run_scenario(DRUGS$LNG, n_missed = 6, position = "start",
                                   challenge_pack = 2, follow_days = 34, return_trace = TRUE)
  q
})

# accumulate runtime across resumed invocations
prev_rt <- if (file.exists(RT_FILE)) readRDS(RT_FILE) else 0
runtime <- prev_rt + as.numeric(Sys.time() - t_start, units = "secs")
meta <- list(date = format(Sys.time(), "%Y-%m-%d %H:%M"), n_main = N_MAIN,
             runtime_s = round(runtime), quick = QUICK,
             n_cells = nrow(escape_grid), accum_note = pk_metrics_note(DRUGS))

saveRDS(list(escape_grid = escape_grid, mt = mt, gc = gc, fi = fi, thr_tbl = thr_tbl,
             pk_tbl = pk_tbl, pu_tbl = pu_tbl, miss_tbl = miss_tbl, di_tbl = di_tbl,
             h2_paired = h2_paired, sens_local = sens_local,
             morris_tbl = morris_tbl, rec_tbl = rec_tbl, pearl_tbl = pearl_tbl,
             qual = qual, meta = meta, drugs_tbl = NULL),
        "results_R/bundle.rds")
writeLines(capture.output(sessionInfo()), "results_R/sessionInfo.txt")

# success: clear checkpoints so a future fresh run recomputes
unlink(file.path(CKPT, "*.rds")); unlink(CKPT, recursive = TRUE)
say("DONE in ", round(runtime), "s (cumulative). bundle.rds written; checkpoints cleared.")
cat(sprintf("DONE in %ds (cumulative). Outputs in results_R/.\n", round(runtime)))
print(mt[, c("drug", "position", "max_tolerated", "mt_lo", "mt_hi")])
