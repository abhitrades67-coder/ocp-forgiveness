# Standalone verification of two code additions:
#  (A) delayed_intake_validation()  -- from the existing full-run late grid
#  (B) paired_position_bootstrap()  -- H2, regenerating LNG/DRSP position matrices
# Run:  Rscript R/_verify_review_additions.R
suppressMessages({
  source("R/pk.R"); source("R/drugs.R"); source("R/hpo.R"); source("R/simulate.R")
  source("R/forgiveness.R"); source("R/parallel_util.R"); source("R/validate.R")
})
DRUGS <- build_drugs(); NS <- 0:7

## (A) Delayed-intake external validation from the existing grid ---------------
grid <- read.csv("results_R/tables/escape_probability_grid.csv", stringsAsFactors = FALSE)
late <- grid[grid$analysis == "late", ]
div <- delayed_intake_validation(late)
cat("=== (A) Delayed-intake external validation (Korver 2005 / Duijkers 2016) ===\n")
print(div, row.names = FALSE)
cat("cross-agent ordering at 24 h (DRSPP < ETG) ok:", attr(div, "ordering_24h_ok"), "\n\n")
write.csv(div, "results_R/tables/validation_delayed_intake.csv", row.names = FALSE)

## (B) Paired H2 bootstrap (same shared population across positions) -----------
N <- 600
cl <- make_model_cluster()
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
h2_all <- list()
for (dk in c("LNG", "DRSP")) {
  subs <- sample_population(N, drug = DRUGS[[dk]], tag = "main")
  mats <- list()
  for (pos in c("start", "mid", "late")) {
    mat <- sapply(NS, function(n) escape_prob_par(cl, dk, n, pos, subs, list(challenge_pack = 2)))
    if (is.null(dim(mat))) mat <- matrix(mat, nrow = N)
    mats[[pos]] <- mat
  }
  pb <- paired_position_bootstrap(mats, NS)
  cat(sprintf("=== (B) H2 paired bootstrap: %s (N=%d) ===\n", dk, N))
  cat("continuous forgiveness index (missed-pills-to-5%% escape), per position:\n")
  print(round(pb$point, 3))
  print(pb$contrasts, row.names = FALSE)
  cat("(diff>0 => first position MORE forgiving; excludes_zero => 95% CI clears 0)\n\n")
  df <- pb$contrasts; df$drug <- dk
  h2_all[[dk]] <- df
}
stopCluster(cl)
write.csv(do.call(rbind, h2_all), "results_R/tables/h2_paired.csv", row.names = FALSE)
cat("Wrote results_R/tables/validation_delayed_intake.csv and h2_paired.csv\n")
