# Parallel helper (Windows PSOCK). Each worker sources the model files.
library(parallel)

make_model_cluster <- function(n_workers = max(1, detectCores() - 1),
                               src_dir = "R") {
  cl <- makeCluster(n_workers, type = "PSOCK")
  src_abs <- normalizePath(src_dir)
  clusterExport(cl, "src_abs", envir = environment())
  clusterEvalQ(cl, {
    suppressMessages({
      source(file.path(src_abs, "pk.R")); source(file.path(src_abs, "drugs.R"))
      source(file.path(src_abs, "hpo.R")); source(file.path(src_abs, "simulate.R"))
      source(file.path(src_abs, "forgiveness.R")); source(file.path(src_abs, "recovery.R"))
      source(file.path(src_abs, "translate.R"))
    })
    DRUGS <- build_drugs()
    TRUE
  })
  cl
}

# Per-subject escape vector over a population. The RK4 integrator is fixed-cost
# per sim (load is already balanced), so we PRESCHEDULE: split the population into
# one chunk per worker (length(cl) dispatches, not one per subject). This cut the
# main-grid time ~2x vs per-task parLapplyLB by removing dispatch overhead.
escape_prob_par <- function(cl, drug_key, n_missed, position, subjects,
                            scen_kwargs = list()) {
  nw <- length(cl)
  chunks <- split(subjects, (seq_along(subjects) - 1L) %% nw)
  res <- parLapply(cl, chunks, function(chunk, dk, nm, pos, sk) {
    drug <- DRUGS[[dk]]
    vapply(chunk, function(s) {
      args <- c(list(drug = drug, n_missed = nm, position = pos,
                     cl_factor = s$cl_factor, vd_factor = s$vd_factor,
                     ee_cl_factor = s$ee_cl_factor, pd = s$pd), sk)
      as.integer(isTRUE(do.call(run_scenario, args)$escaped))
    }, integer(1))
  }, dk = drug_key, nm = n_missed, pos = position, sk = scen_kwargs)
  # reassemble in original subject order
  out <- integer(length(subjects))
  idx <- split(seq_along(subjects), (seq_along(subjects) - 1L) %% nw)
  for (k in seq_along(res)) out[idx[[k]]] <- res[[k]]
  out
}
