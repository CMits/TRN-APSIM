# =============================================================================
# config.R — single source of truth for PSoup repository paths (R side).
#
# Every R script sources this file instead of hard-coding an absolute ROOT.
# Entry-point scripts (run via `Rscript pipeline/xx.R`) begin with a 3-line
# header that locates the repo root from the script's own path and sources
# this file:
#
#   .a <- commandArgs(FALSE); .f <- sub('^--file=', '', grep('^--file=', .a, value = TRUE))
#   ROOT <- if (length(.f)) normalizePath(file.path(dirname(.f[1]), '..'), winslash = '/') else normalizePath(getwd(), winslash = '/')
#   source(file.path(ROOT, 'config.R'))
#
# Scripts that are source()d by others (model/psoup_driver.R, model/randomisers.R,
# model/build_network.R) inherit ROOT + the *_DIR variables from the caller's
# environment, so they need no header. This file also self-locates when sourced
# directly, and honours the PSOUP_ROOT environment variable if set.
# =============================================================================

if (!exists("ROOT") || is.null(ROOT)) {
  env_root <- Sys.getenv("PSOUP_ROOT", unset = NA_character_)
  if (!is.na(env_root) && nzchar(env_root)) {
    ROOT <- normalizePath(env_root, winslash = "/", mustWork = FALSE)
  } else {
    # Self-locate: config.R lives at the repository root.
    of <- NULL
    for (i in seq_len(sys.nframe())) {
      o <- sys.frame(i)$ofile
      if (!is.null(o)) of <- o
    }
    ROOT <- if (!is.null(of)) normalizePath(dirname(of), winslash = "/", mustWork = FALSE)
            else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
}

DATA_DIR           <- file.path(ROOT, "data")
MODEL_DIR          <- file.path(ROOT, "model")
MODEL_BASELINE_DIR <- file.path(MODEL_DIR, "model_baseline")
PIPELINE_DIR       <- file.path(ROOT, "pipeline")
INPUTS_DIR         <- file.path(PIPELINE_DIR, "inputs")
OUT_DIR            <- file.path(PIPELINE_DIR, "out")
FIG_DIR            <- file.path(ROOT, "figures")
DOCS_DIR           <- file.path(ROOT, "docs")

# Ensure output dirs exist for entry scripts that write results.
for (.d in c(OUT_DIR, INPUTS_DIR)) if (!dir.exists(.d)) dir.create(.d, recursive = TRUE, showWarnings = FALSE)

# Resolve a Supplementary Dataset file by S-number + version (long filenames).
dataset_path <- function(number, version = "V3") {
  prefix <- sprintf("Dataset S%d_", number)
  files  <- list.files(DATA_DIR, pattern = "\\.xlsx$", full.names = TRUE)
  hits   <- files[startsWith(basename(files), prefix) &
                    grepl(sprintf("_%s\\.xlsx$", version), basename(files))]
  if (length(hits) == 1) return(hits)
  alt <- files[startsWith(basename(files), prefix)]
  if (length(alt) == 1) return(alt)
  stop(sprintf("dataset_path: could not uniquely resolve S%d (%s) in %s", number, version, DATA_DIR))
}
