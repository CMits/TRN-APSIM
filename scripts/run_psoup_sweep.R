# =============================================================================
# Step 2 (arm A) -- sucrose sweep through the HAND-CURATED PSoup network, in R.
#
# The mirror image of run_sucrose_sweep.py, which does the same job for the
# FlashP arm. Together they make the two arms symmetric: each one propagates the
# paper's genotypes through its OWN network at sucrose {0.5, 1.0, 1.5} and
# records the phenotype node.
#
#   arm A (here)              02_networks/Hand_Curated       -> Sustained_growth
#   arm B (run_sucrose_sweep) Tillering_Sorghum_Benchmark_suc -> Tiller_Number
#
# WHY R. The curated model is canonical in R: PSoup's buildModel() emits
# model_baseline/nextStep.R from Dataset S18_Benchmark_Model_Original.sbgn, and
# those equations carry min() AND-gates and multiplicative necessary-stimulant
# chains. This script runs those equations directly rather than a re-implementation.
#
# This run IS the hand-curated result. It is not an attempt to reproduce the
# published DataSetS1 values -- those came from an earlier R run that no longer
# exists. handcurated_vs_published.csv records the difference so it is visible,
# but nothing gates on it.
#
# Genotypes are read out of DataSetS1_Results_APSIM_Suc.csv, which carries all 26
# `_n` modifier columns per row, so no genotype list is hand-written here.
#
# Writes 01_audit/handcurated_sucrose_response.csv
#        01_audit/handcurated_vs_published.csv
# =============================================================================

ROOT   <- normalizePath(file.path(dirname(sub("^--file=", "",
            grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
PSOUP  <- file.path(ROOT, "02_networks", "Hand_Curated", "psoup_r")
SRC    <- file.path(ROOT, "00_reference", "paper_data", "DataSetS1_Results_APSIM_Suc.csv")
OUTDIR <- file.path(ROOT, "01_audit")

# compile_next_step() / simulate_psoup() come from PSoup's own driver -- the
# solver is not reimplemented here.
source(file.path(PSOUP, "model", "psoup_driver.R"))

# The 32 nodes and 26 modifiers of the curated model, in the published order.
NODES <- c("SL","BRC1_2","CK","Aux_Shoot","ABA","Bud_release","PAT_Bud","Aux_Bud",
           "SUC","Sustained_growth","GA","Decap_signal","CXE","SMXL6_7_8","FHY3_FAR1",
           "SPL9_15","SMXL6_7_8_Transcription","Low_R_FR","MAX2","D14","Perception_SL",
           "D27","MAX3","MAX4","MAX1","LBO","CLAMT","PhyB","High_R_FR","PIFs",
           "HB21_40_53","NCED3")
GENES <- c("SL_n","BRC1_2_n","CK_n","Aux_n","ABA_n","PIN_1_3_4_7_n","Aux_Bud_n","SUC_n",
           "GA_n","Decap_signal_n","CXE_n","SMXL6_7_8_n","FHY3_FAR1_n","SPL9_15_n",
           "MAX2_n","D14_n","D27_n","MAX3_n","MAX4_n","MAX1_n","LBO_n","CLAMT_n",
           "PhyB_n","PIFs_n","HB21_40_53_n","NCED3_n")

PHENOTYPE <- "Sustained_growth"

main <- function() {
  eq <- file.path(PSOUP, "model", "model_baseline", "nextStep.R")
  fn <- compile_next_step(eq)
  cat(sprintf("network : %s\n", eq))
  cat(sprintf("compiled: %d nodes, %d generator repairs\n",
              length(NODES), attr(fn, "n_repairs")))

  d <- read.csv(SRC, check.names = FALSE)
  missing <- setdiff(GENES, names(d))
  if (length(missing)) stop("missing modifier columns: ", paste(missing, collapse = ", "))
  cat(sprintf("genotypes: %d rows, %d unique x %d sucrose levels\n",
              nrow(d), length(unique(d$genotype)), length(unique(d$sucrose_label))))

  n <- nrow(d)
  val <- rep(NA_real_, n); conv <- logical(n); steps <- integer(n)

  for (i in seq_len(n)) {
    g <- setNames(as.numeric(d[i, GENES]), GENES)   # sucrose enters as SUC_n
    r <- simulate_psoup(fn, g, NULL, node_names = NODES)
    if (!is.null(r$state)) val[i] <- r$state[[PHENOTYPE]]
    conv[i]  <- isTRUE(r$converged)
    steps[i] <- r$steps
  }

  dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

  out <- data.frame(
    genotype      = d$genotype,
    sucrose_label = tolower(trimws(d$sucrose_label)),
    SUC_n         = d$SUC_n,
    Sustained_growth = val,
    converged     = conv,
    steps         = steps,
    stringsAsFactors = FALSE)
  f1 <- file.path(OUTDIR, "handcurated_sucrose_response.csv")
  write.csv(out, f1, row.names = FALSE)

  pub <- as.numeric(d$Sustained_growth)
  cmp <- data.frame(
    genotype      = d$genotype,
    sucrose_label = out$sucrose_label,
    SUC_n         = d$SUC_n,
    published_Sustained_growth = pub,
    new_Sustained_growth       = val,
    abs_delta     = abs(val - pub),
    converged     = conv,
    stringsAsFactors = FALSE)
  f2 <- file.path(OUTDIR, "handcurated_vs_published.csv")
  write.csv(cmp, f2, row.names = FALSE)

  # ---- report ----
  nb <- sum(!conv)
  cat(sprintf("\nconverged : %d/%d\n", sum(conv), n))
  if (nb) {
    cat("not converged (values are finite and retained):\n")
    for (i in which(!conv))
      cat(sprintf("    %-34s %-7s -> %.4f\n", d$genotype[i], out$sucrose_label[i], val[i]))
  }
  rho <- suppressWarnings(cor(val, pub, method = "spearman", use = "complete.obs"))
  cat(sprintf("\nvs published DataSetS1 (recorded, not a gate):\n"))
  cat(sprintf("    rows within 1e-6 : %d/%d\n", sum(abs(val - pub) < 1e-6, na.rm = TRUE), n))
  cat(sprintf("    max |delta|      : %.4f\n", max(abs(val - pub), na.rm = TRUE)))
  cat(sprintf("    Spearman rho     : %.4f\n", rho))
  cat(sprintf("\nwrote %s\n      %s\n", f1, f2))
}

main()
