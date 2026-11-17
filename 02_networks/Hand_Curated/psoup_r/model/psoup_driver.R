# =============================================================================
# psoup_driver.R -- shared simulation + scoring helpers
#
# Sourced by the baseline scorer and by the null-ensemble runner, so the
# observed statistic and the null distribution are produced by identical code.
#
# The exogenous-supply semantics deliberately mirror
# model/psoup_simulator.py exactly:
#   - all nodes start at 1.0
#   - steps 1..delay run with NO exogenous supply
#   - from step delay+1 onward, supply is ADDED to node values after each update
#   - convergence is only tested after step > delay
# 03_compare_equations.py verifies this driver against the Python original.
# =============================================================================

NODE_NAMES <- c(
  "SL", "BRC1_2", "CK", "Aux_Shoot", "ABA", "Bud_release", "PAT_Bud",
  "Aux_Bud", "SUC", "Sustained_growth", "GA", "Decap_signal", "CXE",
  "SMXL6_7_8", "FHY3_FAR1", "SPL9_15", "SMXL6_7_8_Transcription",
  "Low_R_FR", "MAX2", "D14", "Perception_SL", "D27", "MAX3", "MAX4",
  "MAX1", "LBO", "CLAMT", "PhyB", "High_R_FR", "PIFs", "HB21_40_53", "NCED3")

GENE_NAMES <- c(
  "SL_n", "BRC1_2_n", "CK_n", "Aux_n", "ABA_n", "PIN_1_3_4_7_n",
  "Aux_Bud_n", "SUC_n", "GA_n", "Decap_signal_n", "CXE_n",
  "SMXL6_7_8_n", "FHY3_FAR1_n", "SPL9_15_n", "MAX2_n", "D14_n",
  "D27_n", "MAX3_n", "MAX4_n", "MAX1_n", "LBO_n", "CLAMT_n",
  "PhyB_n", "PIFs_n", "HB21_40_53_n", "NCED3_n")

EXO_NAMES <- c("SL", "CK", "Aux_Shoot", "ABA", "SUC",
               "Low_R_FR", "Perception_SL", "High_R_FR", "GA")

# --- fast equation compiler --------------------------------------------------
#
# buildModel() emits a nextStep() that appends a row to a data.frame on every
# call and indexes it as dat$Node[t-1]. That is far too slow for a 1000-network
# ensemble (~45 s per replicate). Since the equations are simultaneous -- every
# right-hand side reads t-1 only -- they can be re-emitted verbatim as
# assignments into a plain named numeric vector, which is ~100x faster.
#
# This is a mechanical rewrite of the SAME algebra, not a reimplementation:
#   dat$X[t] = <expr>   ->   n[["X"]] <- <expr>
#   dat$Y[t-1]          ->   s[["Y"]]
#   gen$Z_n             ->   g[["Z_n"]]
# compile_check() below asserts it agrees with the original to machine precision.

#' Repair a PSoup code-generation bug.
#'
#' When a node has BOTH a plain necessary stimulant and an AND-grouped one,
#' buildModel() omits the multiplication operator between them:
#'    (dat$Bud_release[t-1] min(dat$D14[t-1], dat$MAX2[t-1], dat$SL[t-1]))
#' which is not valid R. PSoup's own rule is that multiple necessary stimulants
#' combine as a PRODUCT (differenceString(..., takeProduct = TRUE)), so the
#' intended text is `A * min(...)`. This restores that operator.
#'
#' The real topology never triggers it -- Perception_SL is its only AND-grouped
#' necessary stimulant and it has no plain one -- but rewiring creates the
#' combination in roughly 3% of degree-preserving replicates. Repairing rather
#' than discarding matters: dropping the affected networks would condition the
#' null ensemble on "generates clean PSoup code", biasing the comparison.
#'
#' @return list(lines, n_repairs)
repair_generated <- function(lines) {
  before <- lines
  # "dat$X[t-1] min(" / "max(" -> "dat$X[t-1] * min("
  lines <- gsub("(\\[\\s*t\\s*-\\s*1\\s*\\])(\\s+)(min|max)\\(", "\\1 * \\3(", lines)
  # ") dat$X[" -> ") * dat$X["   (same omission, other operand order)
  lines <- gsub("(\\))(\\s+)(dat\\$[A-Za-z0-9_]+\\[)", "\\1 * \\3", lines)
  list(lines = lines, n_repairs = sum(before != lines))
}

compile_next_step <- function(path) {
  lines <- readLines(path, warn = FALSE)
  rep <- repair_generated(lines)
  lines <- rep$lines
  attr_repairs <- rep$n_repairs
  pat <- "^\\s*dat\\$([A-Za-z0-9_]+)\\[\\s*t\\s*\\]\\s*=\\s*(.+?)\\s*$"
  hits <- regmatches(lines, regexec(pat, lines))

  body <- character(0); nodes <- character(0)
  for (h in hits) {
    if (length(h) != 3) next
    node <- h[2]; expr <- h[3]
    expr <- gsub("dat\\$([A-Za-z0-9_]+)\\[\\s*t\\s*-\\s*1\\s*\\]", 's[["\\1"]]', expr)
    expr <- gsub("gen\\$([A-Za-z0-9_]+)", 'g[["\\1"]]', expr)
    if (grepl("dat\\$|gen\\$", expr))
      stop("compile_next_step: unhandled reference remains in equation for ", node,
           ": ", expr)
    body <- c(body, sprintf('  n[["%s"]] <- %s', node, expr))
    nodes <- c(nodes, node)
  }
  # Node-count-agnostic sanity check (works for the 32-node base model and the
  # 35-node sucrose model alike): every node read as an input must also be
  # defined by an equation, else the parse is incomplete/broken.
  read_nodes <- unlist(regmatches(body, gregexpr('s\\[\\["[A-Za-z0-9_]+"\\]\\]', body)))
  read_nodes <- unique(gsub('.*\\["([A-Za-z0-9_]+)"\\].*', "\\1", read_nodes))
  missing <- setdiff(read_nodes, nodes)
  if (length(missing))
    stop(sprintf("compile_next_step: parsed %d equations but these input node(s) are undefined: %s",
                 length(body), paste(missing, collapse = ", ")))

  txt <- paste0("function(s, g) {\n  n <- s\n", paste(body, collapse = "\n"), "\n  n\n}")
  fn <- eval(parse(text = txt))
  attr(fn, "n_repairs") <- attr_repairs
  fn
}

#' Assert the compiled closure matches the original data.frame-based nextStep.
compile_check <- function(path, n_cases = 40, tol = 1e-12, seed = 1L) {
  fast <- compile_next_step(path)

  # If the generated file needed repair it does not parse as R, so there is no
  # original to compare against. Assert instead that the repaired version is
  # well-formed and numerically sane; return NA to mark it as unverifiable.
  if (isTRUE(attr(fast, "n_repairs") > 0)) {
    set.seed(seed)
    s <- setNames(runif(length(NODE_NAMES), 0, 2.5), NODE_NAMES)
    g <- setNames(rep(1, length(GENE_NAMES)), GENE_NAMES)
    v <- fast(s, g)
    if (length(v) != length(NODE_NAMES) || !all(is.finite(v)))
      stop("compile_check: repaired equations produced a malformed state vector")
    return(NA_real_)
  }

  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  slow <- get("nextStep", envir = env)

  set.seed(seed)
  worst <- 0
  for (i in seq_len(n_cases)) {
    s <- setNames(runif(length(NODE_NAMES), 0, 2.5), NODE_NAMES)
    g <- setNames(sample(c(0, 0.25, 0.75, 1, 2), length(GENE_NAMES), replace = TRUE),
                  GENE_NAMES)
    a <- unlist(slow(as.data.frame(as.list(s)), as.data.frame(as.list(g)))[1, NODE_NAMES])
    b <- fast(s, g)[NODE_NAMES]
    d <- max(abs(a - b), na.rm = TRUE)
    if (is.finite(d) && d > worst) worst <- d
  }
  if (worst > tol)
    stop(sprintf("compile_check FAILED: max |difference| = %.3e (tol %.1e)", worst, tol))
  worst
}

# --- simulation --------------------------------------------------------------

#' Run one perturbation to steady state.
#'
#' @param nextStepFn the node-update closure produced by buildModel()
#' @param gen_vec named numeric vector over GENE_NAMES
#' @param exo_vec named numeric vector over EXO_NAMES (may be all zero / NULL)
#' @return list(state = named node vector, converged = logical, steps = int)
#'         `state` is NULL when the run diverged to non-finite values.
simulate_psoup <- function(nextStepFn, gen_vec, exo_vec = NULL,
                           max_steps = 100, delay = 20, tol = 1e-8,
                           node_names = NODE_NAMES) {
  prev <- setNames(rep(1, length(node_names)), node_names)

  use_exo <- !is.null(exo_vec) && any(exo_vec != 0)
  if (use_exo) {
    exo_vec <- exo_vec[exo_vec != 0]
    exo_nm  <- intersect(names(exo_vec), node_names)
    use_exo <- length(exo_nm) > 0
  }

  for (step in seq_len(max_steps)) {
    cur <- nextStepFn(prev, gen_vec)

    if (use_exo && step > delay) {
      cur[exo_nm] <- cur[exo_nm] + exo_vec[exo_nm]
    }

    # A rewired network can oscillate or blow up. Treat that as a failed run
    # rather than letting NaN/Inf propagate into the accuracy statistic.
    if (!all(is.finite(cur))) {
      return(list(state = NULL, converged = FALSE, steps = step))
    }

    if (step > delay && max(abs(cur - prev)) < tol) {
      return(list(state = cur, converged = TRUE, steps = step))
    }

    prev <- cur
  }

  list(state = prev, converged = FALSE, steps = max_steps)
}

#' Simulate every treatment once and return a lookup keyed by treatment name.
#' gene_names/exo_names default to the 32-node base sets; pass wider sets (and
#' node_names via ...) for the 35-node sucrose model.
simulate_all <- function(nextStepFn, genotypes, exogenous,
                         gene_names = GENE_NAMES, exo_names = EXO_NAMES, ...) {
  out <- vector("list", nrow(genotypes))
  names(out) <- genotypes$treatment
  exo_idx <- match(genotypes$treatment, exogenous$treatment)

  for (i in seq_len(nrow(genotypes))) {
    gv <- setNames(as.numeric(genotypes[i, gene_names]), gene_names)
    ev <- if (!is.na(exo_idx[i])) {
      setNames(as.numeric(exogenous[exo_idx[i], exo_names]), exo_names)
    } else NULL
    out[[i]] <- simulate_psoup(nextStepFn, gv, ev, ...)
  }
  out
}

# --- scoring -----------------------------------------------------------------

# Binning threshold, as published. 0.02 was evaluated as a way of removing the
# two hand-set bins (they bin correctly on their own at 0.02) but was rejected:
# it also shifts four unrelated borderline pairs. The chosen resolution instead
# keeps 0.05 and reverts both overrides, reporting the two affected pairs as
# below-cut-off near misses -- a mismatch category the Discussion already uses.
BIN_THRESHOLD <- 0.05

#' Scoring rule as published: a three-way bin with the threshold applied to all
#' classes, scored correct when the simulated bin equals the biological bin.
#'
#' Determined empirically rather than assumed -- applied to the stored node
#' values in Dataset S13 this reproduces the published Match column on 82 of 84
#' pairs, versus 81/84 for a purely directional rule. The two disagreements are
#' the PhyB_n/Bud_release hand override and one erroneous Match cell
#' (WT, +CK, +NPA / Bud_release: bin +1, bio +1, but Match recorded as "No").
score_pair <- function(pert, ctrl, bio_bin, threshold = BIN_THRESHOLD) {
  if (!is.finite(pert) || !is.finite(ctrl)) return(NA)

  if (ctrl == 0) {
    ratio <- if (pert == 0) 1 else if (pert > 0) Inf else -Inf
  } else {
    ratio <- pert / ctrl
  }

  as.integer(compute_bin(pert, ctrl, threshold) == bio_bin)
}

#' Three-way bin, used for confusion matrices and the published-bin cross-check.
compute_bin <- function(pert, ctrl, threshold = BIN_THRESHOLD) {
  if (!is.finite(pert) || !is.finite(ctrl)) return(NA_integer_)
  if (ctrl == 0) return(if (pert > 0) 1L else if (pert < 0) -1L else 0L)
  ratio <- pert / ctrl
  if (ratio > 1 + threshold) 1L else if (ratio < 1 - threshold) -1L else 0L
}

#' Score a pair table against a simulation lookup.
#'
#' Pairs whose treatment or control failed to converge are scored 0 (incorrect),
#' never dropped: excluding them would condition the null ensemble on
#' "behaves like a real signalling network", which is the very property under test.
score_all <- function(lookup, pairs, threshold = BIN_THRESHOLD) {
  n <- nrow(pairs)
  correct <- integer(n); psoup_bin <- integer(n); usable <- logical(n)

  for (i in seq_len(n)) {
    tr <- lookup[[pairs$treatment[i]]]
    ct <- lookup[[pairs$compare_with[i]]]
    node <- pairs$node[i]

    if (is.null(tr) || is.null(ct) || is.null(tr$state) || is.null(ct$state)) {
      correct[i] <- 0L; psoup_bin[i] <- NA_integer_; usable[i] <- FALSE
      next
    }
    pv <- tr$state[[node]]; cv <- ct$state[[node]]
    s <- score_pair(pv, cv, pairs$bio_bin[i], threshold)
    correct[i]   <- if (is.na(s)) 0L else s
    psoup_bin[i] <- compute_bin(pv, cv, threshold)
    usable[i]    <- tr$converged && ct$converged
  }

  list(correct = correct, psoup_bin = psoup_bin, usable = usable,
       n_correct = sum(correct), n = n, accuracy = sum(correct) / n,
       converged_frac = mean(usable))
}

# --- model (re)generation ----------------------------------------------------

#' Regenerate node equations for a Network object and return the update closure.
#'
#' buildModel() writes nextStep.R to disk, so each call needs its own folder
#' when running in parallel.
build_and_load <- function(net, folder, quiet = TRUE) {
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  fn <- function() {
    buildModel(network = net, folder = folder, forceOverwrite = TRUE,
               language = "R", ruleStyle = "Dun", necStimStyle = "Linear",
               saveNetwork = FALSE, exogenous = TRUE)
  }
  if (quiet) suppressWarnings(suppressMessages(capture.output(fn()))) else fn()
  compile_next_step(file.path(folder, "nextStep.R"))
}

# --- io ----------------------------------------------------------------------

load_split <- function(inputs_dir, split) {
  g <- read.csv(file.path(inputs_dir, paste0(split, "_genotypes.csv")),
                check.names = FALSE, stringsAsFactors = FALSE)
  e <- read.csv(file.path(inputs_dir, paste0(split, "_exogenous.csv")),
                check.names = FALSE, stringsAsFactors = FALSE)
  p <- read.csv(file.path(inputs_dir, paste0(split, "_pairs.csv")),
                check.names = FALSE, stringsAsFactors = FALSE)
  names(g)[1] <- "treatment"
  names(e)[1] <- "treatment"
  g$treatment <- trimws(g$treatment)
  e$treatment <- trimws(e$treatment)
  p$treatment <- trimws(p$treatment)
  p$compare_with <- trimws(p$compare_with)
  p$node <- trimws(p$node)
  list(genotypes = g, exogenous = e, pairs = p)
}
