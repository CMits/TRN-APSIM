# =============================================================================
# 00_build_network.R -- single source of truth for the corrected 32-node network
#
# Loads Dataset S18 (SBGN), applies the two corrections below, and returns a
# PSoup Network object. Sourced by the baseline scorer, the randomiser tests and
# the null runner, so all three operate on an identical starting topology.
#
# The two corrections are different in kind and are kept separate on purpose:
#
#   (1) PARSER PATCH -- restores what the SBGN already says.
#       convertSBGNdiagram() misreads exactly one arc of 56. The file has
#         <arc class="positive influence" source=SUC target=Sustained_growth>
#       but the parsed object records "inhibition". Verified by diffing every
#       arc in the raw XML against the parsed inputs tables: this is the only
#       discrepancy, and no negative arc into Sustained_growth exists anywhere
#       in the file. Left unpatched it inverts sucrose's role:
#         wrong : 2*min(GA,PAT_Bud)/(1+SUC) * Bud_release
#         right : ((SUC + min(GA,PAT_Bud))/2) * Bud_release
#
#   (2) BIOLOGICAL CORRECTION -- changes what the SBGN says, because the
#       diagram is wrong. SL acts on CK indirectly, via
#       Perception_SL -| SMXL6_7_8 -> CK, not through a direct
#       Perception_SL -| CK edge.
#       Yields  CK = 2*((SUC + SMXL6_7_8)/2)/(1 + Aux_Shoot) * CK_n.
#
# With both applied, generateEdgeList() reproduces Network/Edges.csv exactly,
# and the regenerated equations match model/psoup_simulator.py to
# 1.6e-13 across 203 genotype vectors (see 03_compare_equations.py).
# =============================================================================

if (!requireNamespace("PSoup", quietly = TRUE)) stop("PSoup is not installed")

# Repo-relative paths. When sourced by an entry script, DATA_DIR is already in
# scope (from config.R); when run standalone, self-locate and source config.
if (!exists("DATA_DIR")) {
  .a <- commandArgs(FALSE); .f <- sub("^--file=", "", grep("^--file=", .a, value = TRUE))
  .root <- if (length(.f)) normalizePath(file.path(dirname(.f[1]), ".."), winslash = "/") else normalizePath(getwd(), winslash = "/")
  source(file.path(.root, "config.R"))
}
SBGN_PATH <- file.path(DATA_DIR, "Dataset S18_Benchmark_Model_Original.sbgn")

PARSER_PATCHES <- list(
  list(node = "Sustained_growth", from = "SUC",
       was = "inhibition", should_be = "stimulation")
)

BIOLOGICAL_EDITS <- list(
  list(action = "drop", from = "Perception_SL", to = "CK"),
  list(action = "add",  from = "SMXL6_7_8",     to = "CK", influence = "stimulation")
)


.drop_edge <- function(net, from, to) {
  ti <- net@objects$Hormones[[to]]@inputs
  k <- which(ti$Node == from)
  if (length(k) != 1)
    stop(sprintf("drop_edge: %s -> %s matched %d input rows (expected 1)", from, to, length(k)))
  net@objects$Hormones[[to]]@inputs <- ti[-k, , drop = FALSE]
  fo <- net@objects$Hormones[[from]]@outputs
  k2 <- which(fo$Node == to)
  if (length(k2) == 1) net@objects$Hormones[[from]]@outputs <- fo[-k2, , drop = FALSE]
  net
}

.add_edge <- function(net, from, to, influence) {
  ti <- net@objects$Hormones[[to]]@inputs
  if (any(ti$Node == from)) stop(sprintf("add_edge: %s -> %s already present", from, to))
  net@objects$Hormones[[to]]@inputs <- rbind(
    ti, data.frame(Node = from, Coregulator = NA_character_, Influence = influence,
                   Operator = NA_character_, stringsAsFactors = FALSE))
  fo <- net@objects$Hormones[[from]]@outputs
  net@objects$Hormones[[from]]@outputs <- rbind(
    fo, data.frame(Node = to, Coregulator = NA_character_, Influence = influence,
                   Operator = NA_character_, stringsAsFactors = FALSE))
  net
}


#' Load Dataset S18 and return the corrected Network object.
#'
#' Every patch asserts its precondition first. If PSoup's parsing behaviour ever
#' changes, this stops loudly rather than silently producing a different model.
build_corrected_network <- function(sbgn = SBGN_PATH, name = "shootBranching",
                                    quiet = FALSE) {
  stopifnot(file.exists(sbgn))
  net <- PSoup::convertSBGNdiagram(file = sbgn, networkName = name)

  for (p in PARSER_PATCHES) {
    inp <- net@objects$Hormones[[p$node]]@inputs
    idx <- which(inp$Node == p$from)
    if (length(idx) != 1)
      stop(sprintf("parser patch %s -> %s matched %d rows (expected 1)",
                   p$from, p$node, length(idx)))
    if (inp$Influence[idx] != p$was)
      stop(sprintf(paste0("parser patch precondition failed for %s -> %s: expected '%s', ",
                          "found '%s'. PSoup behaviour has changed -- re-verify against ",
                          "the SBGN before proceeding."),
                   p$from, p$node, p$was, inp$Influence[idx]))
    net@objects$Hormones[[p$node]]@inputs$Influence[idx] <- p$should_be
    if (!quiet) message(sprintf("  parser patch : %s -> %s  '%s' => '%s'",
                                p$from, p$node, p$was, p$should_be))
  }

  for (e in BIOLOGICAL_EDITS) {
    if (e$action == "drop") {
      net <- .drop_edge(net, e$from, e$to)
      if (!quiet) message(sprintf("  bio edit     : removed %s -> %s", e$from, e$to))
    } else {
      net <- .add_edge(net, e$from, e$to, e$influence)
      if (!quiet) message(sprintf("  bio edit     : added %s -> %s [%s]", e$from, e$to, e$influence))
    }
  }

  net
}
