"""The paper's `_n` modifier columns mapped onto FlashP network node ids.

Lives in its own module so that consumers which only need the MAP (e.g.
build_apsim_parameters.py) do not have to import run_sucrose_sweep.py, which
pulls in the FlashP validator at module level. Importing the map must not
require a simulator.

Networks do not share a node-naming convention, so the map is keyed by VARIANT:

  benchmark    Tillering_Sorghum_Benchmark(_suc) -- benchmark-style labels
  psoup_style  Shoot_Branching_Sorghum -- PSoup-style labels, as in Hand_Curated

Each arm names its variant in scripts/arms.json.
"""

from __future__ import annotations

# paper `_n` modifier column -> FlashP node id.
# Node labels in this network are benchmark-style by design; the underlying
# sorghum orthologues are recorded in the network's ortholog_crosswalk.json.
_BENCHMARK = {
    "BRC1_2_n": "BRC1",
    "D14_n": "D14",
    "HB21_40_53_n": "GT1",          # HD-Zip I downstream of BRC1/TB1
    "MAX1_n": "MAX1",
    "MAX2_n": "MAX2",
    "MAX3_n": "MAX3",
    "MAX4_n": "MAX4",
    "PIN_1_3_4_7_n": "PIN1",
    "SMXL6_7_8_n": "SMXL678",
    "CLAMT_n": "CLAMT",
    "CXE_n": "CXE",
    "D27_n": "D27",
    "FHY3_FAR1_n": "FHY3FAR1",
    "LBO_n": "LBO",
    "NCED3_n": "NCED3",
    "SPL9_15_n": "SPL9_15",
    "SUC_n": "Sucrose",
}

# Shoot_Branching_Sorghum uses PSoup-style labels (as Hand_Curated does), so four
# entries differ. The other thirteen resolve unchanged.
#
# HB21_40_53_n stays mapped to GT1 even though this network has a native
# HB21_40_53 node: keeping the mapping identical to the benchmark arm means any
# difference between the two FlashP arms is attributable to topology rather than
# to a changed gene map. Measured both ways, the choice moves that genotype
# between rank 11 and 12 of 32 and leaves the top-5/bottom-5 sets identical.
_PSOUP_STYLE = {
    **_BENCHMARK,
    "BRC1_2_n": "BRC1_2",
    "SMXL6_7_8_n": "SMXL6_7_8",
    "FHY3_FAR1_n": "FHY3_FAR1",
    "SUC_n": "SUC",
}

VARIANTS = {"benchmark": _BENCHMARK, "psoup_style": _PSOUP_STYLE}

# paper modifiers with no counterpart node in this network — these act on nodes the
# FlashP model does not represent separately. Recorded so the omission is explicit.
#
# Left identical across variants on purpose. Shoot_Branching_Sorghum does contain
# SL/CK/ABA/GA/Decap_signal/PHYB/PIFS/Aux_Bud/Aux_Shoot, but none of these nine
# modifiers is ever perturbed away from 1.0 by any of the paper's 37 genotypes, so
# mapping them is a numerical no-op; the casing also differs (PHYB vs PhyB_n) and a
# wrong guess would be silent. run_sucrose_sweep.py additionally requires NO_NODE to
# stay a superset of the unmapped `_n` columns or it refuses to run.
_NO_NODE = {"SL_n", "CK_n", "Aux_n", "ABA_n", "Aux_Bud_n", "GA_n",
            "Decap_signal_n", "PhyB_n", "PIFs_n"}

NO_NODE_BY_VARIANT = {k: set(_NO_NODE) for k in VARIANTS}


def gene_to_node(variant: str = "benchmark") -> dict[str, str]:
    try:
        return dict(VARIANTS[variant])
    except KeyError:
        raise SystemExit(
            f"unknown gene-map variant {variant!r} -- "
            f"choose from {', '.join(VARIANTS)}") from None


def no_node(variant: str = "benchmark") -> set[str]:
    gene_to_node(variant)                      # validate the key
    return set(NO_NODE_BY_VARIANT[variant])


# Back-compatible aliases bound to the benchmark variant, so importers that have not
# yet been made arm-aware keep their exact previous behaviour.
GENE_TO_NODE = _BENCHMARK
NO_NODE = _NO_NODE
