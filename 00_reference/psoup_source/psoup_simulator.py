"""
PSoup Simulator — Python port of Script/nextStep.R (32-node version).

Equations are a LINE-BY-LINE port of the R implementation. No modifications,
no sucrose extension. This is the canonical 32-node model used to generate
the original S8 dataset.

The simulation iterates a discrete-time update rule until steady state.
Exogenous supply is applied additively to node values AFTER each step,
starting at step `delay+1`.
"""

from typing import Dict, Optional

# All 32 node names (output variables)
NODE_NAMES = [
    "SL", "BRC1_2", "CK", "Aux_Shoot", "ABA", "Bud_release", "PAT_Bud",
    "Aux_Bud", "SUC", "Sustained_growth", "GA", "Decap_signal", "CXE",
    "SMXL6_7_8", "FHY3_FAR1", "SPL9_15", "SMXL6_7_8_Transcription",
    "Low_R_FR", "MAX2", "D14", "Perception_SL", "D27", "MAX3", "MAX4",
    "MAX1", "LBO", "CLAMT", "PhyB", "High_R_FR", "PIFs", "HB21_40_53",
    "NCED3",
]

# 26 genotype modifiers (genes that can be knocked out/overexpressed)
GENE_NAMES = [
    "SL_n", "BRC1_2_n", "CK_n", "Aux_n", "ABA_n", "PIN_1_3_4_7_n",
    "Aux_Bud_n", "SUC_n", "GA_n", "Decap_signal_n", "CXE_n",
    "SMXL6_7_8_n", "FHY3_FAR1_n", "SPL9_15_n", "MAX2_n", "D14_n",
    "D27_n", "MAX3_n", "MAX4_n", "MAX1_n", "LBO_n", "CLAMT_n",
    "PhyB_n", "PIFs_n", "HB21_40_53_n", "NCED3_n",
]

# 9 exogenous supply channels (applied after delay)
# Note: GA is included here because testing perturbations use +GA even though
# training S6 only had 8 channels. Training runs unaffected (0 in all rows).
EXOGENOUS_NAMES = [
    "SL", "CK", "Aux_Shoot", "ABA", "SUC",
    "Low_R_FR", "Perception_SL", "High_R_FR", "GA",
]


def wt_state() -> Dict[str, float]:
    """Initial state: every node at WT level (1.0)."""
    return {n: 1.0 for n in NODE_NAMES}


def default_genotype() -> Dict[str, float]:
    """Default WT genotype: every modifier = 1.0."""
    return {g: 1.0 for g in GENE_NAMES}


def next_step(s: Dict[str, float], g: Dict[str, float],
              exo: Optional[Dict[str, float]] = None) -> Dict[str, float]:
    """
    Compute the next state given current state s and genotype g.
    Equations are an exact port of Script/nextStep.R lines 7-38.
    Exogenous supply (if provided) is added additively after computing
    the base equations.
    """
    n = {}

    # dat$SL[t] = ((2/(1 + dat$CXE[t-1])) * (gen$SL_n)) * (dat$LBO[t-1])
    n["SL"] = ((2 / (1 + s["CXE"])) * g["SL_n"]) * s["LBO"]

    # dat$BRC1_2[t] = (2 * ((PIFs + SPL9_15)/2) / (1 + (CK + SUC + PhyB)/3)) * BRC1_2_n
    n["BRC1_2"] = (2 * ((s["PIFs"] + s["SPL9_15"]) / 2)
                   / (1 + (s["CK"] + s["SUC"] + s["PhyB"]) / 3)) * g["BRC1_2_n"]

    # dat$CK[t] = (2 * (SUC) / (1 + (Aux_Shoot + Perception_SL)/2)) * CK_n
    n["CK"] = (2 * s["SUC"] / (1 + (s["Aux_Shoot"] + s["Perception_SL"]) / 2)) * g["CK_n"]

    # dat$Aux_Shoot[t] = (2 / (1 + Decap_signal)) * Aux_n
    n["Aux_Shoot"] = (2 / (1 + s["Decap_signal"])) * g["Aux_n"]

    # dat$ABA[t] = NCED3 * ABA_n
    n["ABA"] = s["NCED3"] * g["ABA_n"]

    # dat$Bud_release[t] = 2 * ((SUC + CK)/2) / (1 + (ABA + BRC1_2)/2)
    n["Bud_release"] = 2 * ((s["SUC"] + s["CK"]) / 2) / (1 + (s["ABA"] + s["BRC1_2"]) / 2)

    # dat$PAT_Bud[t] = (2 * ((Aux_Bud + CK + Decap_signal)/3) / (1 + Aux_Shoot)) * PIN_1_3_4_7_n
    n["PAT_Bud"] = (2 * ((s["Aux_Bud"] + s["CK"] + s["Decap_signal"]) / 3)
                    / (1 + s["Aux_Shoot"])) * g["PIN_1_3_4_7_n"]

    # dat$Aux_Bud[t] = 1 * Aux_Bud_n
    n["Aux_Bud"] = 1 * g["Aux_Bud_n"]

    # dat$SUC[t] = Decap_signal * SUC_n
    n["SUC"] = s["Decap_signal"] * g["SUC_n"]

    # dat$Sustained_growth[t] = ((SUC + min(GA, PAT_Bud))/2) * Bud_release
    n["Sustained_growth"] = ((s["SUC"] + min(s["GA"], s["PAT_Bud"])) / 2) * s["Bud_release"]

    # dat$GA[t] = Aux_Bud * GA_n
    n["GA"] = s["Aux_Bud"] * g["GA_n"]

    # dat$Decap_signal[t] = 1 * Decap_signal_n
    n["Decap_signal"] = 1 * g["Decap_signal_n"]

    # dat$CXE[t] = 1 * CXE_n
    n["CXE"] = 1 * g["CXE_n"]

    # dat$SMXL6_7_8[t] = (2 * ((SMXL6_7_8_Transcription + FHY3_FAR1)/2) / (1 + Perception_SL)) * SMXL6_7_8_n
    n["SMXL6_7_8"] = (2 * ((s["SMXL6_7_8_Transcription"] + s["FHY3_FAR1"]) / 2)
                      / (1 + s["Perception_SL"])) * g["SMXL6_7_8_n"]

    # dat$FHY3_FAR1[t] = (2 * High_R_FR / (1 + Low_R_FR)) * FHY3_FAR1_n
    n["FHY3_FAR1"] = (2 * s["High_R_FR"] / (1 + s["Low_R_FR"])) * g["FHY3_FAR1_n"]

    # dat$SPL9_15[t] = (2 / (1 + (FHY3_FAR1 + SMXL6_7_8)/2)) * SPL9_15_n
    n["SPL9_15"] = (2 / (1 + (s["FHY3_FAR1"] + s["SMXL6_7_8"]) / 2)) * g["SPL9_15_n"]

    # dat$SMXL6_7_8_Transcription[t] = 2 / (1 + SMXL6_7_8)
    n["SMXL6_7_8_Transcription"] = 2 / (1 + s["SMXL6_7_8"])

    # dat$Low_R_FR[t] = 1
    n["Low_R_FR"] = 1

    # dat$MAX2[t] = (2 / (1 + (SUC + PhyB)/2)) * MAX2_n
    n["MAX2"] = (2 / (1 + (s["SUC"] + s["PhyB"]) / 2)) * g["MAX2_n"]

    # dat$D14[t] = 1 * D14_n
    n["D14"] = 1 * g["D14_n"]

    # dat$Perception_SL[t] = 1 * min(D14, MAX2, SL)
    n["Perception_SL"] = 1 * min(s["D14"], s["MAX2"], s["SL"])

    # dat$D27[t] = 1 * D27_n
    n["D27"] = 1 * g["D27_n"]

    # dat$MAX3[t] = ((2 * Aux_Shoot / (1 + Perception_SL)) * MAX3_n) * D27
    n["MAX3"] = ((2 * s["Aux_Shoot"] / (1 + s["Perception_SL"])) * g["MAX3_n"]) * s["D27"]

    # dat$MAX4[t] = ((2 * Aux_Shoot / (1 + (Perception_SL + PhyB)/2)) * MAX4_n) * MAX3
    n["MAX4"] = ((2 * s["Aux_Shoot"] / (1 + (s["Perception_SL"] + s["PhyB"]) / 2))
                 * g["MAX4_n"]) * s["MAX3"]

    # dat$MAX1[t] = ((2 / (1 + Perception_SL)) * MAX1_n) * MAX4
    n["MAX1"] = ((2 / (1 + s["Perception_SL"])) * g["MAX1_n"]) * s["MAX4"]

    # dat$LBO[t] = (1 * LBO_n) * CLAMT
    n["LBO"] = (1 * g["LBO_n"]) * s["CLAMT"]

    # dat$CLAMT[t] = (1 * CLAMT_n) * MAX1
    n["CLAMT"] = (1 * g["CLAMT_n"]) * s["MAX1"]

    # dat$PhyB[t] = (2 * High_R_FR / (1 + Low_R_FR)) * PhyB_n
    n["PhyB"] = (2 * s["High_R_FR"] / (1 + s["Low_R_FR"])) * g["PhyB_n"]

    # dat$High_R_FR[t] = 1
    n["High_R_FR"] = 1

    # dat$PIFs[t] = (2 / (1 + High_R_FR)) * PIFs_n
    n["PIFs"] = (2 / (1 + s["High_R_FR"])) * g["PIFs_n"]

    # dat$HB21_40_53[t] = min(BRC1_2, Low_R_FR) * HB21_40_53_n
    n["HB21_40_53"] = min(s["BRC1_2"], s["Low_R_FR"]) * g["HB21_40_53_n"]

    # dat$NCED3[t] = HB21_40_53 * NCED3_n
    n["NCED3"] = s["HB21_40_53"] * g["NCED3_n"]

    # Apply exogenous supply additively
    if exo:
        for node_name, value in exo.items():
            if node_name in n and value != 0:
                n[node_name] += value

    return n


def simulate(genotype: Optional[Dict[str, float]] = None,
             exogenous: Optional[Dict[str, float]] = None,
             max_steps: int = 100,
             delay: int = 20,
             tol: float = 1e-8) -> Dict[str, float]:
    """
    Run the simulation until steady state.

    Args:
        genotype: dict of gene modifier values (missing keys default to 1.0)
        exogenous: dict of exogenous supply (applied after step `delay`)
        max_steps: maximum iterations
        delay: exogenous supply applied only after step `delay`
        tol: convergence tolerance (max abs diff between consecutive states)

    Returns:
        Final state dict (32 node values)
    """
    g = default_genotype()
    if genotype:
        for k, v in genotype.items():
            if k in g:
                g[k] = float(v)

    state = wt_state()

    for step in range(1, max_steps + 1):
        exo = exogenous if step > delay else None
        new_state = next_step(state, g, exo)

        if step > delay:
            max_diff = max(abs(new_state[n] - state[n]) for n in NODE_NAMES)
            if max_diff < tol:
                state = new_state
                break
        state = new_state

    return state


def compute_bin(perturbation_value: float, control_value: float,
                threshold: float = 0.05) -> int:
    """Discretise perturbation vs control into -1 / 0 / +1 bin (5% threshold)."""
    if control_value == 0:
        if perturbation_value > 0:
            return 1
        elif perturbation_value < 0:
            return -1
        else:
            return 0
    ratio = perturbation_value / control_value
    if ratio > 1 + threshold:
        return 1
    elif ratio < 1 - threshold:
        return -1
    else:
        return 0


if __name__ == "__main__":
    # Self-test: WT convergence
    wt = simulate()
    print("WT steady-state values:")
    for name in NODE_NAMES:
        print(f"  {name:30s}: {wt[name]:.6f}")
    # Sanity: all nodes should be ~1.0 for WT
    max_dev = max(abs(wt[n] - 1.0) for n in NODE_NAMES)
    print(f"\nMax deviation from 1.0: {max_dev:.2e}")
    assert max_dev < 1e-6, "WT should converge to 1.0 for all nodes!"
    print("WT convergence test: PASSED")

    # Test a known perturbation: MAX4_n (knockout) should increase Sustained_growth
    max4_ko = simulate(genotype={"MAX4_n": 0})
    print(f"\nMAX4_n knockout Sustained_growth: {max4_ko['Sustained_growth']:.4f}")
    print(f"WT Sustained_growth: {wt['Sustained_growth']:.4f}")
    bin_val = compute_bin(max4_ko["Sustained_growth"], wt["Sustained_growth"])
    print(f"Bin: {bin_val} (expected: +1 for increased branching)")
