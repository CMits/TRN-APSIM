"""A FlashP-topology simulator that uses PSoup's update form.

Why: FlashP's algebraic engine combines inhibitors as a clamped product

    inhibition = min(1 / max(prod(inhibitors), epsilon), K)      epsilon=0.1, K=10

which saturates hard. Once a mutation drives the inhibitor product below 0.1 the
term pins at 10 and stops responding to anything — including sucrose. Combined
with activator_floor=0.01 at the phenotype node, that is what collapsed 22 of 37
genotypes onto an identical sqrt(SUC) curve and flattened 6 more.

PSoup instead uses, for every node:

    node = mean(activators) * 2/(1 + mean(inhibitors)) * gene_modifier

which is smooth and unbounded — no floor, no ceiling. Both forms are WT-normalised
(all inputs 1 -> node 1), which is the invariant FlashP's own validator checks.

This module applies PSoup's form to a FlashP network's topology, so the two can be
compared with topology held constant. Nothing in FlashP_Repo is modified.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple


class PSoupFormSimulator:
    """Fixed-point propagation over a FlashP topology using PSoup's node rule."""

    def __init__(self, equations_path: str | Path,
                 damping: float = 0.7,
                 tol: float = 1e-9,
                 max_iter: int = 500):
        spec = json.loads(Path(equations_path).read_text(encoding="utf-8"))
        self.nodes: List[str] = []
        self.act: Dict[str, List[str]] = {}
        self.inh: Dict[str, List[str]] = {}
        self.src: Dict[str, bool] = {}
        for e in spec["equations"]:
            n = e["n"]
            self.nodes.append(n)
            self.act[n] = list(e.get("a") or [])
            self.inh[n] = list(e.get("inh") or [])
            self.src[n] = bool(e.get("src"))
        self.damping = damping
        self.tol = tol
        self.max_iter = max_iter

    @staticmethod
    def _activation(vals: List[float]) -> float:
        """PSoup: arithmetic mean of activators; 1.0 when there are none."""
        return sum(vals) / len(vals) if vals else 1.0

    @staticmethod
    def _inhibition(vals: List[float]) -> float:
        """PSoup: 2/(1 + mean(inhibitors)); 1.0 when there are none. No clamp."""
        return 2.0 / (1.0 + sum(vals) / len(vals)) if vals else 1.0

    def simulate(self,
                 gene_modifiers: Dict[str, float],
                 exogenous_supply: Dict[str, float] | None = None
                 ) -> Tuple[Dict[str, float], bool, int]:
        exo = exogenous_supply or {}
        state = {n: 1.0 for n in self.nodes}

        converged = False
        it = 0
        for it in range(1, self.max_iter + 1):
            new = {}
            for n in self.nodes:
                g = float(gene_modifiers.get(n, 1.0))
                if self.src[n]:
                    # a source has no upstream regulator: it is just its own level
                    v = 1.0 * g
                else:
                    a = self._activation([state[x] for x in self.act[n]])
                    i = self._inhibition([state[x] for x in self.inh[n]])
                    v = a * i * g
                v += float(exo.get(n, 0.0))
                new[n] = max(v, 0.0)

            delta = max(abs(new[n] - state[n]) for n in self.nodes)
            # damped update, matching FlashP's iteration style
            state = {n: (1 - self.damping) * new[n] + self.damping * state[n]
                     for n in self.nodes}
            if delta < self.tol:
                converged = True
                break

        return state, converged, it
