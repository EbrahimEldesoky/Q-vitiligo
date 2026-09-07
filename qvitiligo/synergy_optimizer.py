"""
Multi-Target Synergy Optimizer for Q-Vitiligo.

Implements a closed-loop evolutionary optimization algorithm that searches for
the optimal extract formulation — a concentration vector of natural compounds —
maximizing a multi-objective fitness function encompassing:

  1. Autoimmune arrest efficacy (JAK1/JAK2 inhibition)
  2. Repigmentation stimulation (MC1R/TYR activation)
  3. Cellular defense (Nrf2/KEAP1 activation)
  4. ADMET / oral bioavailability
  5. Synergy bonuses & toxicity penalties

The optimizer uses a genetic algorithm (GA) with tournament selection,
crossover, and adaptive mutation to evolve a population of candidate
formulations toward the Pareto front.
"""

import math
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from qvitiligo.chemoinformatics import CompoundProfile, ADMETScore
from qvitiligo.quantum_engine import VQEResult


# ---------------------------------------------------------------------------
# Data Structures
# ---------------------------------------------------------------------------

@dataclass
class FormulationCandidate:
    """A candidate extract formulation: a vector of compound concentrations."""
    compound_ids: List[str]
    concentrations_uM: np.ndarray   # concentration of each compound in μM
    fitness: float = 0.0
    fitness_breakdown: Dict[str, float] = field(default_factory=dict)
    generation: int = 0
    rank: int = 0

    def __repr__(self):
        parts = [f"{cid}:{c:.1f}μM" for cid, c in zip(self.compound_ids, self.concentrations_uM)]
        return f"Formulation(fitness={self.fitness:.4f}, {', '.join(parts)})"


@dataclass
class OptimizationResult:
    """Result of the synergy optimization loop."""
    best_formulation: FormulationCandidate
    top_formulations: List[FormulationCandidate]
    generations_run: int
    converged: bool
    convergence_history: List[float]
    wall_time_s: float
    target_coverage: Dict[str, float]
    delivery_route: str


# ---------------------------------------------------------------------------
# Fitness Function
# ---------------------------------------------------------------------------

class SynergyFitnessEvaluator:
    """
    Multi-objective fitness evaluator for vitiligo extract formulations.
    """

    # Weights for each therapeutic arm
    WEIGHT_IMMUNE = 0.35       # Autoimmune arrest (JAK1/JAK2)
    WEIGHT_REPIG = 0.30        # Repigmentation (MC1R/TYR)
    WEIGHT_DEFENSE = 0.15      # Cellular defense (Nrf2)
    WEIGHT_ADMET = 0.10        # Drug-likeness / bioavailability
    WEIGHT_SYNERGY = 0.10      # Synergy bonus

    def __init__(
        self,
        compounds: Dict[str, CompoundProfile],
        binding_energies: Dict[str, Dict[str, float]],
        target_weights: Optional[Dict[str, float]] = None,
    ):
        """
        Args:
            compounds: dict of compound_id -> CompoundProfile
            binding_energies: dict of compound_id -> {target_id -> binding_energy_kcal}
            target_weights: optional override for target importance weights
        """
        self.compounds = compounds
        self.binding_energies = binding_energies
        self.target_weights = target_weights or {
            "JAK1": 0.35, "JAK2": 0.20, "MC1R": 0.25,
            "TYR": 0.10, "Nrf2_KEAP1": 0.10,
        }

    def evaluate(self, candidate: FormulationCandidate) -> float:
        """Compute the multi-objective fitness score for a formulation candidate."""

        immune_score = 0.0
        repig_score = 0.0
        defense_score = 0.0
        admet_score = 0.0
        synergy_score = 0.0

        active_compounds = []

        for i, cid in enumerate(candidate.compound_ids):
            comp = self.compounds.get(cid)
            if comp is None:
                continue
            conc = candidate.concentrations_uM[i]
            if conc <= 0:
                continue

            active_compounds.append((cid, comp, conc))

            # Concentration efficacy curve (sigmoid / Hill equation)
            # EC50 is midpoint of compound's concentration range
            ec50 = (comp.concentration_range_uM[0] + comp.concentration_range_uM[1]) / 2
            hill_coeff = 1.5
            efficacy = (conc ** hill_coeff) / (ec50 ** hill_coeff + conc ** hill_coeff)

            # Get binding energies for this compound
            comp_bindings = self.binding_energies.get(cid, {})

            # Score per target
            for target_id in comp.targets:
                be = comp_bindings.get(target_id, 0.0)
                # More negative binding energy = stronger binding
                binding_strength = max(0.0, min(1.0, abs(be) / 50.0))
                target_score = efficacy * binding_strength * self.target_weights.get(target_id, 0.1)

                if target_id in ("JAK1", "JAK2"):
                    immune_score += target_score
                elif target_id in ("MC1R", "TYR"):
                    repig_score += target_score
                elif target_id == "Nrf2_KEAP1":
                    defense_score += target_score

            # ADMET contribution
            if comp.admet:
                admet_score += comp.admet.bioavailability_score * efficacy * 0.2

        # Synergy bonus: reward diverse multi-target coverage
        targets_covered = set()
        for cid, comp, conc in active_compounds:
            if conc > 0:
                targets_covered.update(comp.targets)
        coverage_ratio = len(targets_covered) / max(len(self.target_weights), 1)
        synergy_score = coverage_ratio

        # Piperine bioenhancer bonus (increases absorption of co-formulated compounds)
        piperine_conc = 0.0
        for cid, comp, conc in active_compounds:
            if comp.name == "Piperine" and conc > 0:
                piperine_conc = conc
        if piperine_conc > 5.0:
            bioenhance_factor = min(1.3, 1.0 + piperine_conc / 100.0)
            admet_score *= bioenhance_factor

        # Toxicity penalty: penalize if any compound exceeds safe range
        tox_penalty = 0.0
        for i, cid in enumerate(candidate.compound_ids):
            comp = self.compounds.get(cid)
            if comp is None:
                continue
            conc = candidate.concentrations_uM[i]
            upper = comp.concentration_range_uM[1]
            if conc > upper * 1.5:
                tox_penalty += 0.2 * ((conc - upper * 1.5) / upper)

        # Weighted total
        fitness = (
            self.WEIGHT_IMMUNE * min(immune_score, 1.0)
            + self.WEIGHT_REPIG * min(repig_score, 1.0)
            + self.WEIGHT_DEFENSE * min(defense_score, 1.0)
            + self.WEIGHT_ADMET * min(admet_score, 1.0)
            + self.WEIGHT_SYNERGY * synergy_score
            - tox_penalty
        )
        fitness = max(0.0, fitness)

        candidate.fitness = round(fitness, 6)
        candidate.fitness_breakdown = {
            "immune_arrest": round(immune_score, 4),
            "repigmentation": round(repig_score, 4),
            "cellular_defense": round(defense_score, 4),
            "admet": round(admet_score, 4),
            "synergy": round(synergy_score, 4),
            "toxicity_penalty": round(tox_penalty, 4),
        }

        return fitness


# ---------------------------------------------------------------------------
# Genetic Algorithm Optimizer
# ---------------------------------------------------------------------------

class SynergyOptimizer:
    """
    Genetic Algorithm (GA) for closed-loop formulation optimization.
    Evolves a population of concentration vectors to maximize multi-target
    synergistic fitness against vitiligo.
    """

    def __init__(
        self,
        compounds: Dict[str, CompoundProfile],
        binding_energies: Dict[str, Dict[str, float]],
        population_size: int = 40,
        elite_fraction: float = 0.15,
        mutation_rate: float = 0.25,
        crossover_rate: float = 0.7,
    ):
        self.compounds = compounds
        self.compound_ids = list(compounds.keys())
        self.n_compounds = len(self.compound_ids)
        self.binding_energies = binding_energies
        self.population_size = population_size
        self.elite_count = max(2, int(population_size * elite_fraction))
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate

        self.evaluator = SynergyFitnessEvaluator(compounds, binding_energies)

    def optimize(
        self,
        max_generations: int = 50,
        convergence_threshold: float = 1e-4,
        convergence_patience: int = 10,
        verbose: bool = True,
    ) -> OptimizationResult:
        """Run the GA optimization loop."""

        t0 = time.time()

        # Initialize population
        population = self._init_population()
        convergence_history = []
        best_fitness_ever = -float("inf")
        stagnation = 0

        for gen in range(max_generations):
            # Evaluate fitness
            for ind in population:
                self.evaluator.evaluate(ind)
                ind.generation = gen

            # Sort by fitness (descending)
            population.sort(key=lambda x: x.fitness, reverse=True)
            best = population[0]

            convergence_history.append(best.fitness)

            if verbose:
                avg_fit = np.mean([ind.fitness for ind in population])
                print(
                    f"  [Gen {gen+1:>3}] Best={best.fitness:.5f}  "
                    f"Avg={avg_fit:.5f}  "
                    f"Immune={best.fitness_breakdown.get('immune_arrest',0):.3f}  "
                    f"Repig={best.fitness_breakdown.get('repigmentation',0):.3f}  "
                    f"Defense={best.fitness_breakdown.get('cellular_defense',0):.3f}"
                )

            # Convergence check
            if best.fitness > best_fitness_ever + convergence_threshold:
                best_fitness_ever = best.fitness
                stagnation = 0
            else:
                stagnation += 1

            if stagnation >= convergence_patience:
                if verbose:
                    print(f"  [CONVERGED] No improvement for {convergence_patience} generations")
                break

            # Selection + Reproduction
            new_population = []

            # Elitism: keep top individuals
            for i in range(self.elite_count):
                elite = FormulationCandidate(
                    compound_ids=list(population[i].compound_ids),
                    concentrations_uM=population[i].concentrations_uM.copy(),
                )
                new_population.append(elite)

            # Fill the rest with crossover + mutation
            while len(new_population) < self.population_size:
                p1 = self._tournament_select(population)
                p2 = self._tournament_select(population)

                if np.random.random() < self.crossover_rate:
                    child = self._crossover(p1, p2)
                else:
                    child = FormulationCandidate(
                        compound_ids=list(p1.compound_ids),
                        concentrations_uM=p1.concentrations_uM.copy(),
                    )

                self._mutate(child)
                self._clamp_concentrations(child)
                new_population.append(child)

            population = new_population

        # Final evaluation
        for ind in population:
            self.evaluator.evaluate(ind)
        population.sort(key=lambda x: x.fitness, reverse=True)

        best = population[0]
        top_5 = population[:5]
        for i, ind in enumerate(top_5):
            ind.rank = i + 1

        # Target coverage analysis
        targets_hit = set()
        for cid, conc in zip(best.compound_ids, best.concentrations_uM):
            if conc > 0.5:
                comp = self.compounds.get(cid)
                if comp:
                    targets_hit.update(comp.targets)
        coverage = {
            t: (1.0 if t in targets_hit else 0.0)
            for t in ["JAK1", "JAK2", "MC1R", "TYR", "Nrf2_KEAP1"]
        }

        # Determine delivery route based on aggregate bioavailability
        avg_bio = np.mean([
            self.compounds[cid].admet.bioavailability_score
            for cid, conc in zip(best.compound_ids, best.concentrations_uM)
            if conc > 0.5 and self.compounds.get(cid) and self.compounds[cid].admet
        ]) if any(c > 0.5 for c in best.concentrations_uM) else 0.5

        route = "Oral (capsule/tablet)" if avg_bio >= 0.5 else "Injectable (subcutaneous)"

        wall_time = time.time() - t0

        return OptimizationResult(
            best_formulation=best,
            top_formulations=top_5,
            generations_run=len(convergence_history),
            converged=(stagnation >= convergence_patience),
            convergence_history=convergence_history,
            wall_time_s=round(wall_time, 2),
            target_coverage=coverage,
            delivery_route=route,
        )

    def _init_population(self) -> List[FormulationCandidate]:
        """Create initial population with random concentrations within safe ranges."""
        population = []
        for _ in range(self.population_size):
            concs = np.zeros(self.n_compounds)
            for i, cid in enumerate(self.compound_ids):
                comp = self.compounds.get(cid)
                if comp:
                    lo, hi = comp.concentration_range_uM
                    concs[i] = np.random.uniform(lo * 0.5, hi)
            population.append(FormulationCandidate(
                compound_ids=list(self.compound_ids),
                concentrations_uM=concs,
            ))
        return population

    def _tournament_select(self, population: List[FormulationCandidate], k: int = 3) -> FormulationCandidate:
        """Tournament selection."""
        indices = np.random.choice(len(population), size=min(k, len(population)), replace=False)
        best = max((population[i] for i in indices), key=lambda x: x.fitness)
        return best

    def _crossover(self, p1: FormulationCandidate, p2: FormulationCandidate) -> FormulationCandidate:
        """Blend crossover (BLX-α) for continuous optimization."""
        alpha = 0.3
        child_concs = np.zeros(self.n_compounds)
        for i in range(self.n_compounds):
            lo = min(p1.concentrations_uM[i], p2.concentrations_uM[i])
            hi = max(p1.concentrations_uM[i], p2.concentrations_uM[i])
            span = hi - lo
            child_concs[i] = np.random.uniform(lo - alpha * span, hi + alpha * span)
        return FormulationCandidate(
            compound_ids=list(self.compound_ids),
            concentrations_uM=child_concs,
        )

    def _mutate(self, ind: FormulationCandidate):
        """Gaussian mutation with adaptive step size."""
        for i in range(self.n_compounds):
            if np.random.random() < self.mutation_rate:
                comp = self.compounds.get(ind.compound_ids[i])
                if comp:
                    lo, hi = comp.concentration_range_uM
                    sigma = (hi - lo) * 0.15
                    ind.concentrations_uM[i] += np.random.normal(0, sigma)

    def _clamp_concentrations(self, ind: FormulationCandidate):
        """Clamp concentrations to valid ranges (allow zero for dropping compounds)."""
        for i in range(self.n_compounds):
            comp = self.compounds.get(ind.compound_ids[i])
            if comp:
                lo, hi = comp.concentration_range_uM
                ind.concentrations_uM[i] = max(0.0, min(hi * 1.5, ind.concentrations_uM[i]))
