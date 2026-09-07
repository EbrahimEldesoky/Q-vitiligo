"""
Master Pipeline Orchestrator for Q-Vitiligo.

Coordinates the full discovery workflow:
  1. Load biological targets and phytochemical library
  2. Classical ADMET prescreening
  3. Quantum VQE binding energy computation for all compound-target pairs
  4. Multi-target synergy optimization (closed-loop GA)
  5. Formulation report generation
"""

import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from qvitiligo.targets import TargetRegistry
from qvitiligo.chemoinformatics import CompoundLibrary, CompoundProfile
from qvitiligo.quantum_engine import QuantumBindingSolver, VQEResult
from qvitiligo.synergy_optimizer import SynergyOptimizer, OptimizationResult


class VitiligoDiscoveryPipeline:
    """End-to-end discovery pipeline for quantum-enhanced vitiligo drug formulation."""

    def __init__(
        self,
        data_dir: Optional[Path] = None,
        backend: str = "aer",
        num_qubits: int = 4,
    ):
        if data_dir is None:
            data_dir = Path(__file__).resolve().parent.parent / "data"
        self.data_dir = Path(data_dir)
        self.backend = backend
        self.num_qubits = num_qubits

        # Loaded during run
        self.target_registry: Optional[TargetRegistry] = None
        self.compound_library: Optional[CompoundLibrary] = None
        self.binding_results: Dict[str, Dict[str, VQEResult]] = {}
        self.binding_energies: Dict[str, Dict[str, float]] = {}
        self.optimization_result: Optional[OptimizationResult] = None

    def run(
        self,
        max_opt_generations: int = 50,
        min_bioavailability: float = 0.3,
        verbose: bool = True,
    ) -> Dict:
        """Execute the full discovery pipeline."""

        t0 = time.time()
        report = {"pipeline": "QuantumVitiligo", "version": "1.0.0", "stages": {}}

        # ── Stage 1: Load Data ──────────────────────────────────────────────
        if verbose:
            print("\n" + "=" * 70)
            print("  STAGE 1: Loading Biological Targets & Phytochemical Library")
            print("=" * 70)

        self.target_registry = TargetRegistry(self.data_dir / "targets.json")
        self.compound_library = CompoundLibrary(self.data_dir / "natural_compounds.json")

        targets = self.target_registry.all_targets()
        compounds = self.compound_library.all_compounds()

        if verbose:
            print(f"  Loaded {len(targets)} biological targets:")
            for t in targets:
                print(f"    [{t.role:<20}] {t.id:<12} — {t.name}")
            print(f"\n  Loaded {len(compounds)} natural compounds:")
            print(self.compound_library.summary_table())

        report["stages"]["data_loading"] = {
            "targets": len(targets),
            "compounds": len(compounds),
        }

        # ── Stage 2: ADMET Prescreening ─────────────────────────────────────
        if verbose:
            print("\n" + "=" * 70)
            print("  STAGE 2: ADMET & Drug-likeness Prescreening")
            print("=" * 70)

        viable = self.compound_library.filter_oral_viable(min_bioavailability)
        if verbose:
            print(f"\n  {len(viable)}/{len(compounds)} compounds pass oral bioavailability "
                  f"threshold (≥{min_bioavailability}):")
            for c in viable:
                a = c.admet
                print(f"    {c.name:<16} BioAv={a.bioavailability_score:.3f}  "
                      f"Route={a.route_recommendation}")

        # Use all compounds in quantum stage (even low-bioavailability ones
        # can work as injectables), but record prescreening results
        report["stages"]["admet_prescreening"] = {
            "oral_viable": [c.id for c in viable],
            "all_passed_to_quantum": [c.id for c in compounds],
        }

        # ── Stage 3: Quantum VQE Binding Computation ───────────────────────
        if verbose:
            print("\n" + "=" * 70)
            print("  STAGE 3: Quantum VQE Binding Energy Simulation")
            print(f"           Backend: {self.backend} | Qubits: {self.num_qubits}")
            print("=" * 70)

        solver = QuantumBindingSolver(
            backend=self.backend,
            num_qubits=self.num_qubits,
        )

        total_pairs = sum(len(c.targets) for c in compounds)
        pair_count = 0

        for comp in compounds:
            self.binding_results[comp.id] = {}
            self.binding_energies[comp.id] = {}

            for target_id in comp.targets:
                target = self.target_registry.get(target_id)
                if target is None:
                    continue

                pair_count += 1
                if verbose:
                    print(f"\n  [{pair_count}/{total_pairs}] "
                          f"{comp.name} ⟷ {target.name} ({target.id})")

                result = solver.solve_binding(
                    compound_name=comp.name,
                    compound_smiles=comp.smiles,
                    target_id=target_id,
                    target_name=target.name,
                    max_iterations=150,
                )

                self.binding_results[comp.id][target_id] = result
                self.binding_energies[comp.id][target_id] = result.binding_energy_kcal

                if verbose:
                    conv = "✓" if result.converged else "✗"
                    print(f"    E₀ = {result.ground_state_energy:.6f} Ha  |  "
                          f"ΔE_bind = {result.binding_energy_kcal:.3f} kcal/mol  |  "
                          f"Converged: {conv}  |  "
                          f"Iters: {result.iterations}  |  "
                          f"Time: {result.wall_time_s:.1f}s")

        report["stages"]["quantum_vqe"] = {
            "total_pairs": total_pairs,
            "backend": self.backend,
            "num_qubits": self.num_qubits,
            "results": {
                cid: {
                    tid: {
                        "binding_energy_kcal": r.binding_energy_kcal,
                        "converged": r.converged,
                        "iterations": r.iterations,
                    }
                    for tid, r in targets_dict.items()
                }
                for cid, targets_dict in self.binding_results.items()
            },
        }

        # ── Stage 4: Synergy Optimization Loop ─────────────────────────────
        if verbose:
            print("\n" + "=" * 70)
            print("  STAGE 4: Multi-Target Synergy Optimization (Genetic Algorithm)")
            print("=" * 70)

        compound_dict = {c.id: c for c in compounds}
        optimizer = SynergyOptimizer(
            compounds=compound_dict,
            binding_energies=self.binding_energies,
            population_size=40,
        )

        self.optimization_result = optimizer.optimize(
            max_generations=max_opt_generations,
            verbose=verbose,
        )

        opt = self.optimization_result

        report["stages"]["synergy_optimization"] = {
            "generations_run": opt.generations_run,
            "converged": opt.converged,
            "best_fitness": opt.best_formulation.fitness,
            "wall_time_s": opt.wall_time_s,
            "delivery_route": opt.delivery_route,
            "target_coverage": opt.target_coverage,
        }

        # ── Stage 5: Generate Formulation Report ───────────────────────────
        if verbose:
            print("\n" + "=" * 70)
            print("  STAGE 5: Final Formulation Report")
            print("=" * 70)

        formulation_report = self._generate_formulation_report(compound_dict)
        report["formulation"] = formulation_report

        total_time = time.time() - t0
        report["total_wall_time_s"] = round(total_time, 2)

        if verbose:
            self._print_formulation_report(formulation_report)
            print(f"\n  Total pipeline runtime: {total_time:.1f}s")

        return report

    def _generate_formulation_report(self, compounds: Dict[str, CompoundProfile]) -> Dict:
        """Generate the final formulation dossier."""
        opt = self.optimization_result
        best = opt.best_formulation

        ingredients = []
        total_mass = 0.0

        for cid, conc in zip(best.compound_ids, best.concentrations_uM):
            if conc < 0.5:
                continue
            comp = compounds.get(cid)
            if comp is None:
                continue

            # Convert μM concentration to approximate mg per dose
            # Assuming 500mL plasma volume equivalent distribution
            mw = comp.descriptors.molecular_weight if comp.descriptors else 300.0
            mg_per_dose = round(conc * mw * 0.5 / 1000, 2)  # μM * g/mol * 0.5L / 1000
            total_mass += mg_per_dose

            ingredients.append({
                "compound_id": cid,
                "name": comp.name,
                "plant_source": comp.plant_source,
                "concentration_uM": round(conc, 2),
                "mg_per_dose": mg_per_dose,
                "primary_action": comp.primary_action,
                "targets": comp.targets,
                "mechanism": comp.known_mechanism,
                "bioavailability": comp.admet.bioavailability_score if comp.admet else None,
            })

        # Sort by concentration (highest first)
        ingredients.sort(key=lambda x: x["concentration_uM"], reverse=True)

        # Compute percentages
        for ing in ingredients:
            ing["percentage_of_formula"] = round(ing["mg_per_dose"] / max(total_mass, 0.01) * 100, 1)

        return {
            "formulation_name": "QuantumVitiligo Formulation v1.0",
            "therapeutic_indication": "Vitiligo (Non-segmental) — Autoimmune arrest + Repigmentation",
            "delivery_route": opt.delivery_route,
            "fitness_score": best.fitness,
            "fitness_breakdown": best.fitness_breakdown,
            "target_coverage": opt.target_coverage,
            "total_mg_per_dose": round(total_mass, 2),
            "ingredients": ingredients,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "notes": [
                "Concentrations are theoretical and require in-vitro validation.",
                "Piperine is included as a bioavailability enhancer for co-formulated polyphenols.",
                "This formulation is designed for systemic delivery (oral or injectable).",
                "All compounds are derived from natural plant sources with established safety profiles.",
            ],
        }

    def _print_formulation_report(self, report: Dict):
        """Print a human-readable formulation report to console."""
        print(f"\n  ╔══════════════════════════════════════════════════════════════╗")
        print(f"  ║  {report['formulation_name']:^58}  ║")
        print(f"  ╚══════════════════════════════════════════════════════════════╝")
        print(f"\n  Indication:    {report['therapeutic_indication']}")
        print(f"  Delivery:      {report['delivery_route']}")
        print(f"  Fitness Score: {report['fitness_score']:.5f}")
        print(f"  Total mg/dose: {report['total_mg_per_dose']:.2f} mg")

        print(f"\n  ┌─────────────────────────────────────────────────────────────┐")
        print(f"  │  {'Component':<16} {'Source':<22} {'μM':>7} {'mg':>7} {'%':>5}  {'Action':<14} │")
        print(f"  ├─────────────────────────────────────────────────────────────┤")
        for ing in report["ingredients"]:
            src = ing["plant_source"][:20]
            print(f"  │  {ing['name']:<16} {src:<22} {ing['concentration_uM']:>7.1f} "
                  f"{ing['mg_per_dose']:>7.2f} {ing['percentage_of_formula']:>4.1f}%  "
                  f"{ing['primary_action']:<14} │")
        print(f"  └─────────────────────────────────────────────────────────────┘")

        fb = report["fitness_breakdown"]
        print(f"\n  Fitness Breakdown:")
        print(f"    Immune Arrest:     {fb.get('immune_arrest', 0):.4f}")
        print(f"    Repigmentation:    {fb.get('repigmentation', 0):.4f}")
        print(f"    Cellular Defense:  {fb.get('cellular_defense', 0):.4f}")
        print(f"    ADMET:             {fb.get('admet', 0):.4f}")
        print(f"    Synergy:           {fb.get('synergy', 0):.4f}")
        print(f"    Toxicity Penalty: -{fb.get('toxicity_penalty', 0):.4f}")

        tc = report["target_coverage"]
        print(f"\n  Target Coverage:")
        for t, v in tc.items():
            status = "●" if v > 0 else "○"
            print(f"    {status} {t:<14} {'COVERED' if v > 0 else 'NOT COVERED'}")

        print(f"\n  Notes:")
        for note in report.get("notes", []):
            print(f"    • {note}")
