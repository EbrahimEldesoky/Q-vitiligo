"""
Tests for the Q-Vitiligo Quantum-Enhanced Discovery Pipeline.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# Ensure the project root is on the path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qvitiligo.targets import TargetRegistry
from qvitiligo.chemoinformatics import (
    CompoundLibrary, compute_descriptors, assess_admet,
)
from qvitiligo.quantum_engine import (
    build_interaction_hamiltonian, QuantumBindingSolver, test_ibm_connection,
)
from qvitiligo.synergy_optimizer import (
    SynergyFitnessEvaluator, SynergyOptimizer, FormulationCandidate,
)


DATA_DIR = Path(__file__).resolve().parent.parent / "data"


# ─── Target Registry Tests ──────────────────────────────────────────────────

class TestTargetRegistry:
    def test_load_targets(self):
        registry = TargetRegistry(DATA_DIR / "targets.json")
        targets = registry.all_targets()
        assert len(targets) >= 5
        assert registry.get("JAK1") is not None
        assert registry.get("MC1R") is not None

    def test_filter_by_role(self):
        registry = TargetRegistry(DATA_DIR / "targets.json")
        immune = registry.get_by_role("autoimmune_arrest")
        assert len(immune) >= 2
        assert all(t.role == "autoimmune_arrest" for t in immune)

    def test_binding_weights_sum(self):
        registry = TargetRegistry(DATA_DIR / "targets.json")
        total = sum(t.binding_weight for t in registry.all_targets())
        assert abs(total - 1.0) < 0.01


# ─── Chemoinformatics Tests ─────────────────────────────────────────────────

class TestChemoinformatics:
    def test_compute_descriptors_piperine(self):
        smiles = "C1CCN(CC1)C(=O)/C=C/C=C/C2=CC3=C(C=C2)OCO3"
        desc, mol = compute_descriptors(smiles)
        assert desc is not None
        assert 200 < desc.molecular_weight < 400
        assert desc.num_heavy_atoms > 15

    def test_admet_scoring(self):
        smiles = "C1CCN(CC1)C(=O)/C=C/C=C/C2=CC3=C(C=C2)OCO3"
        desc, _ = compute_descriptors(smiles)
        admet = assess_admet(desc, "High")
        assert 0 <= admet.oral_absorption <= 1
        assert 0 <= admet.bioavailability_score <= 1
        assert admet.lipinski_violations <= 4

    def test_compound_library_load(self):
        library = CompoundLibrary(DATA_DIR / "natural_compounds.json")
        compounds = library.all_compounds()
        assert len(compounds) >= 8
        for c in compounds:
            assert c.descriptors is not None
            assert c.admet is not None

    def test_filter_by_target(self):
        library = CompoundLibrary(DATA_DIR / "natural_compounds.json")
        jak1_compounds = library.filter_by_target("JAK1")
        assert len(jak1_compounds) >= 2


# ─── Quantum Engine Tests ───────────────────────────────────────────────────

class TestQuantumEngine:
    def test_build_hamiltonian(self):
        smiles = "C1CCN(CC1)C(=O)/C=C/C=C/C2=CC3=C(C=C2)OCO3"
        system = build_interaction_hamiltonian(smiles, "JAK1", num_qubits=4)
        assert system.num_qubits == 4
        assert system.hamiltonian is not None
        assert len(system.hamiltonian) > 0

    def test_vqe_solver_aer(self):
        solver = QuantumBindingSolver(backend="aer", num_qubits=4)
        result = solver.solve_binding(
            compound_name="Piperine",
            compound_smiles="C1CCN(CC1)C(=O)/C=C/C=C/C2=CC3=C(C=C2)OCO3",
            target_id="MC1R",
            target_name="Melanocortin 1 Receptor",
            max_iterations=50,
        )
        assert result.ground_state_energy != 0.0
        assert result.binding_energy_kcal != 0.0
        assert result.num_qubits == 4
        assert result.iterations > 0

    def test_different_compounds_different_energies(self):
        solver = QuantumBindingSolver(backend="aer", num_qubits=4)
        r1 = solver.solve_binding("Piperine",
            "C1CCN(CC1)C(=O)/C=C/C=C/C2=CC3=C(C=C2)OCO3",
            "JAK1", "JAK1", max_iterations=50)
        r2 = solver.solve_binding("Thymoquinone",
            "CC1=CC(=O)C(=CC1=O)C(C)C",
            "JAK1", "JAK1", max_iterations=50)
        # Different molecules should produce different energies
        assert r1.ground_state_energy != r2.ground_state_energy


# ─── Synergy Optimizer Tests ────────────────────────────────────────────────

class TestSynergyOptimizer:
    def _setup(self):
        library = CompoundLibrary(DATA_DIR / "natural_compounds.json")
        compounds = {c.id: c for c in library.all_compounds()}
        # Mock binding energies
        binding_energies = {}
        for cid, comp in compounds.items():
            binding_energies[cid] = {}
            for tid in comp.targets:
                binding_energies[cid][tid] = -np.random.uniform(10, 40)
        return compounds, binding_energies

    def test_fitness_evaluation(self):
        compounds, binding_energies = self._setup()
        evaluator = SynergyFitnessEvaluator(compounds, binding_energies)
        candidate = FormulationCandidate(
            compound_ids=list(compounds.keys()),
            concentrations_uM=np.array([20.0] * len(compounds)),
        )
        fitness = evaluator.evaluate(candidate)
        assert fitness > 0
        assert "immune_arrest" in candidate.fitness_breakdown

    def test_optimizer_convergence(self):
        compounds, binding_energies = self._setup()
        optimizer = SynergyOptimizer(
            compounds, binding_energies,
            population_size=20,
        )
        result = optimizer.optimize(max_generations=15, verbose=False)
        assert result.best_formulation.fitness > 0
        assert result.generations_run > 0
        assert len(result.top_formulations) <= 5


# ─── IBM Connection Test (optional) ─────────────────────────────────────────

class TestIBMConnection:
    @pytest.mark.skipif(True, reason="Requires IBM API key; run manually")
    def test_ibm_connectivity(self):
        assert test_ibm_connection()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
