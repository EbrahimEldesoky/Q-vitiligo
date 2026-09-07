"""
Chemoinformatics module for Q-Vitiligo.
Calculates molecular descriptors, drug-likeness (Lipinski / Veber),
ADMET scoring, oral bioavailability, and compound library management.
"""

import json
import math
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple

import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class MolecularDescriptors:
    """Core physicochemical descriptors for drug-likeness assessment."""
    molecular_weight: float = 0.0
    logp: float = 0.0
    hbd: int = 0           # hydrogen bond donors
    hba: int = 0           # hydrogen bond acceptors
    tpsa: float = 0.0      # topological polar surface area
    rotatable_bonds: int = 0
    aromatic_rings: int = 0
    num_heavy_atoms: int = 0
    num_atoms: int = 0

@dataclass
class ADMETScore:
    """Simplified ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) profile."""
    oral_absorption: float = 0.0     # 0-1 predicted fraction absorbed
    lipinski_violations: int = 0
    veber_pass: bool = True
    bioavailability_score: float = 0.0   # 0-1 composite
    clearance_penalty: float = 0.0       # 0-1 (0=low clearance, ideal)
    toxicity_flag: bool = False
    route_recommendation: str = "oral"   # "oral", "injectable", "topical"

@dataclass
class CompoundProfile:
    """Full compound profile: identity + descriptors + ADMET + target info."""
    id: str
    name: str
    plant_source: str
    smiles: str
    primary_action: str
    known_mechanism: str
    targets: List[str]
    oral_absorption_tier: str
    safety_profile: str
    concentration_range_uM: Tuple[float, float] = (1.0, 100.0)
    descriptors: Optional[MolecularDescriptors] = None
    admet: Optional[ADMETScore] = None
    mol: object = field(default=None, repr=False)  # RDKit Mol object


# ---------------------------------------------------------------------------
# Descriptor Calculation
# ---------------------------------------------------------------------------

def compute_descriptors(smiles: str) -> Tuple[Optional[MolecularDescriptors], Optional[object]]:
    """Compute molecular descriptors from a SMILES string using RDKit."""
    if not HAS_RDKIT:
        return _fallback_descriptors(smiles), None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  [WARN] RDKit could not parse SMILES: {smiles}")
        return _fallback_descriptors(smiles), None

    desc = MolecularDescriptors(
        molecular_weight=Descriptors.ExactMolWt(mol),
        logp=Crippen.MolLogP(mol),
        hbd=Lipinski.NumHDonors(mol),
        hba=Lipinski.NumHAcceptors(mol),
        tpsa=rdMolDescriptors.CalcTPSA(mol),
        rotatable_bonds=Lipinski.NumRotatableBonds(mol),
        aromatic_rings=Descriptors.NumAromaticRings(mol),
        num_heavy_atoms=mol.GetNumHeavyAtoms(),
        num_atoms=mol.GetNumAtoms(),
    )
    return desc, mol


def _fallback_descriptors(smiles: str) -> MolecularDescriptors:
    """Rough heuristic descriptors when RDKit is unavailable."""
    heavy = sum(1 for c in smiles if c.isalpha() and c.isupper())
    return MolecularDescriptors(
        molecular_weight=heavy * 14.0,
        logp=2.0,
        hbd=smiles.count("O") + smiles.count("N"),
        hba=smiles.count("O") + smiles.count("N"),
        tpsa=40.0,
        rotatable_bonds=smiles.count("-") // 2,
        aromatic_rings=smiles.lower().count("c1"),
        num_heavy_atoms=heavy,
        num_atoms=heavy,
    )


# ---------------------------------------------------------------------------
# ADMET & Drug-likeness Assessment
# ---------------------------------------------------------------------------

def assess_admet(desc: MolecularDescriptors, oral_tier: str = "Moderate") -> ADMETScore:
    """
    Score a compound's drug-likeness for systemic vitiligo therapy.
    Uses Lipinski Rule of Five and Veber rules as primary filters,
    then computes a composite bioavailability score.
    """
    # Lipinski Rule of Five violations
    violations = 0
    if desc.molecular_weight > 500:
        violations += 1
    if desc.logp > 5:
        violations += 1
    if desc.hbd > 5:
        violations += 1
    if desc.hba > 10:
        violations += 1

    # Veber rules (rotatable bonds <= 10, TPSA <= 140 Å²)
    veber_pass = (desc.rotatable_bonds <= 10) and (desc.tpsa <= 140)

    # Oral absorption estimate (simplified sigmoid model)
    # Based on TPSA – compounds with TPSA < 60 have ~100% absorption
    oral_abs = 1.0 / (1.0 + math.exp((desc.tpsa - 120) / 30))
    # Penalize high MW
    if desc.molecular_weight > 500:
        oral_abs *= 0.7
    # Boost from user-provided absorption tier
    tier_factor = {"High": 1.0, "Moderate": 0.85, "Low": 0.5}.get(
        oral_tier.split("(")[0].strip(), 0.7
    )
    oral_abs = min(1.0, oral_abs * tier_factor)

    # Hepatic clearance penalty – higher logP = more liver metabolism
    clearance = max(0.0, min(1.0, (desc.logp - 1.0) / 6.0))

    # Composite bioavailability = absorption × (1 - clearance) × lipinski/veber bonuses
    bio = oral_abs * (1.0 - 0.5 * clearance)
    if violations == 0:
        bio *= 1.0
    elif violations == 1:
        bio *= 0.85
    else:
        bio *= 0.60
    if veber_pass:
        bio *= 1.0
    else:
        bio *= 0.80
    bio = round(min(1.0, bio), 4)

    # Route recommendation
    if bio >= 0.5:
        route = "oral"
    elif bio >= 0.25:
        route = "injectable"
    else:
        route = "injectable (nanoformulation)"

    return ADMETScore(
        oral_absorption=round(oral_abs, 4),
        lipinski_violations=violations,
        veber_pass=veber_pass,
        bioavailability_score=bio,
        clearance_penalty=round(clearance, 4),
        toxicity_flag=False,
        route_recommendation=route,
    )


# ---------------------------------------------------------------------------
# Compound Library Management
# ---------------------------------------------------------------------------

class CompoundLibrary:
    """Loads and processes the phytochemical library with full descriptor/ADMET profiles."""

    def __init__(self, data_path: Optional[Path] = None):
        if data_path is None:
            data_path = Path(__file__).resolve().parent.parent / "data" / "natural_compounds.json"
        self.data_path = Path(data_path)
        self.compounds: Dict[str, CompoundProfile] = {}
        self._load()

    def _load(self):
        if not self.data_path.exists():
            raise FileNotFoundError(f"Compounds file not found: {self.data_path}")
        with open(self.data_path, "r", encoding="utf-8") as f:
            raw = json.load(f)

        for item in raw.get("compounds", []):
            conc = item.get("target_concentration_range_uM", [1.0, 100.0])

            profile = CompoundProfile(
                id=item["id"],
                name=item["name"],
                plant_source=item["plant_source"],
                smiles=item["smiles"],
                primary_action=item["primary_action"],
                known_mechanism=item["known_mechanism"],
                targets=item["targets"],
                oral_absorption_tier=item.get("oral_absorption_tier", "Moderate"),
                safety_profile=item.get("safety_profile", "Unknown"),
                concentration_range_uM=(conc[0], conc[1]),
            )

            # Compute descriptors
            desc, mol = compute_descriptors(profile.smiles)
            profile.descriptors = desc
            profile.mol = mol

            # Compute ADMET
            if desc is not None:
                profile.admet = assess_admet(desc, profile.oral_absorption_tier)

            self.compounds[profile.id] = profile

    def get(self, compound_id: str) -> Optional[CompoundProfile]:
        return self.compounds.get(compound_id)

    def all_compounds(self) -> List[CompoundProfile]:
        return list(self.compounds.values())

    def filter_by_target(self, target_id: str) -> List[CompoundProfile]:
        """Return compounds that act on a specific target."""
        return [c for c in self.compounds.values() if target_id in c.targets]

    def filter_oral_viable(self, min_bioavailability: float = 0.4) -> List[CompoundProfile]:
        """Return compounds suitable for oral administration."""
        return [
            c for c in self.compounds.values()
            if c.admet is not None and c.admet.bioavailability_score >= min_bioavailability
        ]

    def summary_table(self) -> str:
        """Render a formatted summary table of all compounds."""
        lines = [
            f"{'ID':<8} {'Name':<16} {'MW':>7} {'LogP':>6} {'TPSA':>6} "
            f"{'Lip.V':>5} {'BioAv':>6} {'Route':<14} {'Targets'}"
        ]
        lines.append("-" * 100)
        for c in self.all_compounds():
            d = c.descriptors
            a = c.admet
            if d and a:
                lines.append(
                    f"{c.id:<8} {c.name:<16} {d.molecular_weight:>7.1f} {d.logp:>6.2f} "
                    f"{d.tpsa:>6.1f} {a.lipinski_violations:>5} {a.bioavailability_score:>6.3f} "
                    f"{a.route_recommendation:<14} {', '.join(c.targets)}"
                )
        return "\n".join(lines)
