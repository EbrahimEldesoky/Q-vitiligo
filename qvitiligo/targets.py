"""
Target management and biological pathway mapping for Vitiligo.
Defines receptor profiles for autoimmune arrest (JAK1/2), repigmentation (MC1R/TYR),
and cellular defense (Nrf2/KEAP1).
"""

import json
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Optional

@dataclass
class TargetProfile:
    id: str
    name: str
    gene: str
    uniprot: str
    pdb_id: str
    role: str  # "autoimmune_arrest", "repigmentation", "cellular_defense"
    pathway: str
    mechanism: str
    active_site_residues: List[str]
    binding_weight: float

class TargetRegistry:
    def __init__(self, data_path: Optional[Path] = None):
        if data_path is None:
            data_path = Path(__file__).resolve().parent.parent / "data" / "targets.json"
        self.data_path = Path(data_path)
        self.targets: Dict[str, TargetProfile] = {}
        self._load()

    def _load(self):
        if not self.data_path.exists():
            raise FileNotFoundError(f"Targets file not found: {self.data_path}")
        with open(self.data_path, "r", encoding="utf-8") as f:
            raw = json.load(f)
        for item in raw.get("targets", []):
            profile = TargetProfile(
                id=item["id"],
                name=item["name"],
                gene=item["gene"],
                uniprot=item["uniprot"],
                pdb_id=item["pdb_id"],
                role=item["role"],
                pathway=item["pathway"],
                mechanism=item["mechanism"],
                active_site_residues=item["active_site_residues"],
                binding_weight=float(item.get("binding_weight", 0.2))
            )
            self.targets[profile.id] = profile

    def get(self, target_id: str) -> Optional[TargetProfile]:
        return self.targets.get(target_id)

    def get_by_role(self, role: str) -> List[TargetProfile]:
        return [t for t in self.targets.values() if t.role == role]

    def all_targets(self) -> List[TargetProfile]:
        return list(self.targets.values())
