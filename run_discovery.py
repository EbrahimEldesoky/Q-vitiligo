#!/usr/bin/env python3
"""
QuantumVitiligo Discovery Runner
================================
CLI entrypoint for the quantum-enhanced natural compound discovery pipeline
for autoimmune vitiligo therapeutics.

Usage:
    python run_discovery.py --backend aer --max-iter 30
    python run_discovery.py --backend ibm --max-iter 15 --qubits 4
"""

import argparse
import json
import sys
from pathlib import Path

try:
    from quantum_vitiligo.pipeline import VitiligoDiscoveryPipeline
except ImportError:
    from qvitiligo.pipeline import VitiligoDiscoveryPipeline


def main():
    parser = argparse.ArgumentParser(
        description="QuantumVitiligo: Quantum-Augmented Phytopharmacological Optimization Engine for Vitiligo",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Run with local Aer simulator:
    python run_discovery.py --backend aer --max-iter 30

  Run with IBM Quantum cloud:
    python run_discovery.py --backend ibm --max-iter 15 --qubits 4

  Save results to file:
    python run_discovery.py --backend aer --output results/discovery_report.json
        """,
    )
    parser.add_argument(
        "--target", default="vitiligo",
        help="Disease target (default: vitiligo)",
    )
    parser.add_argument(
        "--backend", choices=["aer", "ibm"], default="aer",
        help="Quantum backend: 'aer' (local simulator) or 'ibm' (IBM Quantum cloud)",
    )
    parser.add_argument(
        "--qubits", type=int, default=4,
        help="Number of qubits for VQE simulation (default: 4)",
    )
    parser.add_argument(
        "--max-iter", type=int, default=30,
        help="Maximum generations for the optimization loop (default: 30)",
    )
    parser.add_argument(
        "--min-bio", type=float, default=0.3,
        help="Minimum bioavailability score for oral prescreening (default: 0.3)",
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Output JSON file path for the discovery report",
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Suppress verbose console output",
    )

    args = parser.parse_args()

    print(r"""
    ==================================================================================
      ___  _   _    _    _   _ _____ _   _ __  __  __     _____ _____ ___ _     ___ ____   ___  
     / _ \| | | |  / \  | \ | |_   _| | | |  \/  | \ \   / /_ _|_   _|_ _| |   |_ _/ ___| / _ \ 
    | | | | | | | / _ \ |  \| | | | | | | | |\/| |  \ \ / / | |  | |  | || |    | | |  _ | | | |
    | |_| | |_| |/ ___ \| |\  | | | | |_| | |  | |   \ V /  | |  | |  | || |___ | | |_| || |_| |
     \__\_\\___//_/   \_\_| \_| |_|  \___/|_|  |_|    \_/  |___| |_| |___|_____|___\____| \___/ 

    Quantum-Augmented Phytopharmacological Discovery Engine for Vitiligo
    Hybrid Classical Chemoinformatics & VQE Quantum Electronic Structure
    ==================================================================================
    """)

    print(f"  Configuration:")
    print(f"    Target Disease:  {args.target}")
    print(f"    Quantum Backend: {args.backend}")
    print(f"    Qubits:          {args.qubits}")
    print(f"    Max Generations: {args.max_iter}")
    print(f"    Min Bioavail:    {args.min_bio}")
    print()

    # Run pipeline
    pipeline = VitiligoDiscoveryPipeline(
        backend=args.backend,
        num_qubits=args.qubits,
    )

    report = pipeline.run(
        max_opt_generations=args.max_iter,
        min_bioavailability=args.min_bio,
        verbose=not args.quiet,
    )

    # Save report
    output_path = args.output
    if output_path is None:
        output_dir = Path(__file__).resolve().parent / "results"
        output_dir.mkdir(exist_ok=True)
        output_path = output_dir / "discovery_report.json"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert numpy arrays for JSON serialization
    def _serialize(obj):
        if hasattr(obj, "tolist"):
            return obj.tolist()
        if hasattr(obj, "__dict__"):
            return str(obj)
        return obj

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=_serialize, ensure_ascii=False)

    print(f"\n  [INFO] Discovery report saved to: {output_path}")
    print(f"  [INFO] Pipeline completed successfully.")

    return 0


if __name__ == "__main__":
    sys.exit(main())

