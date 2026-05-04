#!/usr/bin/env python3
"""Protein structure analysis toolkit.

Implements:
- Primary structure analysis
- Secondary structure analysis
- Tertiary structure analysis
- Quaternary structure analysis
- Binding pocket analysis
- Tunnel prediction

Usage:
  python protein_structure_analysis.py --pdb input.pdb --out analysis.json
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Polypeptide import is_aa
from scipy.spatial import cKDTree


AA_MASS = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10, "C": 121.16,
    "Q": 146.15, "E": 147.13, "G": 75.07, "H": 155.16, "I": 131.17,
    "L": 131.17, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13,
    "S": 105.09, "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
}


@dataclass
class Pocket:
    center: Tuple[float, float, float]
    volume: float
    points: int
    depth_score: float


def read_structure(path: str):
    parser = PDBParser(QUIET=True)
    return parser.get_structure("protein", path)


def atom_coordinates(structure) -> np.ndarray:
    coords = [atom.coord for atom in structure.get_atoms() if atom.element != "H"]
    return np.array(coords, dtype=np.float64)


def extract_sequences(structure) -> Dict[str, str]:
    ppb = PPBuilder()
    seqs: Dict[str, str] = {}
    for model in structure:
        for chain in model:
            peptides = ppb.build_peptides(chain)
            seq = "".join(str(p.get_sequence()) for p in peptides)
            if seq:
                seqs[chain.id] = seq
        break
    return seqs


def analyze_primary(sequences: Dict[str, str]) -> Dict:
    chain_results = {}
    total_len = 0
    combined = ""
    for chain, seq in sequences.items():
        comp = {aa: seq.count(aa) for aa in AA_MASS}
        mw = sum(AA_MASS.get(aa, 0.0) for aa in seq) - 18.015 * max(len(seq) - 1, 0)
        chain_results[chain] = {
            "length": len(seq),
            "molecular_weight_da": round(mw, 2),
            "composition": comp,
        }
        total_len += len(seq)
        combined += seq

    hydrophobic = set("AVILMFWY")
    hydrophobic_frac = round(sum(1 for aa in combined if aa in hydrophobic) / max(len(combined), 1), 4)
    return {
        "chains": chain_results,
        "total_residues": total_len,
        "hydrophobic_fraction": hydrophobic_frac,
    }


def classify_secondary(phi: float | None, psi: float | None) -> str:
    if phi is None or psi is None:
        return "coil"
    if -160 <= phi <= -20 and -80 <= psi <= -10:
        return "helix"
    if -180 <= phi <= -60 and 90 <= psi <= 180:
        return "sheet"
    return "coil"


def analyze_secondary(structure) -> Dict:
    ppb = PPBuilder()
    per_chain = {}
    for model in structure:
        for chain in model:
            helix = sheet = coil = 0
            for peptide in ppb.build_peptides(chain):
                phipsi = peptide.get_phi_psi_list()
                for phi, psi in phipsi:
                    cls = classify_secondary(phi, psi)
                    if cls == "helix":
                        helix += 1
                    elif cls == "sheet":
                        sheet += 1
                    else:
                        coil += 1
            total = helix + sheet + coil
            if total > 0:
                per_chain[chain.id] = {
                    "helix": helix,
                    "sheet": sheet,
                    "coil": coil,
                    "fractions": {
                        "helix": round(helix / total, 4),
                        "sheet": round(sheet / total, 4),
                        "coil": round(coil / total, 4),
                    },
                }
        break
    return {"method": "phi/psi heuristic", "chains": per_chain}


def analyze_tertiary(structure, coords: np.ndarray) -> Dict:
    centroid = coords.mean(axis=0)
    rg = math.sqrt(np.mean(np.sum((coords - centroid) ** 2, axis=1)))

    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=8.0)
    contacts = len(pairs)

    residues = [res for res in structure.get_residues() if is_aa(res, standard=True)]
    ca_count = sum(1 for res in residues if "CA" in res)

    return {
        "radius_of_gyration": round(float(rg), 4),
        "compactness_contacts_8A": contacts,
        "residue_count": len(residues),
        "ca_atoms": ca_count,
    }


def analyze_quaternary(structure) -> Dict:
    first_model = next(structure.get_models())
    chain_coords = {}
    for chain in first_model:
        coords = [atom.coord for atom in chain.get_atoms() if atom.element != "H"]
        if coords:
            chain_coords[chain.id] = np.array(coords)

    interfaces = []
    ids = sorted(chain_coords.keys())
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            a, b = ids[i], ids[j]
            tree = cKDTree(chain_coords[a])
            neighbors = tree.query_ball_point(chain_coords[b], r=5.0)
            contacts = sum(1 for n in neighbors if n)
            if contacts > 0:
                interfaces.append({"chain_a": a, "chain_b": b, "contacts_5A": contacts})

    return {
        "chains": ids,
        "oligomeric_state_guess": f"{len(ids)}-mer" if ids else "unknown",
        "interfaces": interfaces,
    }


def _grid_from_coords(coords: np.ndarray, spacing: float = 1.5, padding: float = 6.0):
    mins = coords.min(axis=0) - padding
    maxs = coords.max(axis=0) + padding
    grid_axes = [np.arange(mins[i], maxs[i] + spacing, spacing) for i in range(3)]
    X, Y, Z = np.meshgrid(*grid_axes, indexing="ij")
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    return points, grid_axes, mins, spacing


def detect_pockets(coords: np.ndarray, spacing: float = 1.5) -> List[Pocket]:
    points, _, _, _ = _grid_from_coords(coords, spacing=spacing)
    tree = cKDTree(coords)
    dists, _ = tree.query(points, k=1)

    candidate_mask = (dists > 2.0) & (dists < 6.5)
    candidates = points[candidate_mask]
    if len(candidates) == 0:
        return []

    # cluster by proximity using KD graph bfs
    ctree = cKDTree(candidates)
    seen = np.zeros(len(candidates), dtype=bool)
    pockets: List[Pocket] = []
    neigh_r = 2.2

    for i in range(len(candidates)):
        if seen[i]:
            continue
        stack = [i]
        seen[i] = True
        cluster_idx = []
        while stack:
            idx = stack.pop()
            cluster_idx.append(idx)
            nbrs = ctree.query_ball_point(candidates[idx], r=neigh_r)
            for n in nbrs:
                if not seen[n]:
                    seen[n] = True
                    stack.append(n)
        if len(cluster_idx) < 20:
            continue
        cluster = candidates[cluster_idx]
        cluster_d = dists[candidate_mask][cluster_idx]
        volume = len(cluster) * (spacing ** 3)
        pockets.append(Pocket(tuple(cluster.mean(axis=0)), float(volume), len(cluster), float(cluster_d.mean())))

    pockets.sort(key=lambda p: (p.volume, p.depth_score), reverse=True)
    return pockets[:8]


def predict_tunnels(coords: np.ndarray, pockets: List[Pocket]) -> List[Dict]:
    if not pockets:
        return []
    protein_center = coords.mean(axis=0)
    tunnels = []
    for p in pockets[:5]:
        direction = np.array(p.center) - protein_center
        norm = np.linalg.norm(direction)
        if norm < 1e-6:
            continue
        unit = direction / norm
        est_length = max(norm - 2.0, 0.0)
        est_bottleneck = max(1.0, 3.5 - p.depth_score * 0.3)
        tunnels.append({
            "from_pocket_center": [round(v, 3) for v in p.center],
            "direction": [round(float(v), 3) for v in unit],
            "estimated_length": round(float(est_length), 3),
            "estimated_bottleneck_radius": round(float(est_bottleneck), 3),
        })
    return tunnels


def run_analysis(pdb_path: str) -> Dict:
    structure = read_structure(pdb_path)
    coords = atom_coordinates(structure)
    sequences = extract_sequences(structure)

    primary = analyze_primary(sequences)
    secondary = analyze_secondary(structure)
    tertiary = analyze_tertiary(structure, coords)
    quaternary = analyze_quaternary(structure)

    pockets = detect_pockets(coords)
    pocket_json = [
        {
            "center": [round(float(v), 3) for v in p.center],
            "volume_A3": round(p.volume, 3),
            "support_points": p.points,
            "depth_score": round(p.depth_score, 3),
        }
        for p in pockets
    ]

    tunnels = predict_tunnels(coords, pockets)

    return {
        "input_pdb": str(pdb_path),
        "primary_structure": primary,
        "secondary_structure": secondary,
        "tertiary_structure": tertiary,
        "quaternary_structure": quaternary,
        "binding_pockets": {
            "method": "grid cavity heuristic",
            "pockets": pocket_json,
        },
        "tunnel_prediction": {
            "method": "pocket-to-surface vector heuristic",
            "tunnels": tunnels,
        },
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True, help="Input PDB file")
    ap.add_argument("--out", default="analysis.json", help="Output JSON path")
    args = ap.parse_args()

    result = run_analysis(args.pdb)
    Path(args.out).write_text(json.dumps(result, indent=2))
    print(f"Saved analysis to {args.out}")


if __name__ == "__main__":
    main()
