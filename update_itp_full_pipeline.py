#!/usr/bin/env python3
"""
update_itp_trustable_dihe_json.py

End-to-end HessFit + QM/MM torsion-scan updater for GROMACS .itp files.

This script starts from the ORIGINAL .itp and performs two stages:

STAGE A — HessFit base update
  1. Update [ atoms ] charges from charge.txt / type_charge.txt.
  2. Update [ bonds ] from HessFit HrmStr1 terms in ff_string.txt.
  3. Update [ angles ] from HessFit HrmBnd1 terms in ff_string.txt.
  4. Update proper [ dihedrals ] from HessFit AmbTrs terms in ff_string.txt.
  5. Keep nonbonded atom types, [ pairs ], and improper dihedrals unchanged.

STAGE B — QM/MM torsion-scan refinement
  1. Read authoritative scan_torsions from dihe_optfile.json.
  2. Map scan index n to n_qm_scan_energy*.csv and n_mm_scan_energy*.csv.
  3. Find the matching proper [ dihedrals ] entry in the Stage-A .itp.
  4. Pair QM and MM scan points by inferred scan-step index, not blindly by row.
  5. Fit QM_relative - MM_relative to a periodic torsion term.
  6. Override only trusted selected proper dihedrals.

Conversion assumptions implemented in Stage A:
  Bonds:
    HessFit HrmStr1 K [kcal mol^-1 A^-2], r0 [A]
    -> GROMACS/GROMOS bond funct 2
       r0_nm = r0_A / 10
       k = K * 4.184 * 100 / r0_nm^2

  Angles:
    HessFit HrmBnd1 K [kcal mol^-1 rad^-2], theta0 [degree]
    -> GROMACS/GROMOS angle funct 2
       k = 2 * K * 4.184 / sin(theta0)^2

  Proper dihedrals:
    HessFit AmbTrs barrier [kcal mol^-1]
    -> GROMACS proper dihedral funct 1
       cp [kJ mol^-1] = barrier * 4.184
       phase and multiplicity are read from the AmbTrs record.

Important trust features:
  - Validates charge, ff_string, scan JSON, scan CSVs, and .itp mappings.
  - Diagnoses missing files, missing points, scan mismatches, high RMSE, and large torsion k.
  - Uses angle-aware QM/MM pairing in relaxed incomplete-scan mode.
  - Selects the most complete available scan-energy CSV when multiple versions exist.
  - Can abort unless Stage A is complete and/or all scan refinements are trusted.

Typical run:
  python update_itp_trustable_dihe_json.py \
      --data ./data \
      --base-itp BK7T.itp \
      --dihe-json dihe_optfile.json \
      --charges type_charge.txt \
      --ff-string ff_string.txt \
      --output-itp BK7T_full_hessfit_scan_refined.itp

Strict final run:
  python update_itp_trustable_dihe_json.py \
      --data ./data \
      --base-itp BK7T.itp \
      --dihe-json dihe_optfile.json \
      --charges type_charge.txt \
      --ff-string ff_string.txt \
      --output-itp BK7T_full_hessfit_scan_refined.itp \
      --require-complete-stage-a \
      --require-all-scans

Relaxed scan refinement allowing up to three missing paired scan points:
  python update_itp_trustable_dihe_json.py \
      --data ./data \
      --base-itp BK7T.itp \
      --dihe-json dihe_optfile.json \
      --charges type_charge.txt \
      --ff-string ff_string.txt \
      --output-itp BK7T_full_hessfit_scan_refined_relaxed.itp \
      --allow-incomplete \
      --expected-points 8 \
      --rmse-max-kj 20
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
from collections import defaultdict
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


HARTREE_TO_KCAL = 627.5094740631
KCAL_TO_KJ = 4.184


# =============================================================================
# Data structures
# =============================================================================

@dataclass
class ChargeEntry:
    label: str
    atom_index: int  # 1-based .itp atom index
    charge: float


@dataclass
class BondTerm:
    labels: Tuple[str, str]
    atoms: Tuple[int, int]
    k_kcal_a2: float
    r0_a: float
    r0_nm: float
    k_gmx: float


@dataclass
class AngleTerm:
    labels: Tuple[str, str, str]
    atoms: Tuple[int, int, int]
    k_kcal_rad2: float
    theta0_deg: float
    k_gmx: float


@dataclass
class DihedralTerm:
    labels: Tuple[str, str, str, str]
    atoms: Tuple[int, int, int, int]
    phase_deg: float
    barrier_kcal: float
    cp_kj: float
    multiplicity: int
    raw_numeric: Tuple[float, ...]


@dataclass
class StageASummary:
    atoms_in_itp: int = 0
    charges_parsed: int = 0
    charges_updated: int = 0
    charge_missing_atom_indices: List[int] = field(default_factory=list)
    charge_duplicate_atom_indices: List[int] = field(default_factory=list)
    charge_invalid_labels: List[str] = field(default_factory=list)

    bond_terms_parsed: int = 0
    bond_lines_in_itp: int = 0
    bond_lines_updated: int = 0
    bond_itp_unmatched: List[str] = field(default_factory=list)
    bond_ff_unmatched: List[str] = field(default_factory=list)
    bond_duplicate_ff_keys: List[str] = field(default_factory=list)
    bond_skipped_wrong_funct: List[str] = field(default_factory=list)

    angle_terms_parsed: int = 0
    angle_lines_in_itp: int = 0
    angle_lines_updated: int = 0
    angle_itp_unmatched: List[str] = field(default_factory=list)
    angle_ff_unmatched: List[str] = field(default_factory=list)
    angle_duplicate_ff_keys: List[str] = field(default_factory=list)
    angle_skipped_wrong_funct: List[str] = field(default_factory=list)

    dihedral_terms_parsed: int = 0
    proper_dihedral_lines_in_itp: int = 0
    proper_dihedral_lines_updated: int = 0
    dihedral_itp_unmatched: List[str] = field(default_factory=list)
    dihedral_ff_unmatched: List[str] = field(default_factory=list)
    dihedral_duplicate_ff_keys: List[str] = field(default_factory=list)
    dihedral_skipped_improper: List[str] = field(default_factory=list)
    dihedral_skipped_wrong_funct: List[str] = field(default_factory=list)

    ff_parse_warnings: List[str] = field(default_factory=list)


@dataclass
class ScanCandidate:
    scan_idx: int
    torsion: Tuple[int, int, int, int]
    central_bond: Tuple[int, int]
    stage_a_phase_deg: float = 180.0
    stage_a_cp_kj: float = 0.0
    stage_a_mult: int = 2
    itp_match_count: int = 0
    itp_match_lines: Tuple[int, ...] = ()
    mapping_issues: List[str] = field(default_factory=list)


@dataclass
class PointPair:
    qm_row: int
    scan_step: int
    qm_angle_deg: float
    expected_angle_deg: float
    angle_error_deg: float
    qm_rel_kcal: float
    mm_rel_kcal: float
    target_kcal: float
    mm_source_row: int


@dataclass
class ScanFit:
    scan_idx: int
    torsion: Tuple[int, int, int, int]
    central_bond: Tuple[int, int]
    qm_file: Optional[str]
    mm_file: Optional[str]
    qm_points: int
    mm_points: int
    paired_points: int
    trusted: bool
    skip_reasons: List[str]
    multiplicity: int
    phase_deg: Optional[float]
    cp_kj: Optional[float]
    cp_kcal: Optional[float]
    rmse_kj: Optional[float]
    rmse_kcal: Optional[float]
    stage_a_phase_deg: float
    stage_a_cp_kj: float
    stage_a_mult: int
    itp_match_count: int
    itp_match_lines: Tuple[int, ...]
    max_angle_error_deg: Optional[float]
    pairing_warnings: List[str]
    log_diagnostics: Dict[str, Any]


# =============================================================================
# General helpers
# =============================================================================

def read_lines(path: Path) -> List[str]:
    return path.read_text(errors="replace").splitlines()


def strip_comment(line: str) -> str:
    return line.split(";", 1)[0]


def existing_comment(line: str) -> str:
    if ";" not in line:
        return ""
    return ";" + line.split(";", 1)[1]


def safe_float(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def version_score(path: Path) -> Tuple[int, float]:
    m = re.search(r"\((\d+)\)(?=\.[^.]+$)", path.name)
    version = int(m.group(1)) if m else 0
    try:
        mtime = path.stat().st_mtime
    except OSError:
        mtime = 0.0
    return version, mtime


def label_to_itp_index(label: str) -> Optional[int]:
    """
    HessFit labels observed in this workflow are H0, C1, N12, ...
    Numeric suffix is zero-based; .itp atom numbering is one-based.
    """
    m = re.search(r"(\d+)$", label.strip())
    if not m:
        return None
    return int(m.group(1)) + 1


def key_bond(atoms: Sequence[int]) -> Tuple[int, int]:
    i, j = (int(x) for x in atoms)
    return tuple(sorted((i, j)))  # type: ignore[return-value]


def key_angle(atoms: Sequence[int]) -> Tuple[int, int, int]:
    a = tuple(int(x) for x in atoms)
    rev = tuple(reversed(a))
    return min(a, rev)  # type: ignore[return-value]


def key_dihedral(atoms: Sequence[int]) -> Tuple[int, int, int, int]:
    a = tuple(int(x) for x in atoms)
    rev = tuple(reversed(a))
    return min(a, rev)  # type: ignore[return-value]


def circular_angle_deg(angle: float) -> float:
    x = (angle + 180.0) % 360.0 - 180.0
    # Preserve +180 if the input is close to +180 rather than -180.
    if abs(x + 180.0) < 1e-10 and angle > 0:
        return 180.0
    return x


def circular_distance_deg(a: float, b: float) -> float:
    return abs(circular_angle_deg(a - b))


def replace_line(lines: List[str], line_no: int, new_line: str) -> None:
    lines[line_no - 1] = new_line


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    keys = sorted({key for row in rows for key in row.keys()})
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def resolve_optional_file(data_dir: Path, value: Optional[str], label: str) -> Optional[Path]:
    if not value:
        return None
    raw = Path(value)
    path = raw if raw.is_absolute() else data_dir / raw
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def choose_first_existing(paths: Iterable[Path]) -> Optional[Path]:
    existing = [p for p in paths if p.exists()]
    if not existing:
        return None
    return sorted(existing, key=version_score, reverse=True)[0]


# =============================================================================
# Input discovery
# =============================================================================

def discover_base_itp(data_dir: Path, explicit: Optional[str]) -> Path:
    p = resolve_optional_file(data_dir, explicit, "Base .itp")
    if p:
        return p
    priority = ["BK7T.itp", "LBAI.itp", "LBAI(2).itp"]
    p = choose_first_existing(data_dir / name for name in priority)
    if p:
        return p
    hits = sorted(data_dir.glob("*.itp"), key=lambda x: x.stat().st_mtime, reverse=True)
    if not hits:
        raise FileNotFoundError(f"No .itp file found in {data_dir}")
    return hits[0]


def discover_dihe_json(data_dir: Path, explicit: Optional[str]) -> Path:
    p = resolve_optional_file(data_dir, explicit, "dihe JSON")
    if p:
        return p
    hits = list(data_dir.glob("dihe_optfile*.json"))
    if not hits:
        raise FileNotFoundError(f"No dihe_optfile*.json found in {data_dir}")
    return sorted(hits, key=version_score, reverse=True)[0]


def discover_aux_file_from_json_or_glob(
    data_dir: Path,
    explicit: Optional[str],
    label: str,
    json_files: Dict[str, Any],
    json_keys: Sequence[str],
    glob_patterns: Sequence[str],
) -> Path:
    p = resolve_optional_file(data_dir, explicit, label)
    if p:
        return p

    for key in json_keys:
        val = json_files.get(key)
        if isinstance(val, str) and val.strip():
            raw = Path(val)
            candidates = [raw] if raw.is_absolute() else [data_dir / raw, data_dir / raw.name]
            p = choose_first_existing(candidates)
            if p:
                return p

    hits: List[Path] = []
    for pattern in glob_patterns:
        hits.extend(data_dir.glob(pattern))
    if not hits:
        raise FileNotFoundError(f"Could not locate {label} in {data_dir}")
    return sorted(hits, key=version_score, reverse=True)[0]


def numeric_row_count(path: Path) -> int:
    count = 0
    for raw in read_lines(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        toks = [t for t in re.split(r"[,\s]+", line) if t]
        if not toks:
            continue
        try:
            [safe_float(t) for t in toks]
            count += 1
        except ValueError:
            continue
    return count


def choose_scan_energy_file(data_dir: Path, scan_idx: int, kind: str, expected_full_points: int) -> Optional[Path]:
    """
    Safer than choosing only by parenthesized version suffix.
    Ranking:
      1. Prefer files with exactly expected_full_points rows.
      2. Then prefer files closer to expected_full_points without exceeding a tie penalty.
      3. Then prefer more numeric rows.
      4. Then prefer higher parenthesized version / mtime.
    """
    hits = list(data_dir.glob(f"{scan_idx}_{kind}_scan_energy*.csv"))
    if not hits:
        return None

    def score(path: Path) -> Tuple[int, int, int, int, float]:
        rows = numeric_row_count(path)
        exact = int(rows == expected_full_points)
        closeness = -abs(expected_full_points - rows)
        not_overlarge = int(rows <= expected_full_points + 2)
        version, mtime = version_score(path)
        return (exact, not_overlarge, closeness, rows, version + mtime * 1e-12)

    return sorted(hits, key=score, reverse=True)[0]


# =============================================================================
# .itp parsing
# =============================================================================

def parse_sections(lines: List[str]) -> Dict[str, List[Tuple[int, str]]]:
    sections: Dict[str, List[Tuple[int, str]]] = defaultdict(list)
    current: Optional[str] = None
    for line_no, line in enumerate(lines, start=1):
        s = line.strip()
        if s.startswith("[") and "]" in s:
            current = s[s.find("[") + 1:s.find("]")].strip().lower()
            continue
        if current:
            sections[current].append((line_no, line))
    return dict(sections)


def parse_atoms(lines: List[str]) -> List[Dict[str, Any]]:
    atoms: List[Dict[str, Any]] = []
    for line_no, line in parse_sections(lines).get("atoms", []):
        toks = strip_comment(line).split()
        if len(toks) < 7:
            continue
        try:
            atoms.append(
                {
                    "line_no": line_no,
                    "nr": int(toks[0]),
                    "type": toks[1],
                    "resnr": toks[2],
                    "residue": toks[3],
                    "atom": toks[4],
                    "cgnr": toks[5],
                    "charge": safe_float(toks[6]),
                    "mass": safe_float(toks[7]) if len(toks) > 7 else None,
                    "tokens": toks,
                    "line": line,
                }
            )
        except Exception:
            continue
    return sorted(atoms, key=lambda x: x["nr"])


def parse_bonds(lines: List[str]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for line_no, line in parse_sections(lines).get("bonds", []):
        toks = strip_comment(line).split()
        if len(toks) < 5:
            continue
        try:
            atoms = (int(toks[0]), int(toks[1]))
            out.append(
                {
                    "line_no": line_no,
                    "atoms": atoms,
                    "key": key_bond(atoms),
                    "funct": int(float(toks[2])),
                    "r0_nm": safe_float(toks[3]),
                    "k": safe_float(toks[4]),
                    "tokens": toks,
                    "line": line,
                }
            )
        except Exception:
            continue
    return out


def parse_angles(lines: List[str]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for line_no, line in parse_sections(lines).get("angles", []):
        toks = strip_comment(line).split()
        if len(toks) < 6:
            continue
        try:
            atoms = (int(toks[0]), int(toks[1]), int(toks[2]))
            out.append(
                {
                    "line_no": line_no,
                    "atoms": atoms,
                    "key": key_angle(atoms),
                    "funct": int(float(toks[3])),
                    "theta0_deg": safe_float(toks[4]),
                    "k": safe_float(toks[5]),
                    "tokens": toks,
                    "line": line,
                }
            )
        except Exception:
            continue
    return out


def is_dihedral_tokens(toks: List[str]) -> bool:
    if len(toks) < 8:
        return False
    try:
        int(toks[0]); int(toks[1]); int(toks[2]); int(toks[3])
        int(float(toks[4])); safe_float(toks[5]); safe_float(toks[6]); int(float(toks[7]))
        return True
    except Exception:
        return False


def parse_dihedrals(lines: List[str]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for line_no, line in parse_sections(lines).get("dihedrals", []):
        toks = strip_comment(line).split()
        if not is_dihedral_tokens(toks):
            continue
        try:
            atoms = (int(toks[0]), int(toks[1]), int(toks[2]), int(toks[3]))
            out.append(
                {
                    "line_no": line_no,
                    "atoms": atoms,
                    "key": key_dihedral(atoms),
                    "funct": int(float(toks[4])),
                    "phase_deg": safe_float(toks[5]),
                    "cp_kj": safe_float(toks[6]),
                    "mult": int(float(toks[7])),
                    "tokens": toks,
                    "line": line,
                }
            )
        except Exception:
            continue
    return out


# =============================================================================
# Charge and ff_string parsing
# =============================================================================

def parse_charges(path: Path) -> Tuple[Dict[int, ChargeEntry], List[str], List[int]]:
    entries: Dict[int, ChargeEntry] = {}
    invalid_labels: List[str] = []
    duplicates: List[int] = []
    for raw in read_lines(path):
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith(";"):
            continue
        toks = line.split()
        if len(toks) < 2:
            continue
        label = toks[0]
        idx = label_to_itp_index(label)
        if idx is None:
            invalid_labels.append(label)
            continue
        try:
            q = safe_float(toks[1])
        except ValueError:
            invalid_labels.append(label)
            continue
        if idx in entries:
            duplicates.append(idx)
        entries[idx] = ChargeEntry(label=label, atom_index=idx, charge=q)
    return entries, invalid_labels, sorted(set(duplicates))


def parse_ff_string(path: Path) -> Tuple[List[BondTerm], List[AngleTerm], List[DihedralTerm], List[str]]:
    bonds: List[BondTerm] = []
    angles: List[AngleTerm] = []
    dihedrals: List[DihedralTerm] = []
    warnings: List[str] = []

    for raw in read_lines(path):
        line = raw.strip()
        if not line or line.startswith("!") or line.startswith("#") or line.startswith(";"):
            continue
        toks = line.split()
        kind = toks[0]
        try:
            if kind == "HrmStr1":
                if len(toks) < 5:
                    warnings.append(f"Short HrmStr1 line: {line}")
                    continue
                labels = (toks[1], toks[2])
                ids = [label_to_itp_index(x) for x in labels]
                if any(x is None for x in ids):
                    warnings.append(f"Could not map HrmStr1 labels: {line}")
                    continue
                k_kcal = safe_float(toks[3])
                r0_a = safe_float(toks[4])
                r0_nm = r0_a / 10.0
                if r0_nm <= 0:
                    warnings.append(f"Invalid HrmStr1 r0: {line}")
                    continue
                k_gmx = k_kcal * KCAL_TO_KJ * 100.0 / (r0_nm ** 2)
                bonds.append(
                    BondTerm(
                        labels=labels,
                        atoms=(int(ids[0]), int(ids[1])),
                        k_kcal_a2=k_kcal,
                        r0_a=r0_a,
                        r0_nm=r0_nm,
                        k_gmx=k_gmx,
                    )
                )

            elif kind == "HrmBnd1":
                if len(toks) < 6:
                    warnings.append(f"Short HrmBnd1 line: {line}")
                    continue
                labels = (toks[1], toks[2], toks[3])
                ids = [label_to_itp_index(x) for x in labels]
                if any(x is None for x in ids):
                    warnings.append(f"Could not map HrmBnd1 labels: {line}")
                    continue
                k_kcal = safe_float(toks[4])
                theta0 = safe_float(toks[5])
                s = math.sin(math.radians(theta0))
                if abs(s) < 1e-12:
                    warnings.append(f"sin(theta0) too small in HrmBnd1: {line}")
                    continue
                k_gmx = 2.0 * k_kcal * KCAL_TO_KJ / (s ** 2)
                angles.append(
                    AngleTerm(
                        labels=labels,
                        atoms=(int(ids[0]), int(ids[1]), int(ids[2])),
                        k_kcal_rad2=k_kcal,
                        theta0_deg=theta0,
                        k_gmx=k_gmx,
                    )
                )

            elif kind == "AmbTrs":
                if len(toks) < 14:
                    warnings.append(f"Short AmbTrs line: {line}")
                    continue
                labels = (toks[1], toks[2], toks[3], toks[4])
                ids = [label_to_itp_index(x) for x in labels]
                if any(x is None for x in ids):
                    warnings.append(f"Could not map AmbTrs labels: {line}")
                    continue
                numeric = tuple(safe_float(x) for x in toks[5:])
                # HessFit format observed:
                # labels + [0, phase, 0, 0, 0.00, barrier, 0.00, 0., multiplicity]
                if len(numeric) < 9:
                    warnings.append(f"AmbTrs numeric vector too short: {line}")
                    continue
                phase = numeric[1]
                barrier = numeric[5]
                mult = int(round(numeric[-1]))
                dihedrals.append(
                    DihedralTerm(
                        labels=labels,
                        atoms=(int(ids[0]), int(ids[1]), int(ids[2]), int(ids[3])),
                        phase_deg=phase,
                        barrier_kcal=barrier,
                        cp_kj=barrier * KCAL_TO_KJ,
                        multiplicity=mult,
                        raw_numeric=numeric,
                    )
                )
        except Exception as exc:
            warnings.append(f"Failed to parse ff_string line: {line} | {exc}")

    return bonds, angles, dihedrals, warnings


# =============================================================================
# Stage A formatting and update
# =============================================================================

def format_atom_line(atom: Dict[str, Any], new_charge: float) -> str:
    toks = list(atom["tokens"])
    toks[6] = f"{new_charge:.6f}"
    mass_txt = f" {safe_float(toks[7]):12.6f}" if len(toks) > 7 else ""
    comment = existing_comment(atom["line"])
    core = (
        f"{int(toks[0]):6d} {toks[1]:<12s} {toks[2]:>6s} {toks[3]:<8s} "
        f"{toks[4]:<8s} {toks[5]:>6s} {new_charge:12.6f}{mass_txt}"
    )
    if comment:
        core += "   " + comment
    return core


def format_bond_line(entry: Dict[str, Any], term: BondTerm) -> str:
    i, j = entry["atoms"]
    comment = existing_comment(entry["line"])
    core = f"{i:6d}{j:6d}{2:6d}{term.r0_nm:14.8f}{term.k_gmx:16.6f}"
    if comment:
        core += "   " + comment
    core += f"   ; HessFit HrmStr1 {'-'.join(term.labels)}"
    return core


def format_angle_line(entry: Dict[str, Any], term: AngleTerm) -> str:
    i, j, k = entry["atoms"]
    comment = existing_comment(entry["line"])
    core = f"{i:6d}{j:6d}{k:6d}{2:6d}{term.theta0_deg:14.8f}{term.k_gmx:16.6f}"
    if comment:
        core += "   " + comment
    core += f"   ; HessFit HrmBnd1 {'-'.join(term.labels)}"
    return core


def format_stage_a_dihedral_line(entry: Dict[str, Any], term: DihedralTerm) -> str:
    i, j, k, l = entry["atoms"]
    comment = existing_comment(entry["line"])
    core = (
        f"{i:6d}{j:6d}{k:6d}{l:6d}{1:6d}"
        f"{term.phase_deg:12.4f}{term.cp_kj:12.5f}{term.multiplicity:6d}"
    )
    if comment:
        core += "   " + comment
    core += f"   ; HessFit AmbTrs {'-'.join(term.labels)}"
    return core


def apply_stage_a(
    original_lines: List[str],
    charges: Dict[int, ChargeEntry],
    charge_invalid_labels: List[str],
    charge_duplicates: List[int],
    bond_terms: List[BondTerm],
    angle_terms: List[AngleTerm],
    dihedral_terms: List[DihedralTerm],
    ff_warnings: List[str],
) -> Tuple[List[str], StageASummary, List[Dict[str, Any]]]:
    out = list(original_lines)
    summary = StageASummary()
    changes: List[Dict[str, Any]] = []

    atoms = parse_atoms(out)
    summary.atoms_in_itp = len(atoms)
    summary.charges_parsed = len(charges)
    summary.charge_invalid_labels = list(charge_invalid_labels)
    summary.charge_duplicate_atom_indices = list(charge_duplicates)

    atom_indices = {a["nr"] for a in atoms}
    summary.charge_missing_atom_indices = sorted(atom_indices - set(charges))

    for atom in atoms:
        idx = atom["nr"]
        if idx not in charges:
            continue
        q = charges[idx]
        replace_line(out, atom["line_no"], format_atom_line(atom, q.charge))
        summary.charges_updated += 1
        changes.append(
            {
                "stage": "A_charge",
                "line_no": atom["line_no"],
                "atom_index": idx,
                "hessfit_label": q.label,
                "old_charge": atom["charge"],
                "new_charge": q.charge,
            }
        )

    bonds = parse_bonds(out)
    angles = parse_angles(out)
    dihedrals = parse_dihedrals(out)

    # Bonds
    summary.bond_terms_parsed = len(bond_terms)
    summary.bond_lines_in_itp = len(bonds)
    bond_map: Dict[Tuple[int, int], List[BondTerm]] = defaultdict(list)
    for term in bond_terms:
        bond_map[key_bond(term.atoms)].append(term)
    summary.bond_duplicate_ff_keys = [
        " ".join(map(str, key)) for key, vals in bond_map.items() if len(vals) > 1
    ]
    itp_bond_keys = {b["key"] for b in bonds}
    for bond in bonds:
        terms = bond_map.get(bond["key"], [])
        if not terms:
            summary.bond_itp_unmatched.append(
                f"line {bond['line_no']}: {' '.join(map(str, bond['atoms']))}"
            )
            continue
        if bond["funct"] != 2:
            summary.bond_skipped_wrong_funct.append(
                f"line {bond['line_no']}: funct={bond['funct']}"
            )
            continue
        term = terms[0]
        replace_line(out, bond["line_no"], format_bond_line(bond, term))
        summary.bond_lines_updated += 1
        changes.append(
            {
                "stage": "A_bond",
                "line_no": bond["line_no"],
                "atoms": " ".join(map(str, bond["atoms"])),
                "old_r0_nm": bond["r0_nm"],
                "new_r0_nm": term.r0_nm,
                "old_k": bond["k"],
                "new_k": term.k_gmx,
            }
        )
    for key, terms in bond_map.items():
        if key not in itp_bond_keys:
            summary.bond_ff_unmatched.append(
                f"{' '.join(map(str, key))} ({'-'.join(terms[0].labels)})"
            )

    # Angles
    summary.angle_terms_parsed = len(angle_terms)
    summary.angle_lines_in_itp = len(angles)
    angle_map: Dict[Tuple[int, int, int], List[AngleTerm]] = defaultdict(list)
    for term in angle_terms:
        angle_map[key_angle(term.atoms)].append(term)
    summary.angle_duplicate_ff_keys = [
        " ".join(map(str, key)) for key, vals in angle_map.items() if len(vals) > 1
    ]
    itp_angle_keys = {a["key"] for a in angles}
    for angle in angles:
        terms = angle_map.get(angle["key"], [])
        if not terms:
            summary.angle_itp_unmatched.append(
                f"line {angle['line_no']}: {' '.join(map(str, angle['atoms']))}"
            )
            continue
        if angle["funct"] != 2:
            summary.angle_skipped_wrong_funct.append(
                f"line {angle['line_no']}: funct={angle['funct']}"
            )
            continue
        term = terms[0]
        replace_line(out, angle["line_no"], format_angle_line(angle, term))
        summary.angle_lines_updated += 1
        changes.append(
            {
                "stage": "A_angle",
                "line_no": angle["line_no"],
                "atoms": " ".join(map(str, angle["atoms"])),
                "old_theta0_deg": angle["theta0_deg"],
                "new_theta0_deg": term.theta0_deg,
                "old_k": angle["k"],
                "new_k": term.k_gmx,
            }
        )
    for key, terms in angle_map.items():
        if key not in itp_angle_keys:
            summary.angle_ff_unmatched.append(
                f"{' '.join(map(str, key))} ({'-'.join(terms[0].labels)})"
            )

    # Proper dihedrals
    summary.dihedral_terms_parsed = len(dihedral_terms)
    proper_itp_dihedrals = [d for d in dihedrals if d["funct"] != 2]
    summary.proper_dihedral_lines_in_itp = len(proper_itp_dihedrals)
    dih_map: Dict[Tuple[int, int, int, int], List[DihedralTerm]] = defaultdict(list)
    for term in dihedral_terms:
        dih_map[key_dihedral(term.atoms)].append(term)
    summary.dihedral_duplicate_ff_keys = [
        " ".join(map(str, key)) for key, vals in dih_map.items() if len(vals) > 1
    ]
    itp_dih_keys = {d["key"] for d in proper_itp_dihedrals}
    for dih in dihedrals:
        if dih["funct"] == 2:
            summary.dihedral_skipped_improper.append(
                f"line {dih['line_no']}: {' '.join(map(str, dih['atoms']))}"
            )
            continue
        terms = dih_map.get(dih["key"], [])
        if not terms:
            summary.dihedral_itp_unmatched.append(
                f"line {dih['line_no']}: {' '.join(map(str, dih['atoms']))}"
            )
            continue
        if dih["funct"] != 1:
            summary.dihedral_skipped_wrong_funct.append(
                f"line {dih['line_no']}: funct={dih['funct']}"
            )
            continue
        term = terms[0]
        replace_line(out, dih["line_no"], format_stage_a_dihedral_line(dih, term))
        summary.proper_dihedral_lines_updated += 1
        changes.append(
            {
                "stage": "A_dihedral",
                "line_no": dih["line_no"],
                "atoms": " ".join(map(str, dih["atoms"])),
                "old_phase_deg": dih["phase_deg"],
                "new_phase_deg": term.phase_deg,
                "old_cp_kj": dih["cp_kj"],
                "new_cp_kj": term.cp_kj,
                "old_mult": dih["mult"],
                "new_mult": term.multiplicity,
            }
        )
    for key, terms in dih_map.items():
        if key not in itp_dih_keys:
            summary.dihedral_ff_unmatched.append(
                f"{' '.join(map(str, key))} ({'-'.join(terms[0].labels)})"
            )

    summary.ff_parse_warnings = list(ff_warnings)
    return out, summary, changes


def stage_a_is_complete(summary: StageASummary) -> bool:
    return (
        summary.charges_updated == summary.atoms_in_itp
        and not summary.charge_missing_atom_indices
        and not summary.charge_duplicate_atom_indices
        and not summary.charge_invalid_labels
        and summary.bond_lines_updated == summary.bond_lines_in_itp
        and not summary.bond_itp_unmatched
        and not summary.bond_skipped_wrong_funct
        and summary.angle_lines_updated == summary.angle_lines_in_itp
        and not summary.angle_itp_unmatched
        and not summary.angle_skipped_wrong_funct
        and summary.proper_dihedral_lines_updated == summary.proper_dihedral_lines_in_itp
        and not summary.dihedral_itp_unmatched
        and not summary.dihedral_skipped_wrong_funct
    )


# =============================================================================
# Stage B scan mapping, CSV reading, angle-aware pairing, fitting
# =============================================================================

def parse_scan_candidates(dihe_json_path: Path, stage_a_lines: List[str], match_reversed: bool) -> List[ScanCandidate]:
    data = json.loads(dihe_json_path.read_text(errors="replace"))
    torsions = data.get("scan_torsions")
    if not isinstance(torsions, list):
        raise ValueError(f"{dihe_json_path} does not contain list-valued scan_torsions")

    dihedrals = parse_dihedrals(stage_a_lines)
    out: List[ScanCandidate] = []

    for idx, item in enumerate(torsions):
        if not isinstance(item, list) or len(item) != 4:
            raise ValueError(f"Bad scan_torsions[{idx}] entry: {item!r}")
        torsion = tuple(int(x) for x in item)
        possible = {torsion}
        if match_reversed:
            possible.add(tuple(reversed(torsion)))

        matches = [d for d in dihedrals if d["funct"] != 2 and d["atoms"] in possible]
        issues: List[str] = []
        phase, cp, mult = 180.0, 0.0, 2
        lines = tuple(d["line_no"] for d in matches)
        if len(matches) == 0:
            issues.append("scan torsion not found uniquely as a proper [ dihedrals ] line after Stage A")
        elif len(matches) > 1:
            issues.append(f"scan torsion maps to multiple proper dihedrals after Stage A: lines {list(lines)}")
            d0 = matches[0]
            phase, cp, mult = d0["phase_deg"], d0["cp_kj"], d0["mult"]
        else:
            d0 = matches[0]
            phase, cp, mult = d0["phase_deg"], d0["cp_kj"], d0["mult"]

        out.append(
            ScanCandidate(
                scan_idx=idx,
                torsion=torsion,  # type: ignore[arg-type]
                central_bond=(torsion[1], torsion[2]),
                stage_a_phase_deg=phase,
                stage_a_cp_kj=cp,
                stage_a_mult=mult,
                itp_match_count=len(matches),
                itp_match_lines=lines,
                mapping_issues=issues,
            )
        )
    return out


def read_numeric_csv(path: Path) -> List[List[float]]:
    rows: List[List[float]] = []
    for raw in read_lines(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        toks = [t for t in re.split(r"[,\s]+", line) if t]
        try:
            vals = [safe_float(t) for t in toks]
        except ValueError:
            continue
        if vals:
            rows.append(vals)
    return rows


def read_qm_scan(path: Path, unit: str) -> Tuple[np.ndarray, np.ndarray]:
    rows = read_numeric_csv(path)
    if not rows or len(rows[0]) < 2:
        raise ValueError(f"{path} lacks numeric angle,energy rows")
    angles = np.array([row[0] for row in rows], dtype=float)
    energy = np.array([row[1] for row in rows], dtype=float)

    if unit == "kcal":
        rel = energy
    elif unit == "kj":
        rel = energy / KCAL_TO_KJ
    elif unit == "hartree":
        rel = (energy - np.min(energy)) * HARTREE_TO_KCAL
    else:
        raise ValueError(f"Unsupported QM energy unit: {unit}")
    return angles, rel


def read_mm_scan(path: Path, unit: str) -> Tuple[np.ndarray, np.ndarray]:
    rows = read_numeric_csv(path)
    if not rows:
        raise ValueError(f"{path} lacks numeric rows")
    if len(rows[0]) >= 2:
        step_raw = np.array([row[0] for row in rows], dtype=float)
        energy = np.array([row[1] for row in rows], dtype=float)
    else:
        step_raw = np.arange(len(rows), dtype=float)
        energy = np.array([row[0] for row in rows], dtype=float)

    if unit == "auto":
        unit = "hartree" if np.nanmax(np.abs(energy)) > 50 else "kcal"

    if unit == "hartree":
        rel = (energy - np.min(energy)) * HARTREE_TO_KCAL
    elif unit == "kcal":
        rel = energy - np.min(energy)
    elif unit == "kj":
        rel = (energy - np.min(energy)) / KCAL_TO_KJ
    else:
        raise ValueError(f"Unsupported MM energy unit: {unit}")
    return step_raw, rel


def mm_steps_to_indices(step_raw: np.ndarray, full_points: int) -> Tuple[List[int], List[str]]:
    warnings: List[str] = []
    rounded = [int(round(x)) for x in step_raw.tolist()]
    if all(abs(x - r) < 1e-6 for x, r in zip(step_raw.tolist(), rounded)):
        unique = len(set(rounded)) == len(rounded)
        if unique and all(0 <= r <= full_points - 1 for r in rounded):
            return rounded, warnings

    warnings.append(
        "MM first column is not a unique 0..N-1 scan-step index; using row order as scan-step index."
    )
    return list(range(len(step_raw))), warnings


def expected_angles_from_start(start_angle: float, full_points: int, step_deg: float) -> List[float]:
    return [circular_angle_deg(start_angle + i * step_deg) for i in range(full_points)]


def best_monotonic_step_assignment(
    observed_angles: Sequence[float],
    expected_angles: Sequence[float],
) -> Tuple[List[int], List[float]]:
    """
    Assign observed QM scan angles to an increasing subset of expected scan-step
    indices by dynamic programming, minimizing total circular angle deviation.

    This correctly handles missing points in the middle of a scan and the
    closure point whose angle may wrap back to the starting angle.
    """
    m = len(observed_angles)
    n = len(expected_angles)
    if m == 0:
        return [], []
    if m > n:
        raise ValueError(f"More observed QM points ({m}) than expected scan grid points ({n})")

    inf = float("inf")
    dp = np.full((m, n), inf)
    prev = np.full((m, n), -1, dtype=int)

    for j in range(n):
        dp[0, j] = circular_distance_deg(observed_angles[0], expected_angles[j])

    for i in range(1, m):
        best_cost = inf
        best_j = -1
        for j in range(n):
            # Update prefix best using preceding expected index.
            candidate_prev = j - 1
            if candidate_prev >= 0 and dp[i - 1, candidate_prev] < best_cost:
                best_cost = dp[i - 1, candidate_prev]
                best_j = candidate_prev
            if best_j >= 0:
                dp[i, j] = best_cost + circular_distance_deg(observed_angles[i], expected_angles[j])
                prev[i, j] = best_j

    end_j = int(np.argmin(dp[m - 1, :]))
    if not math.isfinite(float(dp[m - 1, end_j])):
        raise ValueError("Could not construct monotonic scan-step assignment")

    steps = [0] * m
    steps[-1] = end_j
    for i in range(m - 1, 0, -1):
        steps[i - 1] = int(prev[i, steps[i]])
        if steps[i - 1] < 0:
            raise ValueError("Broken dynamic-programming predecessor in scan-step assignment")

    errors = [
        circular_distance_deg(observed_angles[i], expected_angles[step])
        for i, step in enumerate(steps)
    ]
    return steps, errors


def pair_scan_points(
    qm_angles: np.ndarray,
    qm_rel_kcal: np.ndarray,
    mm_step_raw: np.ndarray,
    mm_rel_kcal: np.ndarray,
    full_points: int,
    step_deg: float,
) -> Tuple[List[PointPair], List[str]]:
    warnings: List[str] = []
    if len(qm_angles) == 0:
        return [], ["No QM points"]
    expected = expected_angles_from_start(float(qm_angles[0]), full_points, step_deg)
    qm_steps, angle_errors = best_monotonic_step_assignment(qm_angles.tolist(), expected)
    mm_steps, mm_warnings = mm_steps_to_indices(mm_step_raw, full_points)
    warnings.extend(mm_warnings)

    mm_map: Dict[int, Tuple[int, float]] = {}
    for row_idx, (step, energy) in enumerate(zip(mm_steps, mm_rel_kcal.tolist())):
        if step in mm_map:
            warnings.append(f"Duplicate MM scan-step index {step}; keeping first occurrence.")
            continue
        mm_map[step] = (row_idx, float(energy))

    pairs: List[PointPair] = []
    missing_mm_steps: List[int] = []
    for q_row, (step, angle, q_e, err) in enumerate(
        zip(qm_steps, qm_angles.tolist(), qm_rel_kcal.tolist(), angle_errors)
    ):
        if step not in mm_map:
            missing_mm_steps.append(step)
            continue
        mm_row, mm_e = mm_map[step]
        pairs.append(
            PointPair(
                qm_row=q_row,
                scan_step=step,
                qm_angle_deg=float(angle),
                expected_angle_deg=float(expected[step]),
                angle_error_deg=float(err),
                qm_rel_kcal=float(q_e),
                mm_rel_kcal=float(mm_e),
                target_kcal=float(q_e - mm_e),
                mm_source_row=int(mm_row),
            )
        )
    if missing_mm_steps:
        warnings.append(f"MM energies missing inferred scan-step indices {sorted(set(missing_mm_steps))}")

    return pairs, warnings


def gaussian_log_summary(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {"exists": False}
    text = path.read_text(errors="replace")
    return {
        "exists": True,
        "normal_termination": "Normal termination" in text,
        "error_termination": "Error termination" in text,
        "optimization_completed_count": text.count("Optimization completed"),
        "optimization_stopped_count": text.count("Optimization stopped"),
        "scf_failure": ("Convergence failure" in text or "Convergence criterion not met" in text),
        "mm_function_not_complete": "MM function not complete" in text,
        "undefined_terms": len(re.findall(r"undefined", text, flags=re.I)),
    }


def scan_log_diagnostics(data_dir: Path, scan_idx: int) -> Dict[str, Any]:
    qm_logs = sorted(data_dir.glob(f"{scan_idx}_qm*.log"))
    mm_logs = sorted(data_dir.glob(f"{scan_idx}_mm_*.log"))
    qm = [gaussian_log_summary(p) | {"file": p.name} for p in qm_logs]
    mm = [gaussian_log_summary(p) | {"file": p.name} for p in mm_logs]
    return {
        "qm_logs_found": len(qm_logs),
        "mm_logs_found": len(mm_logs),
        "qm_logs": qm,
        "mm_logs": mm,
        "failed_mm_logs": [
            x for x in mm if x.get("error_termination") or not x.get("normal_termination", False)
        ],
    }


def fit_periodic_one_term(angle_deg: np.ndarray, target_kcal: np.ndarray, mult: int) -> Tuple[float, float, float]:
    phi = np.deg2rad(angle_deg)
    x = np.column_stack([np.ones(len(phi)), np.cos(mult * phi), np.sin(mult * phi)])
    coeff, *_ = np.linalg.lstsq(x, target_kcal, rcond=None)
    pred = x @ coeff
    a = float(coeff[1])
    b = float(coeff[2])
    cp_kcal = math.hypot(a, b)
    phase_deg = math.degrees(math.atan2(b, a)) % 360.0
    rmse_kcal = float(np.sqrt(np.mean((target_kcal - pred) ** 2)))
    return phase_deg, cp_kcal, rmse_kcal


def fit_scan(
    data_dir: Path,
    candidate: ScanCandidate,
    args: argparse.Namespace,
    output_dir: Path,
) -> ScanFit:
    qm_file = choose_scan_energy_file(data_dir, candidate.scan_idx, "qm", args.full_scan_points)
    mm_file = choose_scan_energy_file(data_dir, candidate.scan_idx, "mm", args.full_scan_points)
    reasons = list(candidate.mapping_issues)
    pairing_warnings: List[str] = []
    logs = scan_log_diagnostics(data_dir, candidate.scan_idx)

    if qm_file is None:
        reasons.append("Missing QM scan-energy CSV")
    if mm_file is None:
        reasons.append("Missing MM scan-energy CSV")
    if qm_file is None or mm_file is None:
        return ScanFit(
            candidate.scan_idx, candidate.torsion, candidate.central_bond,
            qm_file.name if qm_file else None,
            mm_file.name if mm_file else None,
            0, 0, 0, False, reasons, candidate.stage_a_mult,
            None, None, None, None, None,
            candidate.stage_a_phase_deg, candidate.stage_a_cp_kj, candidate.stage_a_mult,
            candidate.itp_match_count, candidate.itp_match_lines,
            None, pairing_warnings, logs,
        )

    try:
        qm_angles, qm_rel = read_qm_scan(qm_file, args.qm_unit)
        mm_raw_steps, mm_rel = read_mm_scan(mm_file, args.mm_unit)
        pairs, pairing_warnings = pair_scan_points(
            qm_angles, qm_rel, mm_raw_steps, mm_rel,
            args.full_scan_points, args.scan_step_deg
        )
    except Exception as exc:
        reasons.append(f"Failed to read or pair scan-energy files: {exc}")
        return ScanFit(
            candidate.scan_idx, candidate.torsion, candidate.central_bond,
            qm_file.name, mm_file.name,
            0, 0, 0, False, reasons, candidate.stage_a_mult,
            None, None, None, None, None,
            candidate.stage_a_phase_deg, candidate.stage_a_cp_kj, candidate.stage_a_mult,
            candidate.itp_match_count, candidate.itp_match_lines,
            None, pairing_warnings, logs,
        )

    qm_count = len(qm_angles)
    mm_count = len(mm_rel)
    paired_count = len(pairs)
    max_angle_error = max((p.angle_error_deg for p in pairs), default=None)

    if qm_count != mm_count:
        reasons.append(f"QM/MM raw point-count mismatch: QM={qm_count}, MM={mm_count}")
    if paired_count < args.expected_points:
        reasons.append(f"Too few paired scan points: {paired_count} < required {args.expected_points}")
    if max_angle_error is not None and max_angle_error > args.angle_match_tolerance_deg:
        reasons.append(
            f"Angle-to-grid assignment error too large: max={max_angle_error:.3f} deg "
            f"> tolerance {args.angle_match_tolerance_deg:.3f} deg"
        )
    if pairing_warnings:
        reasons.extend(f"Pairing warning: {w}" for w in pairing_warnings)

    pair_rows = [asdict(p) for p in pairs]
    write_csv(output_dir / f"scan_{candidate.scan_idx}_paired_points.csv", pair_rows)

    if paired_count < 4:
        reasons.append("Fewer than 4 paired points; cannot fit a three-parameter periodic term.")
        return ScanFit(
            candidate.scan_idx, candidate.torsion, candidate.central_bond,
            qm_file.name, mm_file.name,
            qm_count, mm_count, paired_count, False, reasons, candidate.stage_a_mult,
            None, None, None, None, None,
            candidate.stage_a_phase_deg, candidate.stage_a_cp_kj, candidate.stage_a_mult,
            candidate.itp_match_count, candidate.itp_match_lines,
            max_angle_error, pairing_warnings, logs,
        )

    angle = np.array([p.qm_angle_deg for p in pairs], dtype=float)
    target = np.array([p.target_kcal for p in pairs], dtype=float)

    if args.fit_mult == "original":
        mult = candidate.stage_a_mult
        phase, cp_kcal, rmse_kcal = fit_periodic_one_term(angle, target, mult)
    elif args.fit_mult == "best":
        best = None
        for mult_trial in range(1, args.max_mult + 1):
            phase_t, cp_t, rmse_t = fit_periodic_one_term(angle, target, mult_trial)
            item = (rmse_t, mult_trial, phase_t, cp_t)
            if best is None or item[0] < best[0]:
                best = item
        assert best is not None
        rmse_kcal, mult, phase, cp_kcal = best
    else:
        mult = int(args.fit_mult)
        phase, cp_kcal, rmse_kcal = fit_periodic_one_term(angle, target, mult)

    cp_kj = cp_kcal * KCAL_TO_KJ
    rmse_kj = rmse_kcal * KCAL_TO_KJ

    if rmse_kj > args.rmse_max_kj:
        reasons.append(f"Fit RMSE too high: {rmse_kj:.3f} kJ/mol > limit {args.rmse_max_kj:.3f}")
    if cp_kj > args.k_max_kj:
        reasons.append(f"Fitted torsion cp too large: {cp_kj:.3f} kJ/mol > limit {args.k_max_kj:.3f}")

    trusted = True
    if candidate.mapping_issues:
        trusted = False
    if not args.allow_incomplete and qm_count != args.full_scan_points:
        trusted = False
        reasons.append(
            f"QM point count is {qm_count}; strict mode expects full scan count {args.full_scan_points}"
        )
    if not args.allow_incomplete and mm_count != args.full_scan_points:
        trusted = False
        reasons.append(
            f"MM point count is {mm_count}; strict mode expects full scan count {args.full_scan_points}"
        )
    if paired_count < args.expected_points:
        trusted = False
    if max_angle_error is not None and max_angle_error > args.angle_match_tolerance_deg:
        trusted = False
    if not args.allow_high_rmse and rmse_kj > args.rmse_max_kj:
        trusted = False
    if cp_kj > args.k_max_kj:
        trusted = False

    return ScanFit(
        candidate.scan_idx, candidate.torsion, candidate.central_bond,
        qm_file.name, mm_file.name,
        qm_count, mm_count, paired_count, trusted, reasons, mult,
        phase, cp_kj, cp_kcal, rmse_kj, rmse_kcal,
        candidate.stage_a_phase_deg, candidate.stage_a_cp_kj, candidate.stage_a_mult,
        candidate.itp_match_count, candidate.itp_match_lines,
        max_angle_error, pairing_warnings, logs,
    )


def format_scan_refined_dihedral(entry: Dict[str, Any], fit: ScanFit) -> str:
    i, j, k, l = entry["atoms"]
    comment = existing_comment(entry["line"])
    core = (
        f"{i:6d}{j:6d}{k:6d}{l:6d}{1:6d}"
        f"{fit.phase_deg:12.4f}{fit.cp_kj:12.5f}{fit.multiplicity:6d}"
    )
    if comment:
        core += "   " + comment
    core += (
        f"   ; QM/MM scan-refined idx={fit.scan_idx} "
        f"central={fit.central_bond[0]}-{fit.central_bond[1]} "
        f"RMSE={fit.rmse_kj:.2f} kJ/mol"
    )
    return core


def apply_stage_b(
    stage_a_lines: List[str],
    fits: List[ScanFit],
    match_reversed: bool,
) -> Tuple[List[str], List[Dict[str, Any]]]:
    out = list(stage_a_lines)
    dihedrals = parse_dihedrals(out)

    fit_map: Dict[Tuple[int, int, int, int], ScanFit] = {}
    for fit in fits:
        if not fit.trusted or fit.phase_deg is None or fit.cp_kj is None:
            continue
        fit_map[fit.torsion] = fit
        if match_reversed:
            fit_map[tuple(reversed(fit.torsion))] = fit

    changes: List[Dict[str, Any]] = []
    for entry in dihedrals:
        if entry["funct"] == 2:
            continue
        fit = fit_map.get(entry["atoms"])
        if fit is None:
            continue
        replace_line(out, entry["line_no"], format_scan_refined_dihedral(entry, fit))
        changes.append(
            {
                "stage": "B_scan_dihedral",
                "scan_idx": fit.scan_idx,
                "line_no": entry["line_no"],
                "atoms": " ".join(map(str, entry["atoms"])),
                "old_phase_deg_after_stage_a": entry["phase_deg"],
                "new_phase_deg": fit.phase_deg,
                "old_cp_kj_after_stage_a": entry["cp_kj"],
                "new_cp_kj": fit.cp_kj,
                "old_mult_after_stage_a": entry["mult"],
                "new_mult": fit.multiplicity,
                "rmse_kj": fit.rmse_kj,
                "paired_points": fit.paired_points,
            }
        )
    return out, changes


# =============================================================================
# Optional topology / XYZ diagnostics
# =============================================================================

def parse_topol_bonds(path: Path) -> set:
    lines = [line.strip() for line in read_lines(path) if line.strip()]
    if not lines:
        return set()
    try:
        n = int(lines[0].split()[0])
    except Exception:
        return set()
    bonds = set()
    for line in lines[1:1+n]:
        toks = line.split()
        if len(toks) >= 2:
            try:
                bonds.add(key_bond((int(toks[0]), int(toks[1]))))
            except Exception:
                pass
    return bonds


def infer_element(atom: Dict[str, Any]) -> str:
    mass = atom.get("mass")
    if mass is not None:
        if abs(mass - 1.008) < 0.25:
            return "H"
        if abs(mass - 12.011) < 0.8:
            return "C"
        if abs(mass - 14.007) < 0.8:
            return "N"
        if abs(mass - 15.999) < 0.8:
            return "O"
        if abs(mass - 32.06) < 1.2:
            return "S"
    name = str(atom.get("atom") or atom.get("type") or "").lstrip("0123456789")
    return name[:1].upper() if name else "X"


def parse_xyz_elements(path: Path) -> List[str]:
    lines = read_lines(path)
    if not lines:
        return []
    try:
        n = int(lines[0].split()[0])
    except Exception:
        return []
    out = []
    for line in lines[2:2+n]:
        toks = line.split()
        if toks:
            out.append(toks[0].capitalize())
    return out


def global_diagnostics(data_dir: Path, final_raw_lines: List[str], candidates: List[ScanCandidate]) -> Dict[str, Any]:
    atoms = parse_atoms(final_raw_lines)
    bonds = parse_bonds(final_raw_lines)
    diag: Dict[str, Any] = {
        "final_itp_atom_count": len(atoms),
        "final_itp_bond_count": len(bonds),
        "topol_check": None,
        "xyz_check": None,
        "scan_central_bond_check": [],
    }

    topol = choose_first_existing(data_dir.glob("topol*.txt"))
    if topol:
        topol_bonds = parse_topol_bonds(topol)
        itp_bonds = {b["key"] for b in bonds}
        diag["topol_check"] = {
            "file": topol.name,
            "topol_bond_count": len(topol_bonds),
            "itp_bond_count": len(itp_bonds),
            "topol_missing_in_itp_count": len(topol_bonds - itp_bonds),
            "itp_extra_vs_topol_count": len(itp_bonds - topol_bonds),
            "topol_missing_in_itp_first20": sorted(topol_bonds - itp_bonds)[:20],
            "itp_extra_vs_topol_first20": sorted(itp_bonds - topol_bonds)[:20],
        }
        for cand in candidates:
            k = key_bond(cand.central_bond)
            diag["scan_central_bond_check"].append(
                {
                    "scan_idx": cand.scan_idx,
                    "central_bond": list(cand.central_bond),
                    "in_topol": k in topol_bonds,
                    "in_itp": k in itp_bonds,
                }
            )

    xyz = choose_first_existing(data_dir.glob("*.xyz"))
    if xyz and atoms:
        xyz_elems = parse_xyz_elements(xyz)
        itp_elems = [infer_element(a).capitalize() for a in atoms]
        mismatches = []
        for idx, (itp_el, xyz_el) in enumerate(zip(itp_elems, xyz_elems), start=1):
            if itp_el != xyz_el:
                mismatches.append({"index": idx, "itp": itp_el, "xyz": xyz_el})
        diag["xyz_check"] = {
            "file": xyz.name,
            "xyz_atom_count": len(xyz_elems),
            "itp_atom_count": len(itp_elems),
            "atom_count_match": len(xyz_elems) == len(itp_elems),
            "element_order_mismatch_count": len(mismatches),
            "element_order_mismatches_first20": mismatches[:20],
        }
    return diag


# =============================================================================
# Output
# =============================================================================

def prepend_header(
    raw_lines: List[str],
    base_itp: Path,
    charge_file: Path,
    ff_file: Path,
    dihe_json: Path,
) -> List[str]:
    header = [
        "; -----------------------------------------------------------------------------",
        "; Generated by update_itp_trustable_dihe_json.py",
        f"; Original .itp: {base_itp.name}",
        f"; Charge source: {charge_file.name}",
        f"; HessFit force-string source: {ff_file.name}",
        f"; Scan torsion source: {dihe_json.name}",
        ";",
        "; Stage A:",
        ";   Atom charges updated from charge/type_charge text file.",
        ";   Bonds updated from HessFit HrmStr1.",
        ";   Angles updated from HessFit HrmBnd1.",
        ";   Proper dihedrals updated from HessFit AmbTrs.",
        ";   Original nonbonded atom types, [ pairs ], and improper dihedrals kept unchanged.",
        ";",
        "; Stage A conversion assumptions:",
        ";   Bonds: k = K*4.184*100/r0_nm^2, r0_nm = r0_A/10.",
        ";   Angles: k = 2*K*4.184/sin(theta0)^2.",
        ";   Proper dihedrals: cp_kJ = AmbTrs barrier_kcal*4.184.",
        ";",
        "; Stage B:",
        ";   Selected proper dihedrals optionally replaced by trusted QM/MM scan fits",
        ";   defined by dihe_optfile.json scan_torsions.",
        "; -----------------------------------------------------------------------------",
        "",
    ]
    return header + raw_lines


def write_stage_a_summary(path: Path, summary: StageASummary) -> None:
    rows = []
    for key, value in asdict(summary).items():
        if isinstance(value, list):
            rows.append({"metric": key, "value": " | ".join(map(str, value))})
        else:
            rows.append({"metric": key, "value": value})
    write_csv(path, rows)


def write_scan_mapping(path: Path, candidates: List[ScanCandidate]) -> None:
    rows = []
    for c in candidates:
        rows.append(
            {
                "scan_idx": c.scan_idx,
                "torsion": " ".join(map(str, c.torsion)),
                "central_bond": " ".join(map(str, c.central_bond)),
                "stage_a_phase_deg": c.stage_a_phase_deg,
                "stage_a_cp_kj": c.stage_a_cp_kj,
                "stage_a_mult": c.stage_a_mult,
                "itp_match_count": c.itp_match_count,
                "itp_match_lines": " ".join(map(str, c.itp_match_lines)),
                "mapping_issues": " | ".join(c.mapping_issues),
            }
        )
    write_csv(path, rows)


def write_scan_fit_results(path: Path, fits: List[ScanFit]) -> None:
    rows = []
    for f in fits:
        row = asdict(f)
        row["torsion"] = " ".join(map(str, f.torsion))
        row["central_bond"] = " ".join(map(str, f.central_bond))
        row["itp_match_lines"] = " ".join(map(str, f.itp_match_lines))
        row["skip_reasons"] = " | ".join(f.skip_reasons)
        row["pairing_warnings"] = " | ".join(f.pairing_warnings)
        row.pop("log_diagnostics", None)
        rows.append(row)
    write_csv(path, rows)


def write_report(
    path: Path,
    args: argparse.Namespace,
    base_itp: Path,
    charge_file: Path,
    ff_file: Path,
    dihe_json: Path,
    summary: StageASummary,
    stage_a_complete: bool,
    candidates: List[ScanCandidate],
    fits: List[ScanFit],
    stage_a_changes: List[Dict[str, Any]],
    stage_b_changes: List[Dict[str, Any]],
    global_diag: Dict[str, Any],
) -> None:
    trusted = [f.scan_idx for f in fits if f.trusted]
    skipped = [f.scan_idx for f in fits if not f.trusted]
    lines: List[str] = []
    lines.append("Full HessFit + QM/MM scan .itp update report")
    lines.append("=" * 76)
    lines.append("")
    lines.append("Inputs")
    lines.append(f"  data directory: {args.data}")
    lines.append(f"  original .itp: {base_itp.name}")
    lines.append(f"  charge file: {charge_file.name}")
    lines.append(f"  ff_string file: {ff_file.name}")
    lines.append(f"  dihe JSON: {dihe_json.name}")
    lines.append("")
    lines.append("Stage A — HessFit base update")
    lines.append(f"  complete by strict Stage-A criteria: {stage_a_complete}")
    lines.append(f"  charges updated: {summary.charges_updated}/{summary.atoms_in_itp}")
    lines.append(f"  bonds updated: {summary.bond_lines_updated}/{summary.bond_lines_in_itp}")
    lines.append(f"  angles updated: {summary.angle_lines_updated}/{summary.angle_lines_in_itp}")
    lines.append(
        f"  proper dihedrals updated from AmbTrs: "
        f"{summary.proper_dihedral_lines_updated}/{summary.proper_dihedral_lines_in_itp}"
    )
    lines.append(f"  missing charge atom indices: {summary.charge_missing_atom_indices}")
    lines.append(f"  duplicate charge atom indices: {summary.charge_duplicate_atom_indices}")
    lines.append(f"  invalid charge labels: {summary.charge_invalid_labels}")
    lines.append(f"  unmatched ITP bonds: {len(summary.bond_itp_unmatched)}")
    lines.append(f"  unmatched FF bond terms: {len(summary.bond_ff_unmatched)}")
    lines.append(f"  unmatched ITP angles: {len(summary.angle_itp_unmatched)}")
    lines.append(f"  unmatched FF angle terms: {len(summary.angle_ff_unmatched)}")
    lines.append(f"  unmatched ITP proper dihedrals: {len(summary.dihedral_itp_unmatched)}")
    lines.append(f"  unmatched FF dihedral terms: {len(summary.dihedral_ff_unmatched)}")
    lines.append(f"  ff_string parse warnings: {len(summary.ff_parse_warnings)}")
    lines.append("")
    lines.append("Stage B — QM/MM scan refinement")
    lines.append(f"  full scan points expected for strict mode: {args.full_scan_points}")
    lines.append(f"  minimum paired points required: {args.expected_points}")
    lines.append(f"  scan-step grid: {args.scan_step_deg:.6f} degree")
    lines.append(f"  trusted scans: {trusted}")
    lines.append(f"  skipped scans: {skipped}")
    lines.append(f"  Stage-B updated torsion lines: {len(stage_b_changes)}")
    lines.append("")
    for fit in fits:
        status = "TRUSTED/UPDATED" if fit.trusted else "SKIPPED"
        lines.append(
            f"  scan {fit.scan_idx}: D {' '.join(map(str, fit.torsion))}; "
            f"central {fit.central_bond[0]}-{fit.central_bond[1]} -> {status}"
        )
        lines.append(
            f"    files: QM={fit.qm_file}, MM={fit.mm_file}; "
            f"raw points QM/MM={fit.qm_points}/{fit.mm_points}; paired={fit.paired_points}"
        )
        lines.append(
            f"    .itp matches after Stage A: {fit.itp_match_count}; "
            f"lines={list(fit.itp_match_lines)}"
        )
        if fit.phase_deg is not None:
            lines.append(
                f"    fit: phase={fit.phase_deg:.4f} deg, cp={fit.cp_kj:.5f} kJ/mol, "
                f"mult={fit.multiplicity}, RMSE={fit.rmse_kj:.4f} kJ/mol"
            )
        if fit.max_angle_error_deg is not None:
            lines.append(f"    max angle-grid assignment error: {fit.max_angle_error_deg:.6f} deg")
        for reason in fit.skip_reasons:
            lines.append(f"    warning: {reason}")
        failed_mm = fit.log_diagnostics.get("failed_mm_logs", [])
        if failed_mm:
            lines.append(f"    failed/non-normal MM logs detected: {len(failed_mm)}")
    lines.append("")
    lines.append("Optional global checks")
    lines.append(f"  final .itp atoms: {global_diag.get('final_itp_atom_count')}")
    lines.append(f"  final .itp bonds: {global_diag.get('final_itp_bond_count')}")
    if global_diag.get("topol_check"):
        t = global_diag["topol_check"]
        lines.append(f"  topol file: {t['file']}")
        lines.append(f"  topol bonds missing in final .itp: {t['topol_missing_in_itp_count']}")
        lines.append(f"  final .itp extra bonds vs topol: {t['itp_extra_vs_topol_count']}")
    else:
        lines.append("  topol check: not run; no topol*.txt found")
    if global_diag.get("xyz_check"):
        x = global_diag["xyz_check"]
        lines.append(f"  xyz file: {x['file']}")
        lines.append(f"  xyz/.itp atom-count match: {x['atom_count_match']}")
        lines.append(f"  element-order mismatch count: {x['element_order_mismatch_count']}")
    else:
        lines.append("  xyz check: not run; no *.xyz found")
    lines.append("")
    lines.append("Output interpretation")
    lines.append(
        "  Final .itp = original .itp + Stage-A HessFit charges/bonds/angles/proper dihedrals "
        "+ trusted Stage-B scan refinements."
    )
    lines.append(
        "  Nonbonded atom types, [ pairs ], and improper dihedrals are not intentionally modified."
    )
    lines.append(
        "  In relaxed incomplete-scan mode, QM/MM points are paired by inferred scan-step index "
        "from the angle grid rather than blindly by CSV row."
    )
    path.write_text("\n".join(lines) + "\n")


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Full original-.itp -> HessFit base -> QM/MM scan-refined .itp updater."
    )
    parser.add_argument("--data", default="./data", help="Directory containing input data.")
    parser.add_argument("--base-itp", default=None, help="Original .itp inside --data or absolute path.")
    parser.add_argument("--dihe-json", default=None, help="dihe_optfile.json inside --data or absolute path.")
    parser.add_argument("--charges", default=None, help="type_charge.txt / charge.txt inside --data or absolute path.")
    parser.add_argument("--ff-string", default=None, help="ff_string.txt inside --data or absolute path.")
    parser.add_argument("--out-dir", default=None, help="Output directory. Default: ./data/itp_update_output")
    parser.add_argument("--output-itp", default="full_hessfit_scan_refined.itp")

    parser.add_argument("--qm-unit", default="kcal", choices=["kcal", "kj", "hartree"])
    parser.add_argument("--mm-unit", default="hartree", choices=["hartree", "kcal", "kj", "auto"])

    parser.add_argument("--full-scan-points", type=int, default=11,
                        help="Full expected scan point count, normally 11 for S 10 36.0.")
    parser.add_argument("--expected-points", type=int, default=11,
                        help="Minimum paired QM/MM points required to update a scan.")
    parser.add_argument("--scan-step-deg", type=float, default=36.0,
                        help="Expected torsion scan step in degrees.")
    parser.add_argument("--angle-match-tolerance-deg", type=float, default=2.0,
                        help="Maximum allowed QM angle deviation from inferred scan grid.")
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="Allow scans with fewer than full QM/MM points if paired-point threshold is met.")

    parser.add_argument("--fit-mult", default="original",
                        help="'original', 'best', or a specific integer multiplicity.")
    parser.add_argument("--max-mult", type=int, default=6,
                        help="Maximum multiplicity considered for --fit-mult best.")
    parser.add_argument("--rmse-max-kj", type=float, default=15.0)
    parser.add_argument("--k-max-kj", type=float, default=250.0)
    parser.add_argument("--allow-high-rmse", action="store_true")

    parser.add_argument("--require-complete-stage-a", action="store_true",
                        help="Abort unless all Stage-A charges/bonds/angles/proper dihedrals are updated cleanly.")
    parser.add_argument("--require-all-scans", action="store_true",
                        help="Abort unless every requested scan torsion is trusted.")
    parser.add_argument("--no-match-reversed", action="store_true",
                        help="Do not match reversed dihedral atom order in the .itp.")
    parser.add_argument("--no-backup", action="store_true")
    args = parser.parse_args()

    data_dir = Path(args.data).resolve()
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")
    out_dir = Path(args.out_dir).resolve() if args.out_dir else data_dir / "itp_update_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    base_itp = discover_base_itp(data_dir, args.base_itp)
    dihe_json = discover_dihe_json(data_dir, args.dihe_json)
    dihe_json_data = json.loads(dihe_json.read_text(errors="replace"))
    json_files = dihe_json_data.get("files", {}) if isinstance(dihe_json_data, dict) else {}

    charge_file = discover_aux_file_from_json_or_glob(
        data_dir, args.charges, "charge file", json_files,
        ["atom2type", "charge_file", "charges"],
        ["type_charge*.txt", "charge*.txt"],
    )
    ff_file = discover_aux_file_from_json_or_glob(
        data_dir, args.ff_string, "ff_string file", json_files,
        ["force_file", "ff_string", "forcefield"],
        ["ff_string*.txt"],
    )

    original_lines = read_lines(base_itp)
    charges, charge_invalid, charge_duplicates = parse_charges(charge_file)
    bond_terms, angle_terms, dihedral_terms, ff_warnings = parse_ff_string(ff_file)

    stage_a_lines, stage_a_summary, stage_a_changes = apply_stage_a(
        original_lines,
        charges,
        charge_invalid,
        charge_duplicates,
        bond_terms,
        angle_terms,
        dihedral_terms,
        ff_warnings,
    )
    stage_a_complete = stage_a_is_complete(stage_a_summary)

    if args.require_complete_stage_a and not stage_a_complete:
        # Write diagnostics before aborting.
        write_stage_a_summary(out_dir / "stage_A_hessfit_base_summary.csv", stage_a_summary)
        write_csv(out_dir / "stage_A_parameter_changes.csv", stage_a_changes)
        (out_dir / "manifest.json").write_text(json.dumps({
            "status": "aborted",
            "reason": "require-complete-stage-a requested, but Stage A diagnostics are not clean",
        }, indent=2) + "\n")
        raise SystemExit(
            f"Aborted: Stage A is incomplete. See {out_dir / 'stage_A_hessfit_base_summary.csv'}"
        )

    match_reversed = not args.no_match_reversed
    candidates = parse_scan_candidates(dihe_json, stage_a_lines, match_reversed=match_reversed)
    fits = [fit_scan(data_dir, candidate, args, out_dir) for candidate in candidates]

    if args.require_all_scans and any(not fit.trusted for fit in fits):
        global_diag = global_diagnostics(data_dir, stage_a_lines, candidates)
        write_stage_a_summary(out_dir / "stage_A_hessfit_base_summary.csv", stage_a_summary)
        write_csv(out_dir / "stage_A_parameter_changes.csv", stage_a_changes)
        write_scan_mapping(out_dir / "scan_definition_mapping.csv", candidates)
        write_scan_fit_results(out_dir / "scan_fit_results.csv", fits)
        write_report(
            out_dir / "diagnostics_report.txt",
            args, base_itp, charge_file, ff_file, dihe_json,
            stage_a_summary, stage_a_complete, candidates, fits,
            stage_a_changes, [], global_diag,
        )
        (out_dir / "manifest.json").write_text(json.dumps({
            "status": "aborted",
            "reason": "require-all-scans requested, but at least one scan is untrusted",
            "trusted_scans": [f.scan_idx for f in fits if f.trusted],
            "skipped_scans": [f.scan_idx for f in fits if not f.trusted],
        }, indent=2) + "\n")
        raise SystemExit(
            f"Aborted: untrusted scans remain {[f.scan_idx for f in fits if not f.trusted]}. "
            f"See {out_dir / 'diagnostics_report.txt'}"
        )

    final_raw_lines, stage_b_changes = apply_stage_b(stage_a_lines, fits, match_reversed=match_reversed)
    final_lines = prepend_header(final_raw_lines, base_itp, charge_file, ff_file, dihe_json)
    out_itp = out_dir / args.output_itp
    out_itp.write_text("\n".join(final_lines) + "\n")

    if not args.no_backup:
        shutil.copy2(base_itp, out_dir / (base_itp.name + ".bak"))

    global_diag = global_diagnostics(data_dir, final_raw_lines, candidates)

    write_stage_a_summary(out_dir / "stage_A_hessfit_base_summary.csv", stage_a_summary)
    write_csv(out_dir / "stage_A_parameter_changes.csv", stage_a_changes)
    write_scan_mapping(out_dir / "scan_definition_mapping.csv", candidates)
    write_scan_fit_results(out_dir / "scan_fit_results.csv", fits)
    write_csv(out_dir / "stage_B_scan_parameter_changes.csv", stage_b_changes)
    write_report(
        out_dir / "diagnostics_report.txt",
        args, base_itp, charge_file, ff_file, dihe_json,
        stage_a_summary, stage_a_complete, candidates, fits,
        stage_a_changes, stage_b_changes, global_diag,
    )

    manifest = {
        "status": "completed",
        "output_itp": str(out_itp),
        "input": {
            "data_dir": str(data_dir),
            "original_itp": base_itp.name,
            "charge_file": charge_file.name,
            "ff_string_file": ff_file.name,
            "dihe_json": dihe_json.name,
        },
        "stage_a": {
            "complete": stage_a_complete,
            "charges_updated": stage_a_summary.charges_updated,
            "bonds_updated": stage_a_summary.bond_lines_updated,
            "angles_updated": stage_a_summary.angle_lines_updated,
            "proper_dihedrals_updated": stage_a_summary.proper_dihedral_lines_updated,
        },
        "stage_b": {
            "trusted_scans": [f.scan_idx for f in fits if f.trusted],
            "skipped_scans": [f.scan_idx for f in fits if not f.trusted],
            "scan_refined_dihedrals_updated": len(stage_b_changes),
        },
        "outputs": {
            "final_itp": str(out_itp),
            "diagnostics_report": str(out_dir / "diagnostics_report.txt"),
            "stage_A_summary": str(out_dir / "stage_A_hessfit_base_summary.csv"),
            "stage_A_changes": str(out_dir / "stage_A_parameter_changes.csv"),
            "scan_mapping": str(out_dir / "scan_definition_mapping.csv"),
            "scan_fit_results": str(out_dir / "scan_fit_results.csv"),
            "stage_B_changes": str(out_dir / "stage_B_scan_parameter_changes.csv"),
        },
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print("Completed full original-.itp -> HessFit base -> QM/MM scan-refined .itp update.")
    print(f"  Output ITP: {out_itp}")
    print("")
    print("Stage A:")
    print(f"  Complete by strict Stage-A criteria: {stage_a_complete}")
    print(f"  Charges updated: {stage_a_summary.charges_updated}/{stage_a_summary.atoms_in_itp}")
    print(f"  Bonds updated: {stage_a_summary.bond_lines_updated}/{stage_a_summary.bond_lines_in_itp}")
    print(f"  Angles updated: {stage_a_summary.angle_lines_updated}/{stage_a_summary.angle_lines_in_itp}")
    print(
        "  Proper dihedrals updated from AmbTrs: "
        f"{stage_a_summary.proper_dihedral_lines_updated}/{stage_a_summary.proper_dihedral_lines_in_itp}"
    )
    print("")
    print("Stage B:")
    print(f"  Trusted scans: {[f.scan_idx for f in fits if f.trusted]}")
    print(f"  Skipped scans: {[f.scan_idx for f in fits if not f.trusted]}")
    print(f"  Scan-refined dihedrals updated: {len(stage_b_changes)}")
    print(f"  Diagnostics: {out_dir / 'diagnostics_report.txt'}")


if __name__ == "__main__":
    main()
