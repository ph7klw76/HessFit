from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
from collections import defaultdict, deque
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple


def strip_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def try_float(x: str) -> Optional[float]:
    try:
        return float(x)
    except Exception:
        return None


def try_int_float(x: str) -> Optional[int]:
    try:
        return int(float(x))
    except Exception:
        return None


def infer_element_from_mass_or_name(mass: Optional[float], atom_name: str, atom_type: str = "") -> str:
    if mass is not None:
        if abs(mass - 1.008) < 0.25:
            return "H"
        if abs(mass - 12.011) < 0.7:
            return "C"
        if abs(mass - 14.007) < 0.7:
            return "N"
        if abs(mass - 15.999) < 0.7:
            return "O"
        if abs(mass - 32.06) < 1.0:
            return "S"
        if abs(mass - 30.974) < 1.0:
            return "P"

    text = atom_name.strip() or atom_type.strip()
    text = text.lstrip("0123456789")
    if not text:
        return "X"

    two = text[:2].capitalize()
    if two in {"Cl", "Br", "Si", "Na", "Li", "Mg", "Ca", "Fe", "Zn", "Cu", "Mn", "Co", "Ni", "Al", "Se"}:
        return two

    return text[0].upper()


def parse_itp_sections(itp_path: Path) -> List[Tuple[str, List[str]]]:
    sections: List[Tuple[str, List[str]]] = []
    current: Optional[str] = None

    for raw in itp_path.read_text(errors="replace").splitlines():
        line = raw.strip()

        if line.startswith("[") and "]" in line:
            current = line.split("]", 1)[0].strip("[]").strip().lower()
            sections.append((current, []))
            continue

        if current:
            sections[-1][1].append(raw)

    return sections


def parse_atoms(sections: List[Tuple[str, List[str]]]) -> List[Dict]:
    atoms: List[Dict] = []

    for name, lines in sections:
        if name != "atoms":
            continue

        for raw in lines:
            s = strip_comment(raw)
            if not s:
                continue

            p = s.split()
            if len(p) < 7:
                continue

            try:
                nr = int(p[0])
                atom_type = p[1]
                atom_name = p[4] if len(p) >= 5 else f"A{nr}"
                charge = float(p[6]) if len(p) >= 7 else 0.0
                mass = float(p[7]) if len(p) >= 8 else None
            except ValueError:
                continue

            atoms.append({
                "nr": nr,
                "type": atom_type,
                "atom": atom_name,
                "charge": charge,
                "mass": mass,
                "element": infer_element_from_mass_or_name(mass, atom_name, atom_type),
            })

    atoms.sort(key=lambda x: x["nr"])
    return atoms


def parse_bonds(sections: List[Tuple[str, List[str]]]) -> List[Dict]:
    bonds: List[Dict] = []

    for name, lines in sections:
        if name != "bonds":
            continue

        for raw in lines:
            s = strip_comment(raw)
            if not s:
                continue

            p = s.split()
            if len(p) < 2:
                continue

            try:
                i = int(p[0])
                j = int(p[1])
                funct = int(float(p[2])) if len(p) >= 3 else None
                b0 = float(p[3]) if len(p) >= 4 else None
                kb = float(p[4]) if len(p) >= 5 else None
            except ValueError:
                continue

            bonds.append({"i": i, "j": j, "funct": funct, "b0": b0, "kb": kb, "raw": raw.rstrip()})

    return bonds


def parse_dihedrals(sections: List[Tuple[str, List[str]]], proper_functs: Set[int]) -> List[Dict]:
    dihedrals: List[Dict] = []
    section_counter = 0

    for name, lines in sections:
        if name != "dihedrals":
            continue

        section_counter += 1

        for raw in lines:
            s = strip_comment(raw)
            if not s:
                continue

            p = s.split()
            if len(p) < 5:
                continue

            try:
                ai, aj, ak, al = map(int, p[:4])
                funct = int(float(p[4]))
            except ValueError:
                continue

            if funct not in proper_functs:
                continue

            dihedrals.append({
                "ai": ai,
                "aj": aj,
                "ak": ak,
                "al": al,
                "funct": funct,
                "phi0": try_float(p[5]) if len(p) > 5 else None,
                "cp": try_float(p[6]) if len(p) > 6 else None,
                "mult": try_int_float(p[7]) if len(p) > 7 else None,
                "section_index": section_counter,
                "raw": raw.rstrip(),
            })

    return dihedrals


def parse_xyz(xyz_path: Path) -> List[Dict]:
    lines = xyz_path.read_text(errors="replace").splitlines()
    if not lines:
        raise ValueError(f"Empty XYZ file: {xyz_path}")

    nat = int(lines[0].split()[0])
    atoms: List[Dict] = []

    for idx, raw in enumerate(lines[2:2 + nat], start=1):
        p = raw.split()
        if len(p) < 4:
            raise ValueError(f"Bad XYZ atom line {idx + 2}: {raw}")

        atoms.append({
            "nr": idx,
            "element": p[0],
            "coord": (float(p[1]), float(p[2]), float(p[3])),
        })

    if len(atoms) != nat:
        raise ValueError(f"XYZ atom count mismatch: header says {nat}, parsed {len(atoms)}")

    return atoms


def build_graph(natoms: int, bonds: Sequence[Dict]) -> Dict[int, Set[int]]:
    graph: Dict[int, Set[int]] = {i: set() for i in range(1, natoms + 1)}
    for b in bonds:
        i, j = b["i"], b["j"]
        if i in graph and j in graph:
            graph[i].add(j)
            graph[j].add(i)
    return graph


def has_path_excluding_edge(graph: Dict[int, Set[int]], start: int, goal: int, edge: Tuple[int, int]) -> bool:
    a, b = edge
    q = deque([start])
    seen = {start}

    while q:
        u = q.popleft()
        if u == goal:
            return True

        for v in graph[u]:
            if (u == a and v == b) or (u == b and v == a):
                continue
            if v not in seen:
                seen.add(v)
                q.append(v)

    return False


def detect_ring_edges_and_atoms(graph: Dict[int, Set[int]], bonds: Sequence[Dict]) -> Tuple[Set[Tuple[int, int]], Set[int]]:
    ring_edges: Set[Tuple[int, int]] = set()

    for b in bonds:
        i, j = b["i"], b["j"]
        edge = tuple(sorted((i, j)))
        if has_path_excluding_edge(graph, i, j, edge):
            ring_edges.add(edge)

    ring_atoms = set()
    for i, j in ring_edges:
        ring_atoms.add(i)
        ring_atoms.add(j)

    return ring_edges, ring_atoms


def is_heavy(element: str) -> bool:
    return element.capitalize() != "H"


def dihedral_angle_deg(p0, p1, p2, p3) -> float:
    def sub(a, b):
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

    def dot(a, b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

    def cross(a, b):
        return (
            a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0],
        )

    def norm(a):
        return math.sqrt(dot(a, a))

    def scale(a, s):
        return (a[0]*s, a[1]*s, a[2]*s)

    def normalize(a):
        n = norm(a)
        if n == 0:
            return (0.0, 0.0, 0.0)
        return scale(a, 1.0 / n)

    b0 = sub(p0, p1)
    b1 = sub(p2, p1)
    b2 = sub(p3, p2)

    b1n = normalize(b1)
    v = sub(b0, scale(b1n, dot(b0, b1n)))
    w = sub(b2, scale(b1n, dot(b2, b1n)))

    x = dot(v, w)
    y = dot(cross(b1n, v), w)
    return math.degrees(math.atan2(y, x))


def score_dihedral(
    d: Dict,
    atoms: List[Dict],
    ring_atoms: Set[int],
    ring_edges: Set[Tuple[int, int]],
    bridge_central_edges: Set[Tuple[int, int]],
) -> int:
    ai, aj, ak, al = d["ai"], d["aj"], d["ak"], d["al"]
    central = tuple(sorted((aj, ak)))

    score = 0

    if central in bridge_central_edges:
        score += 1000

    for outer in (ai, al):
        if is_heavy(atoms[outer - 1]["element"]):
            score += 100
        else:
            score -= 200

        if outer in ring_atoms:
            score += 50

    if aj in ring_atoms:
        score += 30
    if ak in ring_atoms:
        score += 30

    if central in ring_edges:
        score -= 1000

    cp = d.get("cp")
    if cp is not None:
        score += min(int(abs(cp)), 100)

    return score


def recommend_torsions(
    atoms: List[Dict],
    xyz_atoms: List[Dict],
    bonds: List[Dict],
    dihedrals: List[Dict],
    max_per_bond: int = 1,
    include_nonring_rotors: bool = False,
) -> Tuple[List[Dict], List[Dict], Dict]:
    natoms = len(atoms)
    graph = build_graph(natoms, bonds)
    ring_edges, ring_atoms = detect_ring_edges_and_atoms(graph, bonds)
    bond_edges = {tuple(sorted((b["i"], b["j"]))) for b in bonds}

    bridge_edges: Set[Tuple[int, int]] = set()
    for b in bonds:
        edge = tuple(sorted((b["i"], b["j"])))
        i, j = edge

        if edge in ring_edges:
            continue

        if i in ring_atoms and j in ring_atoms:
            if is_heavy(atoms[i - 1]["element"]) and is_heavy(atoms[j - 1]["element"]):
                bridge_edges.add(edge)

        if include_nonring_rotors:
            if is_heavy(atoms[i - 1]["element"]) and is_heavy(atoms[j - 1]["element"]):
                bridge_edges.add(edge)

    candidates: List[Dict] = []
    grouped: Dict[Tuple[int, int], List[Dict]] = defaultdict(list)

    for d in dihedrals:
        ai, aj, ak, al = d["ai"], d["aj"], d["ak"], d["al"]
        if any(x < 1 or x > natoms for x in [ai, aj, ak, al]):
            continue

        central = tuple(sorted((aj, ak)))

        if central not in bond_edges:
            continue
        if central in ring_edges:
            continue
        if central not in bridge_edges:
            continue

        try:
            angle = dihedral_angle_deg(
                xyz_atoms[ai - 1]["coord"],
                xyz_atoms[aj - 1]["coord"],
                xyz_atoms[ak - 1]["coord"],
                xyz_atoms[al - 1]["coord"],
            )
        except Exception:
            angle = float("nan")

        item = {
            "scan_torsion": [ai, aj, ak, al],
            "central_bond": [aj, ak],
            "central_bond_sorted": list(central),
            "central_bond_atom_names": f"{atoms[aj - 1]['atom']}-{atoms[ak - 1]['atom']}",
            "central_bond_atom_types": f"{atoms[aj - 1]['type']}-{atoms[ak - 1]['type']}",
            "central_bond_elements": f"{atoms[aj - 1]['element']}-{atoms[ak - 1]['element']}",
            "outer_atom_names": f"{atoms[ai - 1]['atom']} ... {atoms[al - 1]['atom']}",
            "outer_atom_elements": f"{atoms[ai - 1]['element']} ... {atoms[al - 1]['element']}",
            "current_dihedral_deg": round(angle, 3),
            "itp_funct": d.get("funct"),
            "itp_phi0_deg": d.get("phi0"),
            "itp_cp": d.get("cp"),
            "itp_mult": d.get("mult"),
            "score": score_dihedral(d, atoms, ring_atoms, ring_edges, bridge_edges),
            "raw_itp_line": d.get("raw", ""),
        }

        candidates.append(item)
        grouped[central].append(item)

    recommended: List[Dict] = []
    for central, items in grouped.items():
        items_sorted = sorted(items, key=lambda x: x["score"], reverse=True)
        recommended.extend(items_sorted[:max_per_bond])

    recommended = sorted(recommended, key=lambda x: (x["central_bond_sorted"], -x["score"]))

    meta = {
        "ring_atom_count": len(ring_atoms),
        "ring_edge_count": len(ring_edges),
        "bridge_central_bond_count": len(bridge_edges),
        "proper_dihedral_count": len(dihedrals),
        "candidate_count": len(candidates),
        "recommended_count": len(recommended),
    }

    return recommended, candidates, meta


def check_atom_order(atoms: List[Dict], xyz_atoms: List[Dict]) -> List[Dict]:
    mismatches = []

    n = min(len(atoms), len(xyz_atoms))
    for idx in range(n):
        itp_e = atoms[idx]["element"].capitalize()
        xyz_e = xyz_atoms[idx]["element"].capitalize()

        if itp_e != xyz_e:
            mismatches.append({
                "index": idx + 1,
                "itp_element": itp_e,
                "xyz_element": xyz_e,
                "itp_atom_name": atoms[idx]["atom"],
                "itp_atom_type": atoms[idx]["type"],
            })

    return mismatches


def update_dihe_json(
    json_file: Path,
    scan_torsions: List[List[int]],
    yes: bool = False,
    dry_run: bool = False,
    backup: bool = True,
) -> None:
    if not json_file.exists():
        raise FileNotFoundError(f"JSON file not found: {json_file}")

    data = json.loads(json_file.read_text())
    old = data.get("scan_torsions", None)

    data["scan_torsions"] = scan_torsions

    print(f"\nUpdating JSON file: {json_file}")
    if old is None:
        print("  Existing scan_torsions: none")
    else:
        print(f"  Existing scan_torsions: {old}")
    print(f"  New scan_torsions:      {scan_torsions}")

    if dry_run:
        print("\nDry run only. No JSON file was changed.")
        print(json.dumps(data, indent=2))
        return

    if old is not None and old != scan_torsions and not yes:
        answer = input("\nOverwrite existing scan_torsions? Type YES to continue: ").strip()
        if answer != "YES":
            print("Aborted. JSON file was not changed.")
            return

    if backup:
        backup_path = json_file.with_suffix(json_file.suffix + ".bak")
        shutil.copy2(json_file, backup_path)
        print(f"  Backup written: {backup_path}")

    json_file.write_text(json.dumps(data, indent=2) + "\n")
    print("  JSON updated successfully.")


def write_outputs(prefix: str, atoms, xyz_atoms, recommended, candidates, mismatches, meta) -> None:
    json_path = Path(f"{prefix}_scan_torsions.json")
    csv_path = Path(f"{prefix}_torsion_refinement_candidates.csv")
    check_path = Path(f"{prefix}_atom_order_check.txt")

    json_data = {
        "scan_torsions": [x["scan_torsion"] for x in recommended],
        "notes": {
            "meaning": "Each list [i, j, k, l] means Gaussian ModRedundant D i j k l; central scanned bond is j-k.",
            "selection": "Recommended torsions are proper .itp dihedrals whose central bond is a non-ring heavy-heavy bridge between ring atoms.",
        },
        "meta": meta,
    }
    json_path.write_text(json.dumps(json_data, indent=2) + "\n")

    fieldnames = [
        "scan_torsion",
        "central_bond",
        "central_bond_atom_names",
        "central_bond_atom_types",
        "central_bond_elements",
        "outer_atom_names",
        "outer_atom_elements",
        "current_dihedral_deg",
        "itp_funct",
        "itp_phi0_deg",
        "itp_cp",
        "itp_mult",
        "score",
        "raw_itp_line",
    ]

    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for c in sorted(candidates, key=lambda x: (x["central_bond_sorted"], -x["score"])):
            row = dict(c)
            row["scan_torsion"] = " ".join(map(str, row["scan_torsion"]))
            row["central_bond"] = " ".join(map(str, row["central_bond"]))
            row.pop("central_bond_sorted", None)
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    lines = [
        f"ITP atoms: {len(atoms)}",
        f"XYZ atoms: {len(xyz_atoms)}",
        f"Element-order mismatches: {len(mismatches)}",
    ]

    if not mismatches:
        lines.append("Atom order matches by element sequence.")
    else:
        lines.append("Mismatches:")
        for m in mismatches:
            lines.append(
                f"  {m['index']}: ITP {m['itp_element']} ({m['itp_atom_name']}/{m['itp_atom_type']}) "
                f"vs XYZ {m['xyz_element']}"
            )

    lines.append("")
    lines.append("Meta:")
    for k, v in meta.items():
        lines.append(f"  {k}: {v}")

    check_path.write_text("\n".join(lines) + "\n")

    print(f"\nWrote {json_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {check_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recommend QM torsion-scan dihedrals from a GROMACS .itp and matching .xyz, and optionally update dihe_optfile.json."
    )
    parser.add_argument("itp", help="Input GROMACS .itp file")
    parser.add_argument("xyz", help="Input XYZ file with matching atom order")
    parser.add_argument("--prefix", default=None, help="Output prefix. Default: stem of the .itp file")
    parser.add_argument("--max-per-bond", type=int, default=1, help="Representative scan torsions per central bond. Default: 1")
    parser.add_argument("--proper-functs", default="1,3,4,5,9", help="Comma-separated proper dihedral function numbers. Default: 1,3,4,5,9")
    parser.add_argument("--include-nonring-rotors", action="store_true", help="Also recommend non-ring heavy-heavy rotatable central bonds.")
    parser.add_argument("--json-only", action="store_true", help="Print only the scan_torsions JSON block.")
    parser.add_argument("--update-json", default=None, help="Path to dihe_optfile.json to update automatically.")
    parser.add_argument("--yes", action="store_true", help="Overwrite existing scan_torsions without prompting.")
    parser.add_argument("--dry-run", action="store_true", help="Preview JSON update without writing.")
    parser.add_argument("--no-backup", action="store_true", help="Do not create .bak file before updating JSON.")

    args = parser.parse_args()

    itp_path = Path(args.itp)
    xyz_path = Path(args.xyz)
    prefix = args.prefix or itp_path.with_suffix("").name
    proper_functs = {int(x.strip()) for x in args.proper_functs.split(",") if x.strip()}

    sections = parse_itp_sections(itp_path)
    atoms = parse_atoms(sections)
    bonds = parse_bonds(sections)
    dihedrals = parse_dihedrals(sections, proper_functs)
    xyz_atoms = parse_xyz(xyz_path)

    if len(atoms) != len(xyz_atoms):
        raise RuntimeError(f"Atom count mismatch: ITP has {len(atoms)}, XYZ has {len(xyz_atoms)}")

    mismatches = check_atom_order(atoms, xyz_atoms)
    if mismatches:
        print("WARNING: Atom element order mismatches were found.")
        print("Do not trust scan_torsions until the atom-order report is inspected.")

    recommended, candidates, meta = recommend_torsions(
        atoms=atoms,
        xyz_atoms=xyz_atoms,
        bonds=bonds,
        dihedrals=dihedrals,
        max_per_bond=args.max_per_bond,
        include_nonring_rotors=args.include_nonring_rotors,
    )

    scan_torsions = [x["scan_torsion"] for x in recommended]
    json_block = {"scan_torsions": scan_torsions}

    if args.json_only:
        print(json.dumps(json_block, indent=2))
        return

    print("Atom order check:")
    print(f"  ITP atoms: {len(atoms)}")
    print(f"  XYZ atoms: {len(xyz_atoms)}")
    print(f"  Element-order mismatches: {len(mismatches)}")

    print("\nRecommended scan_torsions JSON block:")
    print(json.dumps(json_block, indent=2))

    print("\nRecommended torsion details:")
    for r in recommended:
        print(
            f"  D {' '.join(map(str, r['scan_torsion']))} "
            f"| central {r['central_bond_atom_names']} ({r['central_bond_atom_types']}) "
            f"| angle {r['current_dihedral_deg']} deg "
            f"| cp {r['itp_cp']} | score {r['score']}"
        )

    write_outputs(prefix, atoms, xyz_atoms, recommended, candidates, mismatches, meta)

    if args.update_json:
        update_dihe_json(
            json_file=Path(args.update_json),
            scan_torsions=scan_torsions,
            yes=args.yes,
            dry_run=args.dry_run,
            backup=not args.no_backup,
        )


if __name__ == "__main__":
    main()
