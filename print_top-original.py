import argparse
from pathlib import Path
from typing import List, Tuple


BOHR_TO_ANG = 0.529177210903


def infer_element(atom_name: str) -> str:
    """
    Infer element from a PDB atom name if columns 77-78 are missing.

    This is a fallback only. Proper PDB files should contain the element field.
    """
    name = atom_name.strip()

    # Remove leading digits, e.g. "1H" -> "H"
    name = name.lstrip("0123456789")

    if not name:
        return "X"

    # Common two-letter elements. Add more if needed.
    two_letter = {
        "Cl", "Br", "Si", "Na", "Li", "Mg", "Ca", "Fe", "Zn", "Cu", "Mn",
        "Co", "Ni", "Al", "Ag", "Au", "Pt", "Pd", "Sn", "Se", "As"
    }

    if len(name) >= 2:
        candidate = name[0].upper() + name[1].lower()
        if candidate in two_letter:
            return candidate

    return name[0].upper()


def read_pdb_atoms(pdb_path: Path, bohr_to_ang: bool = False) -> List[Tuple[str, float, float, float]]:
    atoms = []
    scale = BOHR_TO_ANG if bohr_to_ang else 1.0

    with pdb_path.open("r", errors="replace") as f:
        for line_number, line in enumerate(f, start=1):
            if not line.startswith(("ATOM", "HETATM")):
                continue

            try:
                x = float(line[30:38]) * scale
                y = float(line[38:46]) * scale
                z = float(line[46:54]) * scale
            except ValueError as exc:
                raise ValueError(
                    f"Could not parse coordinates on line {line_number}: {line.rstrip()}"
                ) from exc

            element = ""
            if len(line) >= 78:
                element = line[76:78].strip()

            if not element:
                atom_name = line[12:16]
                element = infer_element(atom_name)

            atoms.append((element, x, y, z))

    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {pdb_path}")

    return atoms


def write_xyz(atoms: List[Tuple[str, float, float, float]], xyz_path: Path, comment: str) -> None:
    with xyz_path.open("w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(comment + "\n")
        for element, x, y, z in atoms:
            f.write(f"{element:<2s} {x:14.6f} {y:14.6f} {z:14.6f}\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert PDB coordinates to XYZ format for HessFit."
    )
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument(
        "xyz_file",
        nargs="?",
        help="Output XYZ file. Default: input filename with .xyz extension",
    )
    parser.add_argument(
        "--bohr-to-ang",
        action="store_true",
        help="Multiply coordinates by 0.529177210903 before writing XYZ. Use only if input coordinates are in bohr.",
    )

    args = parser.parse_args()

    pdb_path = Path(args.pdb_file).expanduser().resolve()
    if args.xyz_file:
        xyz_path = Path(args.xyz_file).expanduser().resolve()
    else:
        xyz_path = pdb_path.with_suffix(".xyz")

    atoms = read_pdb_atoms(pdb_path, bohr_to_ang=args.bohr_to_ang)

    unit_note = "converted from bohr to Angstrom" if args.bohr_to_ang else "coordinates in Angstrom from PDB"
    comment = f"{pdb_path.name}; {unit_note}"

    write_xyz(atoms, xyz_path, comment)

    print(f"Wrote: {xyz_path}")
    print(f"Atoms: {len(atoms)}")
    print("Note: XYZ format contains coordinates only; bonds/topology must come from topol.txt or another topology file.")


if __name__ == "__main__":
    main()
