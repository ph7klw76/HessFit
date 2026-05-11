#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

HARTREE_TO_KCAL = 627.5094740631

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-t", "--type", dest="ftype", default="qm", help="qm or mm")
    p.add_argument("-f", "--file", dest="file_log", required=True, help="Gaussian log file")
    p.add_argument("-o", "--output", dest="output", required=True, help="Output CSV")
    return p.parse_args()

def read_text(path):
    return Path(path).read_text(errors="replace").splitlines()

def check_normal_termination(lines, filename):
    if not any("Normal termination" in line for line in lines):
        bad = []
        for key in ["Error termination", "Convergence failure", "Convergence criterion not met"]:
            if any(key in line for line in lines):
                bad.append(key)
        msg = ", ".join(bad) if bad else "no Normal termination found"
        raise RuntimeError(
            f"{filename} is not a completed Gaussian scan log: {msg}. "
            "Do not extract scan energies from failed logs."
        )

def parse_qm_scan(lines):
    """
    Robust Gaussian relaxed-scan parser.

    It detects scan-coordinate lines such as:
      ! D4    D(10,2,3,5)      0.0815      Scan !

    It pairs each scan angle with the final SCF energy before the next scan angle.
    Energies are returned relative to the minimum in kcal/mol.
    """
    scan_pat = re.compile(
        r"!\s*D\d+\s+D\([^)]*\)\s+"
        r"([-+]?(?:\d+\.\d*|\d+|\.\d+)(?:[DEde][-+]?\d+)?)\s+Scan"
    )
    energy_pat = re.compile(
        r"SCF Done:\s+E\([^)]+\)\s+=\s+([-+]?\d+\.\d+)"
    )

    records = []
    current_angle = None
    last_energy = None

    def d_to_float(x):
        return float(x.replace("D", "E").replace("d", "E"))

    for line in lines:
        em = energy_pat.search(line)
        if em:
            last_energy = float(em.group(1))

        sm = scan_pat.search(line)
        if sm:
            angle = d_to_float(sm.group(1))

            if current_angle is None:
                current_angle = angle
                continue

            # Same scan coordinate printed again; do not create a new point.
            if abs(angle - current_angle) < 1.0e-6:
                continue

            # New scan point begins: finalize previous point.
            if last_energy is not None:
                records.append((current_angle, last_energy))

            current_angle = angle
            last_energy = None

    if current_angle is not None and last_energy is not None:
        records.append((current_angle, last_energy))

    # Remove consecutive duplicates, keeping the last energy for the angle.
    cleaned = []
    for angle, energy in records:
        if cleaned and abs(cleaned[-1][0] - angle) < 1.0e-6:
            cleaned[-1] = (angle, energy)
        else:
            cleaned.append((angle, energy))

    if not cleaned:
        raise RuntimeError(
            "Parsed zero scan points. Check that the log contains a relaxed ModRedundant scan "
            "and lines containing both 'SCF Done' and 'Scan'."
        )

    emin = min(e for _, e in cleaned)
    return [(angle, (energy - emin) * HARTREE_TO_KCAL) for angle, energy in cleaned]

def main():
    args = parse_args()
    lines = read_text(args.file_log)

    if args.ftype.lower() != "qm":
        raise RuntimeError("This patched log2scan.py only handles -t qm. Use get_mm_energy.py for MM logs.")

    check_normal_termination(lines, args.file_log)
    rows = parse_qm_scan(lines)

    with open(args.output, "w") as f:
        for angle, rel_kcal in rows:
            f.write(f"{angle:.8f},{rel_kcal:.10f}\n")

    print(f"Wrote {args.output}")
    print(f"Scan points: {len(rows)}")
    print("Energy unit: relative kcal/mol")

if __name__ == "__main__":
    main()
