#!/usr/bin/env python
# coding: utf-8
"""
Submit POSCAR files to a magnetic atomate workflow in ascending magnetic-site
count order.

This script scans the current directory for POSCAR-like files, counts magnetic
sites using the same default magnetic-moment table used by the random-AFM
generator, sorts structures from small to large magnetic complexity, and submits
FireWorks workflows to the configured LaunchPad.

The default settings are intended for the random-AFM workflow:
    - one FM reference
    - five random AFM trial configurations

Before using this script, copy
``random_afm_hex_trig_tet_magnetism_soc_mae.py`` to the active atomate
``magnetism_soc.py`` path and place ``random_afm_generator.py`` in the same
workflow package.

Use --dry-run first to check sorting and tags before submitting jobs.
"""

import argparse
from glob import glob

from atomate.vasp.powerups import add_modify_incar, add_tags
from atomate.vasp.workflows.base.magnetism_soc import MagneticOrderingsSOCWF
from fireworks import LaunchPad
from monty.serialization import dumpfn
from pymatgen.core import Structure


DEFAULT_MAGMOMS = {
    "Co": 5,
    "Co3+": 0.6,
    "Co4+": 1,
    "Cr": 5,
    "Fe": 5,
    "Mn": 5,
    "Mn3+": 4,
    "Mn4+": 3,
    "Mo": 5,
    "Ni": 5,
    "V": 5,
    "W": 5,
    "Ce": 5,
    "Eu": 10,
    "Ti3+": 1.73,
    "V3+": 2.83,
    "Cr3+": 3.88,
    "Cr2+": 4.9,
    "Mn2+": 5.92,
    "Fe3+": 5.92,
    "Co2+": 3.88,
    "Ni2+": 2.83,
    "Cu": 1.73,
    "Cu2+": 1.73,
    "Pr": 3.58,
    "Nd": 3.62,
    "Pm": 2.68,
    "Sm": 0.85,
    "Gd": 7.94,
    "Tb": 9.72,
    "Dy": 10.65,
    "Ho": 10.6,
    "Er": 9.58,
    "Tm": 7.56,
    "Yb": 4.54,
    "Ru": 2.2,
    "Os": 2.2,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Submit POSCAR files sorted by magnetic atom count."
    )
    parser.add_argument(
        "--pattern",
        default="*POSCAR",
        help="Glob pattern for input structures. Default: *POSCAR",
    )
    parser.add_argument(
        "--num-afm",
        type=int,
        default=5,
        help="Number of random AFM trial configurations. Default: 5",
    )
    parser.add_argument(
        "--num-samples",
        type=int,
        default=2,
        help="SOC angular sampling density for hex/trig/tet workflows. Default: 2",
    )
    parser.add_argument(
        "--nbands-factor",
        type=float,
        default=2,
        help="SOC NBANDS scaling factor. Default: 2",
    )
    parser.add_argument(
        "--kpoints-factor",
        type=float,
        default=4,
        help="SOC k-point reciprocal-density scaling factor. Default: 4",
    )
    parser.add_argument(
        "--standard-cell",
        choices=("primitive", "conventional"),
        default="primitive",
        help="Standardized cell used for the FM SOC branch. Default: primitive",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible AFM generation. Default: None",
    )
    parser.add_argument(
        "--ismear",
        type=int,
        default=0,
        help="ISMEAR passed through updated_user_incar_settings. Default: 0",
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=0.05,
        help="SIGMA passed through updated_user_incar_settings. Default: 0.05",
    )
    parser.add_argument(
        "--error-file",
        default="error_dict.json",
        help="Output JSON file for structures that fail workflow creation.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print sorted structures without submitting workflows.",
    )
    return parser.parse_args()


def count_magnetic_atoms(structure, magmoms=DEFAULT_MAGMOMS):
    return sum(1 for site in structure if site.specie.symbol in magmoms)


def tag_from_filename(filename):
    if filename.endswith(".POSCAR"):
        return filename[: -len(".POSCAR")]
    return filename


def load_and_sort_structures(pattern):
    struct_and_tag = []
    for filename in sorted(glob(pattern)):
        structure = Structure.from_file(filename)
        tag = tag_from_filename(filename)
        magnetic_atom_count = count_magnetic_atoms(structure)
        struct_and_tag.append((structure, tag, magnetic_atom_count, filename))

    return sorted(struct_and_tag, key=lambda item: (item[2], item[3]))


def build_workflow(structure, tag, args):
    wf_orderings = MagneticOrderingsSOCWF(
        structure,
        random_afm_kwargs={
            "num_afm": args.num_afm,
            "seed": args.seed,
        },
    ).get_wf(
        num_samples=args.num_samples,
        nbands_factor=args.nbands_factor,
        kpoints_factor=args.kpoints_factor,
        num_orderings_hard_limit=args.num_afm + 1,
        standard_cell=args.standard_cell,
        updated_user_incar_settings={
            "ISMEAR": args.ismear,
            "SIGMA": args.sigma,
        },
    )
    wf = add_tags(wf_orderings, [tag])
    return add_modify_incar(wf)


def main():
    args = parse_args()
    sorted_structures = load_and_sort_structures(args.pattern)

    if not sorted_structures:
        raise ValueError(f"No input structures matched pattern: {args.pattern}")

    for _, tag, magnetic_atom_count, filename in sorted_structures:
        print(f"{filename}: tag={tag}, magnetic_atoms={magnetic_atom_count}")

    if args.dry_run:
        return

    lpad = LaunchPad.auto_load()
    error_dict = {}

    for structure, tag, magnetic_atom_count, filename in sorted_structures:
        try:
            wf = build_workflow(structure, tag, args)
            lpad.add_wf(wf)
            print(f"Submitted {tag} ({magnetic_atom_count} magnetic atoms)")
        except Exception as exc:
            print(f"Error with {tag}: {exc}")
            error_dict[tag] = {
                "filename": filename,
                "magnetic_atom_count": magnetic_atom_count,
                "error": str(exc),
            }

    dumpfn(error_dict, args.error_file)


if __name__ == "__main__":
    main()
