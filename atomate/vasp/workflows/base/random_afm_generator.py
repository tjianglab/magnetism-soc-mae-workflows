# coding: utf-8
"""
Random FM/AFM magnetic-configuration generator for large magnetic cells.

This helper is designed for cases where pymatgen's symmetry-based
MagneticStructureEnumerator becomes impractical because the structure contains
many magnetic atoms or many magnetic species. It creates one FM configuration
and several random zero-net-moment AFM trial configurations by splitting the
sites of each magnetic species into spin-up and spin-down groups.

Typical use:
    from random_afm_generator import MyCustomClass

    generator = MyCustomClass(structure, num_afm=6, seed=123)
    structures, magmoms = generator.get_fm_afm_struct()

The class name MyCustomClass is kept for compatibility with earlier local
workflow scripts that imported ``from random_afm import MyCustomClass``.
"""

import random

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


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


class MyCustomClass:
    """Generate one FM and multiple random AFM structures.

    Args:
        structure: input pymatgen Structure.
        num_afm: number of random AFM trial configurations to generate.
        magmoms: optional element-to-moment dictionary. Defaults to
            DEFAULT_MAGMOMS.
        seed: optional random seed for reproducible AFM configurations.
        symprec: symmetry tolerance used by SpacegroupAnalyzer.
        afm_standard_cell: standardized cell used for AFM trials. The default
            is "primitive".
        fm_standard_cell: standardized cell used for the FM reference. The
            default is "conventional".
        odd_count_supercell_axis: axis doubled when any magnetic species has an
            odd site count. Use "a", "b", "c", or "auto". For hexagonal,
            trigonal, and tetragonal cells, "auto" usually doubles c because
            a and b are the most similar lattice lengths.
    """

    def __init__(
        self,
        structure: Structure,
        num_afm=5,
        magmoms=None,
        seed=None,
        symprec=0.1,
        afm_standard_cell="primitive",
        fm_standard_cell="conventional",
        odd_count_supercell_axis="auto",
    ):
        self.original_structure = structure
        self.num_afm = num_afm
        self.magmoms = dict(DEFAULT_MAGMOMS if magmoms is None else magmoms)
        self.rng = random.Random(seed)
        self.symprec = symprec
        self.afm_standard_cell = afm_standard_cell
        self.fm_standard_cell = fm_standard_cell
        self.odd_count_supercell_axis = odd_count_supercell_axis

        self.structure = self._standardize_structure(structure, afm_standard_cell)
        if len(self.structure.sites) == 0:
            raise ValueError("After applying symmetry, the structure is empty.")

    def print_struct_sites(self):
        """Return sites for compatibility with older workflow scripts."""
        return self.structure.sites

    def get_fm_afm_struct(self):
        """Return ``([fm, afm1, ...], [fm_magmom, afm1_magmom, ...])``."""
        afm_base_structure = self._make_even_magnetic_counts(self.structure)
        afm_structures = []
        afm_magmoms = []

        for _ in range(self.num_afm):
            magmom = self._make_random_afm_magmom(afm_base_structure)
            afm_structure = afm_base_structure.copy()
            afm_structure.add_site_property("magmom", magmom)
            afm_structures.append(afm_structure)
            afm_magmoms.append(magmom)

        fm_structure = self._standardize_structure(self.original_structure, self.fm_standard_cell)
        fm_magmom = self._make_fm_magmom(fm_structure)
        fm_structure.add_site_property("magmom", fm_magmom)

        return [fm_structure] + afm_structures, [fm_magmom] + afm_magmoms

    def _standardize_structure(self, structure, standard_cell):
        spa = SpacegroupAnalyzer(structure, symprec=self.symprec)
        if standard_cell == "primitive":
            return spa.get_primitive_standard_structure()
        if standard_cell == "conventional":
            return spa.get_conventional_standard_structure()
        if standard_cell == "input":
            return structure.copy()
        raise ValueError(
            "standard_cell must be 'primitive', 'conventional', or 'input', "
            f"not {standard_cell!r}."
        )

    def _make_even_magnetic_counts(self, structure):
        counts = self._magnetic_site_indices_by_element(structure)
        has_odd_count = any(len(indices) % 2 for indices in counts.values())
        if not has_odd_count:
            return structure.copy()

        afm_structure = structure.copy()
        afm_structure.make_supercell(self._supercell_matrix_for_odd_counts(structure))
        return afm_structure

    def _supercell_matrix_for_odd_counts(self, structure):
        axis = self.odd_count_supercell_axis
        if axis == "auto":
            axis = self._auto_supercell_axis(structure)

        matrices = {
            "a": [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
            "b": [[1, 0, 0], [0, 2, 0], [0, 0, 1]],
            "c": [[1, 0, 0], [0, 1, 0], [0, 0, 2]],
        }
        if axis not in matrices:
            raise ValueError("odd_count_supercell_axis must be 'a', 'b', 'c', or 'auto'.")
        return matrices[axis]

    def _auto_supercell_axis(self, structure):
        a, b, c = structure.lattice.a, structure.lattice.b, structure.lattice.c
        diff_ab = abs(a - b)
        diff_bc = abs(b - c)
        diff_ca = abs(c - a)
        if diff_ab <= diff_bc and diff_ab <= diff_ca:
            return "c"
        if diff_bc <= diff_ab and diff_bc <= diff_ca:
            return "a"
        return "b"

    def _magnetic_site_indices_by_element(self, structure):
        element_indices = {}
        for index, site in enumerate(structure.sites):
            element = str(site.specie)
            if element in self.magmoms:
                element_indices.setdefault(element, []).append(index)
        return element_indices

    def _make_fm_magmom(self, structure):
        return [self.magmoms.get(str(site.specie), 0.0) for site in structure.sites]

    def _make_random_afm_magmom(self, structure):
        magmom = [0.0] * len(structure.sites)
        for element, indices in self._magnetic_site_indices_by_element(structure).items():
            shuffled_indices = list(indices)
            self.rng.shuffle(shuffled_indices)
            half = len(shuffled_indices) // 2
            spin_up = set(shuffled_indices[:half])
            spin_down = set(shuffled_indices[half:])
            moment = self.magmoms[element]

            for index in spin_up:
                magmom[index] = moment
            for index in spin_down:
                magmom[index] = -moment

        return magmom
