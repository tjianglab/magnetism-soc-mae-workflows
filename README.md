# Magnetism SOC MAE Workflows

Atomate VASP workflow templates for high-throughput screening of magnetic
materials. These workflows enumerate FM/AFM magnetic orderings, run structural
relaxation and static calculations, and add spin-orbit coupling (SOC)
calculations for magnetic anisotropy energy (MAE) analysis.

## Scope

This repository currently contains two workflow variants:

- `hex_trig_tet_magnetism_soc_mae.py`
  - For hexagonal, trigonal, and tetragonal crystal systems.
  - Supports `standard_cell="primitive"` and `standard_cell="conventional"` for
    the FM SOC branch.
  - Samples symmetry-relevant SOC magnetization directions using `SAXIS`.

- `orthorhombic_magnetism_soc_mae.py`
  - For orthorhombic crystal systems.
  - Defaults to `standard_cell="conventional"`.
  - Compares SOC energies along the three principal axes: `[100]`, `[010]`, and
    `[001]`.

- `random_afm_hex_trig_tet_magnetism_soc_mae.py`
  - For large hexagonal, trigonal, and tetragonal magnetic cells where
    symmetry-based magnetic-ordering enumeration is impractical.
  - Uses `random_afm_generator.py` to generate one FM reference and multiple
    random zero-net-moment AFM trial configurations.
  - The FM/AFM energy difference can be used as input for magnetic stability
    analysis or approximate Curie-temperature trend estimates.

- `random_afm_generator.py`
  - Standalone pymatgen helper that splits each magnetic species into spin-up
    and spin-down groups.
  - If a magnetic species has an odd number of sites, it doubles one cell axis
    before creating AFM configurations.

## Intended Use

These files are designed as drop-in workflow templates for an atomate VASP
installation. In the working atomate environment, copy the desired workflow file
to:

```text
atomate/vasp/workflows/base/magnetism_soc.py
```

For example:

```bash
cp atomate/vasp/workflows/base/hex_trig_tet_magnetism_soc_mae.py \
   atomate/vasp/workflows/base/magnetism_soc.py
```

Then use the workflow class in the same way as the customized
`magnetism_soc.py` workflow.

For the random-AFM workflow, also make sure `random_afm_generator.py` is
importable by Python. The simplest option is to place it in the same atomate
workflow package:

```text
atomate/vasp/workflows/base/random_afm_generator.py
```

Alternatively, keep it elsewhere and add that directory to `PYTHONPATH`.

## Example

```python
from atomate.vasp.workflows.base.magnetism_soc import MagneticOrderingsSOCWF

wf = MagneticOrderingsSOCWF(structure).get_wf(
    num_samples=7,
    nbands_factor=2,
    kpoints_factor=2,
    standard_cell="primitive",
)
```

For the orthorhombic workflow:

```python
wf = MagneticOrderingsSOCWF(structure).get_wf(
    nbands_factor=2,
    kpoints_factor=2,
    standard_cell="conventional",
)
```

For the random-AFM workflow:

```python
from atomate.vasp.workflows.base.magnetism_soc import MagneticOrderingsSOCWF

wf = MagneticOrderingsSOCWF(
    structure,
    random_afm_kwargs={
        "num_afm": 6,
        "seed": 123,
    },
).get_wf(
    num_samples=7,
    nbands_factor=2,
    kpoints_factor=2,
    standard_cell="primitive",
)
```

## Workflow Summary

1. Enumerate candidate collinear magnetic orderings using pymatgen, or generate
   random AFM candidates for large magnetic cells.
2. Run VASP relaxation and static calculations for each ordering.
3. Attach SOC calculations to the FM-derived branch.
4. Use the resulting SOC total-energy differences to analyze MAE.
5. Store magnetic ordering metadata through atomate/FireWorks database tasks.

## Requirements

These workflows assume a working atomate VASP environment with:

- atomate
- FireWorks
- pymatgen
- VASP execution configured through atomate
- a configured atomate database file

## Notes

- The workflows provide the calculation setup for MAE. Final MAE values should
  be obtained by comparing the converged SOC total energies for different
  `SAXIS` directions.
- K-point density, `NBANDS` scaling, and INCAR settings should be convergence
  tested for production studies.
- These workflows were developed for high-throughput permanent-magnet screening
  where both magnetic ordering stability and SOC anisotropy are important.

## License

This project is released under the MIT License.
