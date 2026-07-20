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

## Workflow Summary

1. Enumerate candidate collinear magnetic orderings using pymatgen.
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
