# Test Marks

This test suite uses custom `pytest` marks so changes can be validated with a focused subset of tests instead of always running the full regression module.

## Mark Meanings

- `slow`
  End-to-end regression tests that take noticeable time to run.

- `smoke`
  Very fast checks such as import or basic collection sanity.

- `validation`
  Input validation and expected failure/`SystemExit` coverage.

- `psi4`
  Tests that generate Psi4 inputs.

- `qchem`
  Tests that generate Q-Chem inputs.

- `nwchem`
  Tests that generate NWChem inputs.

- `amber`
  Tests that exercise the Amber preparation workflow.

- `charmm`
  Tests that exercise the CHARMM preparation workflow.

- `sapt`
  SAPT or FISAPT-style workflows.

- `hf`
  Hartree-Fock supermolecular workflows.

- `brc`
  Balanced redistribution charge schemes, including BRC variants like `BRC2` and `BRCD`.

- `dz`
  Distributed-zeroing charge schemes, including `DZ1`, `DZ2`, and `DZ3`.

- `z`
  Zeroing charge schemes, including `Z1`, `Z2`, and `Z3`.

- `see`
  Static electrostatic embedding charge scheme.

- `template`
  Template-driven QM region conversion workflow.

## Common Commands

Run only very fast checks:

```bash
pytest -m smoke
```

Run only input validation tests:

```bash
pytest -m validation
```

Run only Psi4 + CHARMM + HF coverage:

```bash
pytest -m "psi4 and charmm and hf"
```

Run only Amber SAPT workflows:

```bash
pytest -m "amber and sapt"
```

Run only template conversion coverage:

```bash
pytest -m template
```

Skip long-running tests:

```bash
pytest -m "not slow"
```

Run only a charge scheme family:

```bash
pytest -m dz
pytest -m z
pytest -m brc
pytest -m see
```

## Notes

- Marks are registered in `pyproject.toml`.
- Multiple marks can be combined with `and`, `or`, and `not`.
- The main regression module is `sparcleqc/tests/test_sparcle_qc.py`.
