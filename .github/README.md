# Homeostatic Economics — numerical demonstrator (V1.1)

Reference engine and reproducible package accompanying the article:

> Benjamin Albouy, *Mise en calcul d'une architecture biophysique–monétaire sous contraintes : un démonstrateur numérique exploratoire* (*Computational instantiation of a constrained biophysical–monetary architecture: an exploratory numerical demonstrator*), V1.1, November 2026.

| Document | File |
|---|---|
| Article, French original (authoritative text) | [`docs/EH-Article-demonstrateur-V1.1.pdf`](../docs/EH-Article-demonstrateur-V1.1.pdf) |
| Article, English translation (for information) | [`docs/EH-Demonstrator-Article-V1.1-EN.pdf`](../docs/EH-Demonstrator-Article-V1.1-EN.pdf) |
| Package user guide (French) | [`README.md`](../README.md) |
| Changes from V1 to V1.1 (French) | [`CHANGELOG_V1.1.md`](../CHANGELOG_V1.1.md) |

The demonstrator computes one particular instance of a coupled biophysical–monetary dynamical system under viability constraints (three biophysical indicators aggregated by a geometric mean, a money stock with primary issuance modulated by a regulation coefficient and a demurrage term, and a possibly delayed prescribed feedback), integrated by an explicit Euler scheme over a finite horizon. It belongs to the Homeostatic Economics research programme: https://vteh.eu

## Contents

| Item | Role |
|---|---|
| `eh_reproductible.py` | Reference engine. One run regenerates the ten tables, the CSV files and the figures of the article. |
| `csv/` | Ten reference tables, T1 to T10. |
| `figures/` | The eight figures of the article, as PNG (300 dpi) and vector PDF. |
| `execution.log`, `EXECUTION_REFERENCE.txt` | Standard output, environment and fingerprint of the reference run. |
| `controle/` | Independent check with higher-order integrators, its dependencies and expected output. |
| `tests/` | Engine tests and package non-regression test, with expected outputs. |
| `MANIFEST_SHA256.txt` | SHA-256 fingerprints of all package files. |

SHA-256 of the engine, identical to the one printed in §13 of the article:

    ffc2339bcb2b61adb5a0c71f6d30f063636d1015a3db547e36e7bd62f89bd9b5

## Requirements

Python 3.8 or later. The engine, its tables, its CSV files and the tests use the standard library only. Optional: `matplotlib` for the figures; `numpy` and `scipy` for the check script.

    python3 -m pip install matplotlib
    python3 -m pip install -r controle/requirements.txt

Reference environment: Python 3.12.3, matplotlib 3.10.8; check script: numpy 2.5.4, scipy 1.18.1.

## Quick start

**1. Check package integrity**

    sha256sum -c MANIFEST_SHA256.txt

Expected: `OK` for every file.

**2. Regenerate tables and figures** (the engine writes next to itself, so run a copy in an empty directory)

    mkdir essai && cp eh_reproductible.py essai/ && cd essai
    python3 eh_reproductible.py > execution.log

Expected: exit code 0; ten CSV files and a log byte-for-byte identical to those of the package; sixteen files in `figures/`.

**3. Non-regression test** (step 2 and the comparison, automated)

    python3 tests/test_non_regression.py

Expected: `BILAN : dix CSV et journal identiques octet pour octet`. Figure rendering may differ outside the reference environment without any change in the plotted data; the test reports this for information.

**4. Engine tests**

    python3 tests/tests_moteur.py

Expected: `BILAN : 111 tests conformes, 0 non conformes`.

**5. Higher-order integrator check**

    python3 controle/controle_ordre_superieur.py

Expected: `BILAN : tous les contrôles sont CONFORMES`. The script integrates the same equations outside the reference engine (DOP853, rtol = 1e-12, and the method of steps for delays).

Program outputs, identifiers and comments are in French.

## Reproducibility

The engine has no source of randomness: two runs in the same environment produce byte-identical files. Across environments, CSV files and the log depend on double-precision floating-point arithmetic, and figure rendering depends on the matplotlib version and available fonts; the plotted data are the same. The repository disables line-ending conversion (`.gitattributes`) so that the fingerprints in `MANIFEST_SHA256.txt` remain valid after cloning on any system.

## Licence

Code (`*.py`, `controle/requirements.txt`): MIT — see [`LICENSE`](../LICENSE). Data, figures, documents and articles: CC BY 4.0. Details in [`LICENSES.md`](../LICENSES.md).

## Citation

Please cite the article; GitHub's "Cite this repository" button uses [`CITATION.cff`](../CITATION.cff).
