# Computational supplements for the \{4,5\}-PBD block-number spectrum

This repository contains the computational supplements for the paper

> The complete block-number spectrum of \{4,5\}-pairwise balanced designs
> and Zarankiewicz numbers in the $p=4$ range for sufficiently large orders

The code checks the exact arithmetic behind the spectrum theorem, verifies
finite packing modules and certificates used in the transition analysis, and
generates the finite residue tables reported in the paper.

## Scope

The scripts verify finite and arithmetic parts of the arguments. They do not
replace the asymptotic dense-decomposition theorems used to complete the
bounded leaves for sufficiently large orders. In particular, a successful run
is a reproducibility check for the displayed calculations and certificates,
not an independent proof of the asymptotic decomposition input.

## Requirements

- Python 3.10 or later
- `ortools` for the SAT/CP-SAT searches and the scripts that import it

The remaining scripts use only the Python standard library. A typical setup is

```text
python -m pip install ortools
```

## Main checks

Run from the repository root:

```text
python scripts/compute_pbd_spectrum.py --csv data/spectrum_summary.csv --residue-table > data/residue_table.txt
python scripts/verify_complete_spectrum.py
python scripts/verify_transition_deficits.py
python scripts/verify_half_integral_signature_spectrum.py
python scripts/verify_dual_endpoint.py
python scripts/verify_endpoint_modules.py
python scripts/verify_dual_modules.py
python scripts/verify_dual_fixed_deficits.py
python scripts/verify_first_deficits.py
python scripts/verify_residue_bookkeeping.py
```

Scripts whose names begin with `search_`, together with
`verify_first_window_modules.py`, `analyze_dual_transition.py`, and the SAT
drivers, perform the more expensive finite searches. Their comments identify
the certificate or module family being checked.

## Data

The `data/` directory contains generated arithmetic summaries only. No private
manuscript files, temporary compilation products, cached solver files, or
external copyrighted PDFs are included.

## Authors

Yong Zhang and Guangzhou Chen
