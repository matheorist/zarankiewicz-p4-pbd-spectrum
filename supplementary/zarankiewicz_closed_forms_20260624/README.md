# Supplementary computations for the (v=2) DGH relaxation

This directory contains the exact-arithmetic scripts and small result summaries
associated with the paper on closed forms for the (v=2) DGH relaxation in the
Zarankiewicz problem (z(m,n;3,3)).

The material covers four computational components:

- migrated supports across all residue classes;
- the progression bands used in the closed-form formulas;
- exact spot checks against the full finite LP;
- comparisons with one-cut relaxations, including two infinite gap families.

The scripts use rational arithmetic for the displayed identities. They verify
the algebraic formulas and the finite checks described in the manuscript; they
do not replace the written proofs or the graph-design realization problem.

## Requirements

- Python 3.10 or later
- SymPy 1.12 or later

Install the only nonstandard dependency with:

```text
python -m pip install "sympy>=1.12"
```

## Reproduce the main checks

Run the commands from this directory:

```text
python scripts/all_residue_migration_proof.py
python scripts/all_residue_exact_spotcheck.py
python scripts/companion_certificate_verify.py
python scripts/progression_r4_subprogression_verify.py
python scripts/progression_allbands_proof.py
python scripts/progression_deterministic_index_verify.py
python scripts/progression_allbands_lower.py
python scripts/single_cut_gain_proof.py
python scripts/one_cut_global_proof.py
python scripts/one_cut_unbounded_proof.py
```

The file `data/verification_summary.txt` records the successful local run
summaries, while `data/exact_spotchecks.txt` records five exact candidate-versus-
full-LP comparisons. The scripts are intentionally kept separate from the
manuscript source and its compilation products.

## Scope and provenance

No private review notes, temporary compilation files, screenshots, cached
solver files, or external PDFs are included. The numerical summaries are small
derived outputs, not an independent experimental data set.
