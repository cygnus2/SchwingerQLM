# Spin-1 Quantum-Link Schwinger Model

Exact-diagonalization and matrix-product-state implementations of the
spin-1 quantum-link Schwinger model.

## Model

The physical states satisfy

\[
G_n = E_n-E_{n-1}-N_n+q_n=0, \qquad q_n=\frac{1-(-1)^n}{2}.
\]

## Implementations

- `EDcodes/`: exact gauge-invariant basis construction and sparse ED.
- `MPS/`: TeNPy-based MPS/DMRG implementation.

## Boundary conditions

The ED implementation supports:

- OBC with \(L-1\) dynamical links and fixed external fluxes.
- PBC with \(L\) matter sites and \(L\) links.

## Requirements

Python, SciPy, Jupyter, TeNPy, needed. C++, MKL, and Boost programmes exist, but not required.

## Running the ED notebook

```bash
jupyter lab EDcodes/spin1_qlm_ed.ipynb
