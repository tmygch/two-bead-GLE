# Toy6

Simulation and theoretical-analysis codes for the six-bead elastic-network example in the main text of

**“Exact Generalized Langevin Dynamics of Pair Coordinates in Elastic Networks”**

The network is

```text
          6
          |
1 -- 2 -- 3 -- 4
     \   /
       5
```

## Files

- `toy6.cpp`  
  Brownian-dynamics simulation of the six-bead elastic network and numerical estimation of the memory kernel.

- `run.py`  
  Runs the tagged-pair cases used for the six-bead example.

- `theory_toy6.py`  
  Calculates the exact memory kernel from the matrix expression derived in the manuscript.

## Requirements

### Simulation

- C++17-compatible compiler
- `g++`

### Theoretical calculation

- Python 3
- NumPy
- SciPy

For example,

```bash
python3 -m pip install numpy scipy
```

## Simulation

For a single tagged-pair case, compile with

```bash
g++ -O3 -std=c++17 -DMMM1=0 -o TN=0.out toy6.cpp
```

and run with

```bash
./TN=0.out > a.dat
```

Alternatively, run the tagged-pair cases used in the main-text example with

```bash
python3 run.py
```

The run script creates a separate directory for each calculation and copies the C++ source file into that directory before compilation.

## Exact theoretical memory kernel

For the tagged pair (i, j) = (2, 6), run

```bash
python3 theory_toy6.py
```

The script evaluates

```math
\mu(t)
=
\tilde{\mathbf l}{}'_{i\bullet}
\cdot
e^{-\tilde L{}' t}\tilde L{}'^{-}
\cdot
\tilde{\mathbf l}{}'_{\bullet i}.
```

Here the superscript minus denotes the generalized inverse used in the manuscript.

The calculated data are written to files such as

```text
Toy6_N6_mu(t)_i002_j006_genth.dat
Toy6_N6_mu(s)_i002_j006_genth.dat
```

No plotting package is required. The figures in the manuscript were prepared separately.

## Reproducibility note

The original simulations were performed using the Numerical Recipes random-number routines `ran2` and `gasdev`.

These routines are not redistributed in this repository. The public version of `toy6.cpp` instead uses the C++ standard-library random-number generator. Therefore, the same nominal seed does not generate trajectories that are bitwise identical to those from the original simulations. The stochastic model and numerical integration scheme are otherwise unchanged.

## Citation

If you use this code, please cite the corresponding paper:

> T. Miyaguchi *et al.*, “Exact Generalized Langevin Dynamics of Pair Coordinates in Elastic Networks.”
