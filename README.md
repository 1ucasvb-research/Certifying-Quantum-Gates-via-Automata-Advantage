# Certifying Quantum Gates via Automata Advantage

Code used in the paper "[Certifying Quantum Gates via Automata Advantage](https://arxiv.org/abs/2510.09575)", by Anna Schroeder, Lucas B. Vieira, Jan Nöller, Nikolai Miklin, and Mariami Gachechiladze.

## QSQ robustness
The `QSQ_robustness_Sx.cpp` is the code used to generate the Quantum System Quizzing robustness data. Requires [Armadillo](https://arma.sourceforge.net/) and [tinyformat by Chris Foster](https://github.com/c42f/tinyformat) (included). Compiles directly with `g++` and `-std=c++17 -larmadillo` flags. This code is licensed LGPL v2.1

### Usage

Some parameters are hardcoded, but easy to modify at the top of the file:
- `DENSITY_PLOT_W`, `DENSITY_PLOT_H`: determine the width and height of the plot. Original values used were `4000` and `3000`, respectively.
- `MAX_PFAIL`, `MAX_INFID`: determine the range of probability of failure and infidelity covered by the plot. Original values `0.35f` and `0.15f`, respectively.
- `INFIDELITY_OPT_ATTEMPTS`, `INFIDELITY_OPT_ITERS`, `INFIDELITY_OPT_POWER`: control the unitary gauge optimization.
- `PROGRESS_PRINT`, `PROGRESS_SAVE`: how often to print and save progress.

After compiling, run with: `QSQ_robustness_Sx <L> <num_samples> <filename>`, where:
- `L`: integer specifying the size of the S-gate QSQ protocol, i.e., the sequences S^(2k), for k=0,...,L.
- `num_samples`: total number of noisy channels to sample for the plot.
- `filename`: a suffix to name the output file.

The output is a binary file named ``data_{DENSITY_PLOT_W}x{DENSITY_PLOT_H}_L{L}_{filename}.bin`` according to the settings above.
The file contains raw data of an array of `float`s of size `DENSITY_PLOT_W` x `DENSITY_PLOT_H`. Example Mathematica code to read it:

```Mathematica
BinaryReadList["data_4000x3000_L3_name.bin", "Real32"];
data = Image[ArrayReshape[data, {3000, 4000}]]
```

The code is single-threaded for simplicity, but can be ran in multiple instances with different target output files if necessary.
