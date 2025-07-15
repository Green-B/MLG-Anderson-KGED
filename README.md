# MLG Anderson KGED

This code calculates the Kubo-Greenwood conductivity of monolayer graphene (MLG) with Anderson disorder.
We use a tight-binding model for MLG and obtain its eigenenergies and velocity matrix elements through exact diagonzlation.
By substituting the Hamiltonian and velocity matrix elements in `mlganderson-new-1-diag.py`, this code can be adapted to other systems as well.

Each of the essential Python scripts here except the JuPyter notebook `mlganderson-new-4-calculate-sigma.ipynb` is somewhat computationally intensive and is intended to be used in HPC settings with command line arguments as described in the codes.
The shell scripts listed below are examples of how the essential scripts may be submitted to an HPC cluster.
The notebook can readily be used on a standard computer with arguments set in the notebook during use.
Files named in the pattern `mlganderson-new-#-...` are used to generate and process data and are in principle all that is necessary to obtain results.
Files named in the pattern `mlganderson-new-x#-...` are for convenience due to adding data after a first complete run, making use of data from a partially-completed run, and similar circumstances.

## Overview - theory

With $\hbar=1$ and $e=1$, the real part of the longitudinal conductivity as obtained from the Kubo formula is

$$
\begin{aligned}
\Re\left( \sigma_{xx}(\omega) \right)
& = \Re\left( \lim_{\eta \to 0} -\frac{i}{V} \sum_{n,m} \frac{f(E_n)-f(E_m)}{E_n-E_m} \frac{\left|v_{nm}\right|^2}{E_n-E_m+\omega+i\eta} \right) \\
& = \frac{\pi}{V} \sum_{n,m} \frac{f(E_n)-f(E_m)}{E_n-E_m} \left|v_{nm}\right|^2 \delta(E_n-E_m+\omega) \\
& = \frac{\pi}{V} \sum_{n,m} \left( \int_{-\infty}^{\infty} dE \int_{-\infty}^{\infty} dE' \delta(E-E_n)\delta(E'-E_m) \right) \frac{f(E_n)-f(E_m)}{E_n-E_m} \left|v_{nm}\right|^2 \delta(E_n-E_m+\omega) \\
& = \frac{\pi}{V} \int_{-\infty}^{\infty} dE \int_{-\infty}^{\infty} dE' \left( \sum_{n,m} \left|v_{nm}\right|^2  \delta(E-E_n)\delta(E'-E_m) \right) \frac{f(E)-f(E')}{E-E'} \delta(E-E'+\omega) \\
& = \frac{\pi}{V} \int_{-\infty}^{\infty} dE \int_{-\infty}^{\infty} dE' J\left(E,E'\right) \frac{f(E)-f(E')}{E-E'} \delta(E-E'+\omega)
\end{aligned}
$$

where we have defined

$$
J\left(E,E'\right) = \sum_{n,m} |v_{nm}|^2  \delta(E-E_n)\delta(E'-E_m)
$$

Now transform variables to $s=E-E'$, $\varepsilon = \frac{E+E'}{2}$ with the Jacobian $J = \frac{\partial s}{\partial E}\frac{\partial \varepsilon}{\partial E'} - \frac{\partial s}{\partial E'} \frac{\partial \varepsilon}{\partial E} = 1$ and the previous expression becomes

$$
\begin{aligned}
\Re\left( \sigma_{xx}(\omega) \right)
& = \frac{\pi}{V} \int_{-\infty}^{\infty} d\varepsilon \int_{-\infty}^{\infty} ds J\left(\varepsilon+\frac{s}{2},\varepsilon-\frac{s}{2}\right) \frac{f(\varepsilon+\frac{s}{2})-f(\varepsilon-\frac{s}{2})}{s} \delta(s+\omega) \\
& = \frac{\pi}{V} \int_{-\infty}^{\infty} d\varepsilon J\left(\varepsilon+\frac{s}{2},\varepsilon-\frac{s}{2}\right) \frac{f(\varepsilon+\frac{s}{2})-f(\varepsilon-\frac{s}{2})}{s}
\end{aligned}
$$

The factor $J\left(\varepsilon+\frac{s}{2},\varepsilon-\frac{s}{2}\right)$ is the subject of `mlganderson-new-2-calc-j.py` and `mlganderson-new-3-combine-j.py`, where the Dirac deltas are treated by binning them into divisions in the space of $\varepsilon$ and $s$.
The integration over $\varepsilon$ is carried out in `mlganderson-new-4-calculate-sigma.ipynb`.

## Overview - code

### Essential scripts:

Each of the essential scripts requires the output files of the preceding script and no other input files.

- `mlganderson-new-1-diag.py` creates models of disordered MLG and calculates their energies and x-direction velocity matrix elements. y-direction velocity is available in the code as well.
	- Command line arguments:
		1. Anderson disorder strength $w$ such that site energies are in $[-w,w]$
		2. One-third the number of unit cells along each axis of the graphene lattice (i.e., $N$ such that the lattice is $3N×3N$; multiples of three ensure that the K-points are sampled)
		3. Number of k-points $N_k$ along each axis; if $N_k>1$ then periodic boundary conditions are used, while if $N_k=1$ open boundary conditions are used
		4. Number of disorder configurations to simulate
	- Files produced:
		1. `en_v2/e-v-w#-n#-pbc#u#v#-rep#.pkl` - energies and velocity matrix elements with the file labeled by disorder, system size, periodic boundary conditions, and configuration index
		2. `en_v2/e-min-max-w#-n#.pkl.pkl` - analytic bounds on minimum and maximumm energies with the file labeled by disorder and system size

- `mlganderson-new-2-calc-j.py` calculates a factor $J$ for a single configuration or realization of disordered graphene. $J$ is based on energies and velocity matrix elements and contains all information about the system necessary to calculate the conductivity later.
	- Command line arguments:
		1. Anderson disorder strength $w$ such that site energies are in $[-w,w]$
		2. Number of bins to use in histogramming energy averages
		3. Number of bins to use in histogramming energy differences
	- Files produced:
		1. `jes_individual/jes-w%s-ne%s-ns%s.pkl` - $J$ for this configuration with the file labeled by disorder, the two numbers of histogram bins, and configuration index

- `mlganderson-new-3-combine-j.py` combines the results of many configurations so as to average over them.
	- Command line arguments:
		1. Anderson disorder strength $w$ such that site energies are in $[-w,w]$
		2. Number of bins to use in histogramming energy averages
		3. Number of bins to use in histogramming energy differences
	- Files produced:
		1. `jes-combined-w%s-ne%s-ns%s-n%s-nc%s.pkl` - $J$ combined for many configurations with the file labeled by disorder, the two numbers of histogram bins, the number of configurations aggregated in this data

- `mlganderson-new-4-calculate-sigma.ipynb` uses the averaged $J$ to calculate the conductivity.
	- Command line arguments: None, input should be given in the notebook
	- Files produced: None; figures and numerical data are produced in the notebook

### Utility scripts:

Each of the utility scripts inputs and outputs files of the same type, e.g., a script that takes a `jes-combined` file as input will produce a `jes-combined` file as output.

- `mlganderson-new-x1-jes-reduce-ne-ns.py` combines bins in histogram used to calculate $J$, allowing for smoother averaging at the cost of lower resolution in frequency in the calculated conductivity.
	- Command line arguments:
		1. Anderson disorder strength $w$ such that site energies are in $[-w,w]$
		2. Number of bins previously used in histogramming energy averages
		3. Number of bins previously used in histogramming energy differences
		4. Total number of unit cells in the supercell
		5. Number of disorder configurations previously simulated
		6. Factor by which to reduce the number of energy average bins by combining bins in groups of as many
		7. Factor by which to reduce the number of energy difference bins by combining bins in groups of as many
	- Input and output files:
		1. `jes-combined-w%s-ne%s-ns%s-n%s-nc%s.pkl` - $J$ combined for many configurations with the file labeled by disorder, the two numbers of histogram bins, the number of configurations aggregated in this data

- `mlganderson-new-x2-jes-repshift.py` shifts the indexing of MLG disorder configurations. `mlganderson-new-3-combine-j.py` is designed to work on data from a single run of `mlganderson-new-2-calc-j.py` which labels output configurations starting with 1. Therefore, this script makes it easy to combine results of multiple runs of `mlganderson-new-2-calc-j.py` with `mlganderson-new-3-combine-j.py`.
	- Command line arguments:
		1. Anderson disorder strength $w$ such that site energies are in $[-w,w]$
		2. Number of bins previously used in histogramming energy averages
		3. Number of bins previously used in histogramming energy differences
		4. Number by which to shift disorder configuration/repetition numbers
	- Input and output files:
		1. `jes_individual/jes-w%s-ne%s-ns%s.pkl` - $J$ for this configuration with the file labeled by disorder, the two numbers of histogram bins, and configuration index

- `mlganderson-new-x3-combine-jes-combined.py` allows two copies of the averaged $J$ outputted by `mlganderson-new-3-combine-j.py` to be combined. It can serve as another way to combine results from multiple runs.
	- Command line arguments:
		1. Number of configurations used in jes-combined file #1
		2. User-written suffix in filename identifying jes-combined file #1
		3. Number of configurations used in jes-combined file #2
		4. User-written suffix in filename identifying jes-combined file #2
		4. User-written suffix in filename to identify output jes-combined file that combines files #1 and #2
	- Input and output files:
		1. `jes-combined-w%s-ne%s-ns%s-n%s-nc%s.pkl` - $J$ combined for many configurations with the file labeled by disorder, the two numbers of histogram bins, the number of configurations aggregated in this data

### Shell scripts:

- `mlganderson-sub-to-nodes.sh`
	- Command line arguments: 
		1. Anderson disorder strength $w$ such that site energies are in $[-w,w]$
		2. One-third the number of unit cells along each axis of the graphene lattice (i.e., $N$ such that the lattice is $3N×3N$; multiples of three ensure that the K-points are sampled)
		3. Number of disorder configurations to be simulated per node submission
		4. Number of node submissions to schedule

- `mlganderson-node-script.sh`
	- Command-line arguments: None; this script is run on nodes with input provided by `mlganderson-sub-to-nodes.sh`
