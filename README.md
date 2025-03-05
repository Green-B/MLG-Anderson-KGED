# MLG Anderson KGED

This code calculates the Kubo-Greenwood conductivity of monolayer graphene (MLG) with Anderson disorder. We use a tight-binding model for MLG and obtain its eigenenergies and velocity matrix elements through exact diagonzlation.

Each of the scripts here except the JuPyter notebook `mlganderson-new-4-calculate-sigma.ipynb` is designed to be used in HPC settings with command line arguments as described in the codes. The notebook can readily be used on a standard computer with arguments set in the notebook during use. Files named in the pattern `mlganderson-new-#-...` are used to generate and process data and are in principle all that is necessary to obtain results. Files named in the pattern `mlganderson-new-x#-...` are for convenience due to adding data after a first complete run, making use of data from a partially-completed run, and other such circumstances.

Essential scripts:
- `mlganderson-new-1-diag.py` creates models of disordered MLG and calculates their energies and velocity matrix elements.
- `mlganderson-new-2-calc-j.py` calculates a factor $J$ for a single configuration or realization of disordered graphene. $J$ is based on energies and velocity matrix elements and contains all information about the system necessary to calculate the conductivity later.
- `mlganderson-new-2-combine-j.py` combines the results of many configurations so as to average over them.
- `mlganderson-new-4-calculate-sigma.ipynb` uses the averaged $J$ to calculate the conductivity.

Utility scripts:
- `mlganderson-new-x1-jes-reduce-ne-ns.py` combines bins in histogram used to calculate $J$, allowing for smoother averaging at the cost of lower resolution in frequency in the calculated conductivity.
- `mlganderson-new-x2-jes-repshift.py` shifts the indexing of MLG disorder configurations. `mlganderson-new-2-combine-j.py` is designed to work on data from a single run of `mlganderson-new-2-calc-j.py` which labels output configurations starting with 1. Therefore, this script makes it easy to combine results of multiple runs of `mlganderson-new-2-calc-j.py` with `mlganderson-new-2-combine-j.py`.
- `mlganderson-new-x3-combine-jes-combined.py` allows two copies of the averaged $J$ outputted by `mlganderson-new-2-combine-j.py` to be combined. It can serve as another way to combine results from multiple runs.
