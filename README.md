
# M-MFEx : Multiconfigurational Hamiltonian for Frenkel Exciton

The package uses OpenMolcas (v24.06 tag: 12-g13ccb99e8) outputs to create a Frenkel Matrix (limited to the first excitation manifold). With the Hamiltonian, it solves the time-independent Schr&ouml;dinger Equation to compute the photo-physical properties of the system.


## Installation Guide

It is suggested to create a virtual environment using `venv` or `conda`. The package has been tested with Python 3.13.0 (the pinned `numpy` and `scipy` versions require Python 3.10 or newer).

### `venv`
* Git Clone the repository running
```
git clone https://github.com/m02n1ng/mmfex.git
```
* inside the repositopry folder, create and activate the python environment
```
cd mmfex
python -m venv .venv
source .venv/bin/activate
```
* then run to install the dependencies
```
pip install -e .
```

### `conda`

* Git Clone the repository running
```
git clone https://github.com/m02n1ng/mmfex.git
```
* create and activate the conda environment
```
conda create -n mmfex python=3.13
conda activate mmfex
```
* then install the dependencies
```
pip install -e .
```

## Testing

To check that your installation reproduces the reference results (a Frenkel dimer, see `testing/`), install the test extra and run pytest:
```
pip install -e ".[test]"
pytest
```

The package can be used also without any installation, if the module `actions` is in the current folder.

## User Guide

#### 1. Create a directory `mother_directory` that contains all the couplings subdirectories
Once created the folder, copy all the `.xyz` and the `.JobMix` files of all the monomers inside the `mother_directory`. Please notice that the script works correctly only if the monomer data are named `mon_<PROGRESSIVE NUMBER>.xyz` and `mon_<PROGRESSIVE NUMBER>.JobMix`. It is important that the monomer data have a progressive labeling.  
Copy also an Openmolcas input script for the coupling named `coupling.input`.  

Then one can run the script `coupling_script.py` to automatically create the subdirectories and run the OpenMolcas couplings inputs (notice that the jobs will run in parallel). 

#### 2. Check the end of the calculations
To have an update of the calculations use the script `check_completion.py`.

#### 3. Elaborate the data
Once all the coupling calculations are ended one can use the package to calculate the Frenkel Hamiltonian.  
Within the `mother_directory`, use the `actions` module to extract the data and build the Frenkel-Davydov Hamiltonian. See the `testing`folder for some examples.

If the `setup.py` has not been installed, use the `mmfex_gen.py` to generate an output file with the excitonic properties of the system. To make it work, copy and paste the `actions` module inside `mother_directory` to make it available.

## How to cite

M-MFEx is described in a forthcoming OpenMolcas publication. Citation details
will be added here upon publication. In the meantime, citation metadata for the
software is in [`CITATION.cff`](CITATION.cff).


