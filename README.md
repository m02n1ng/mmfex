
# M-MFEx : Multiconfigurational Hamiltonian for Frenkel Exciton

The package uses OpenMolcas (v24.06 tag: 12-g13ccb99e8) outputs to create a Frenkel Matrix (limited to the first excitation manifold). With the Hamiltonian, it solves the time-undependent Schr&ouml;dinger Equation to compute the photo-physical properties of the system.


## Installation Guide

It is suggested to create a virtual environment by using `venv` or `conda`.

### `venv`
* Git Clone the repository running
```
git clone https://github.com/m02n1ng/Frenkel_Hamiltonian.git
```
* inside the repositopry folder, create and activate the python environment
```
python -m venv mmfex
source mmfex/bin/activate
```
* then run to install the dependencies
```
pip install -e .
```

### `conda`

* Git Clone the repository running
```
git clone https://github.com/m02n1ng/Frenkel_Hamiltonian.git
```
* create and activate the conda environment
```
conda create -n mmfex
conda activate mmfex
```
* then install the dependencies
```
pip install -e .
```


The package can be used also without any installation, if the module `actions` is in the current folder.

## User Guide

#### 1. Create a directory `mother_directory` that will contains all the couplings subdirectories
Once created the folder, copy all the _**insert here files needed to run couplings**_ files of all the monomers. It is important that the monomer data have a prograssive labeling.  
Then one can run the _**insert name of the code that creates and run couplings inputs**_ to automatically create the subdirectories and run the OpenMolcas couplings inputs.

#### 2. Check the end of the calculations
_**THERE SHOULD BE AN INPUT FILE MADE BY LUCA**_

#### 3. Elaborate the data
Within the `mother_directory`, one can use the `actions` module to extract the data and build the Frenkel-Davydov Hamiltonian. See the `testing`folder for an example.



