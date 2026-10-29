#!/usr/bin/env python

from actions import *
# from mmfex.actions import *

m = <INITIAL MONOMER>
n = <LAST MONOMER>
sigma = <TYPE NUMBER HERE>  # width of the gaussian for the spectrum, in eV;
                            # sqrt(2) times the standard deviation, FWHM ~ 1.665*sigma

total_roots = <TYPE NUMBER HERE> # total number of states used in the CASSCF calculation;
dim = <TYPE NUMBER HERE> # number of states considered in the Frenkel Hamiltonian per monomer;

C, E, TDM = collect(first_monomer=m, 
                    last_monomer=n, 
                    computed_roots=total_roots,
                    wanted_roots=dim)

Frenkel = build(E=E,
                C=C,
                units='eV')

eigenvalues, eigenvectors, oscillators, TDM_EXC = compute(Frenkel=Frenkel,
                                                          TDM=TDM,
                                                          units='eV',
                                                          printed_output=True)

x, y = get_spectrum(eigenvalues,
                    oscillators,
                    sigma,
                    csv = True)
