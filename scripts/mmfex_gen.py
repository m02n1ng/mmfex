#!/usr/bin/env python

from actions import *
import csv

m = <INITIAL MONOMER>
n = <LAST MONOMER>
sigma = <TYPE NUMBER HERE> # standard deviation of the gaussian function for the spectrum

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
                                                          printed_output=True)

x, y = get_spectrum(eigenvalues,
                    oscillators,
                    sigma,
                    csv = True)
