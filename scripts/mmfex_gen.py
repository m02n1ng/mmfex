#!/usr/bin/env python
# mmfex: Frenkel-Davydov exciton Hamiltonian from OpenMolcas outputs
# Copyright (C) 2026  Roberto Cacciari, Luca De Vico
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see
# <https://www.gnu.org/licenses/>.



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
