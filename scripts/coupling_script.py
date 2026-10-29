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

import os
import shutil

m = <NUMBER OF THE FIRST MONOMER> #first monomer
n = <NUMBER OF THE LAST MONOMER> #last monomer

#all the files are into the current work directory
source = os.getcwd()

#the cycle starts here
for i in range(m, n+1):
  alfa = str(i)
  for j in range(i+1, n+1):
    beta = str(j)
    
#   define the destination directory for the files
    dest = source + '/coup_' + alfa + '_' + beta
    
#   create the directory for the coupling
    os.mkdir(dest)
    
#   copy the files from the source to the destination with new names
    shutil.copy(source + '/mon_' + alfa + '.xyz', dest + '/monoA.xyz')
    shutil.copy(source + '/mon_' + beta + '.xyz', dest + '/monoB.xyz')
    shutil.copy(source + '/mon_' + alfa + '.JobMix', dest + '/monoA.JobMix')
    shutil.copy(source + '/mon_' + beta + '.JobMix', dest + '/monoB.JobMix')
    shutil.copy(source + '/coupling.input', dest + '/coupling_' + alfa + '_' + beta + '.input')
    
#   submits the calculation
    os.chdir(dest)
    os.popen(<FILL HERE WITH YOUR USUAL COMMAND LINE RUN>)
    os.chdir(source)
    
    
