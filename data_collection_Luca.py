#!/usr/bin/env python

import os
import time
import shutil
import numpy as np
import re
import scipy.special 
from numpy import linalg as LA

def delta(a, b):
 if a == b:
  return 1
 elif a != b:
  return 0


m = 301 #first monomer
n = 327 #last monomer

A = n - m + 1 # number of molecules 
B = n - m #number of molecules minus one 
s = 2 # number of states considered in FRENKEL

total_s = 5 #number of computed roots
number_coup = int(scipy.special.comb(total_s + total_s, 2))

# C = [A][A_start][A_end][B][B_start][B_end]
C = np.zeros((A, s, s, A, s, s))

# E = [number of monomers][number of states]
E = np.zeros((A, s))

# TDM = [A, total number of states in RASSCF start, total number of states in RASSCF end, space component]
TDM = np.zeros((A,total_s, total_s, 3))
TOT_TDM = int(scipy.special.comb(total_s, 2))

# define the parsing keyword
key = ['Excitonic', 'couplings', 'section:']

#all the files are into the current work directory
source = os.getcwd()

#the cycle starts here
for i in range(m, n+1): 
 alfa = str(i)
 for j in range(i+1, n+1):
  beta = str(j)
# define the destination directory for the files
  dest = os.path.join(source,'coup_'+alfa+'_'+beta)
# define the path for the output file
  out = os.path.join(dest, 'coupling_'+alfa+'_'+beta+'.out')
# go to the parsing keyword
  with open(out, "r") as F:
   output = F.readlines()

 #COUPLING TERMS COLLECTION
 
#find the line number of the keyword 
  iline = 0
  for k  in range(len(output),0, -1): # searches from the end of the output file
   try:
    if output[k].split() == key:
     iline = k
     break
   except:
    pass
  iline += 4 # adds four lines before the start of the couplings
  for l in range(1, number_coup + 1):
   line = output[iline + l].split() 
   indeces = re.split(r'\(|\)', line[0])
   A_start = int(indeces[7])
   A_end = int(indeces[3])
   B_start = int(indeces[5])
   B_end = int(indeces[1])
   if (A_start <= s) and (A_end <= s) and (B_start <= s) and (B_end <= s):
    C[i - m][A_start - 1][A_end - 1][j - m][B_start - 1][B_end - 1] = float(line[2]) # upper triangular
#    if (i-m == 15) and (j-m == 25):
#      print(line)
#      print("A_start - 1 = %i" %(A_start - 1))
#      print("A_end - 1 = %i" %(A_end - 1))
#      print("B_start - 1 = %i" %(B_start - 1))
#      print("B_end - 1 = %i" %(B_end - 1))
#      print("{:e}".format(C[i - m][A_start - 1][A_end - 1][j - m][B_start - 1][B_end - 1]))


  #ENERGY COLLECTION
  
  if i == m: #collect the energies and TDM only from the first monomer interaction term
   
   Esa = 0
   Esb = 0
   for q in range(iline + number_coup + 1, len(output)):
    try:
     if output[q].split()[0] == "State":
      if output[q].split()[4] == "1":
       if Esb < s:
        E[j-m][Esb] = float(output[q].split()[5])
        Esb += 1
      elif  (output[q].split()[4] == "2") and (j == m + 1):
       if Esa < s:
        E[i - m][Esa] = float(output[q].split()[5])
        Esa += 1
     elif output[q].split()[0] == "Frenkel":
      iline = q
      break
    except:
     pass
   #end for q
  # TDM COLLECTION: note it is upper triangular
   jline = 0
   for q in range(iline, len(output)):  ###step for TDMs
    try:
#     if output[q].split()[0] == "Dip1:":
     if (output[q].split()[0] == "Transition") and (output[q].split()[5] == "1:"): # Transition dipole vectors of system 1:
#      jline = q + 1
      jline = q
#      print(jline)
      break
    except:
     pass  
#   TDM_s_start = 0
#   TDM_s_end = 1
   #print(jline)
#   for q in range(jline + 1, len(output)):
   for q in range(jline + 2, len(output)):
#    print(output[q])
#    if output[q].split()[0] == "Dip2:":
    if (output[q].split()[0] == "Transition") and (output[q].split()[5] == "2:"): # Transition dipole vectors of system 2:
     jline = q
#     TDM_s_start = 0
#     TDM_s_end = 1
     break
#    if TDM_s_end == TDM_s_start:
#     TDM[j-m][TDM_s_start][TDM_s_end][0] = 0.0
#     TDM[j-m][TDM_s_start][TDM_s_end][1] = 0.0
#     TDM[j-m][TDM_s_start][TDM_s_end][2] = 0.0
#     TDM_s_end += 1
#    elif TDM_s_end != TDM_s_start:
    else:
     TDM_comps = output[q].split() # from    to          x                    y                    z
     TDM_s_start = int(TDM_comps[0]) - 1
     TDM_s_end   = int(TDM_comps[1]) - 1
     x_comp      = float(TDM_comps[2])
     y_comp      = float(TDM_comps[3])
     z_comp      = float(TDM_comps[4])
     TDM[j-m][TDM_s_start][TDM_s_end][0] = x_comp
     TDM[j-m][TDM_s_start][TDM_s_end][1] = y_comp
     TDM[j-m][TDM_s_start][TDM_s_end][2] = z_comp
#     TDM_s_end += 1
#     if TDM_s_end == total_s: 
#      TDM_s_start += 1
#      TDM_s_end = TDM_s_start + 1
   #end for cycle
   if (j == m + 1):
#    for q in range(jline + 1, jline + 1 + TOT_TDM):
    for q in range(jline + 2, jline + 2 + TOT_TDM):
#     if TDM_s_end == TDM_s_start:
#      TDM[i-m][TDM_s_start][TDM_s_end][0] = 0.0
#      TDM[i-m][TDM_s_start][TDM_s_end][1] = 0.0
#      TDM[i-m][TDM_s_start][TDM_s_end][2] = 0.0
#      TDM_s_end += 1
#     elif TDM_s_end != TDM_s_start:
      TDM_comps = output[q].split() # from    to          x                    y                    z
      TDM_s_start = int(TDM_comps[0]) - 1
      TDM_s_end   = int(TDM_comps[1]) - 1
      x_comp      = float(TDM_comps[2])
      y_comp      = float(TDM_comps[3])
      z_comp      = float(TDM_comps[4])
      TDM[i-m][TDM_s_start][TDM_s_end][0] = x_comp
      TDM[i-m][TDM_s_start][TDM_s_end][1] = y_comp
      TDM[i-m][TDM_s_start][TDM_s_end][2] = z_comp
#      TDM_s_end += 1
#      if TDM_s_end == total_s:
#       TDM_s_start += 1
#       TDM_s_end = TDM_s_start + 1


#print(TDM[26][0][1][:])
#TDM[1][0][1][:] = - TDM[1][0][1][:]
#TDM[8][0][1][:] = - TDM[8][0][1][:]
#TDM[9][0][1][:] = - TDM[9][0][1][:]
#TDM[11][0][1][:] = - TDM[11][0][1][:]
#TDM[12][0][1][:] = - TDM[12][0][1][:]
#TDM[13][0][1][:] = - TDM[13][0][1][:]
#TDM[15][0][1][:] = - TDM[15][0][1][:]
#TDM[17][0][1][:] = - TDM[17][0][1][:]
#TDM[19][0][1][:] = - TDM[19][0][1][:]
#TDM[20][0][1][:] = - TDM[20][0][1][:]
#TDM[25][0][1][:] = - TDM[25][0][1][:]
#TDM[26][0][1][:] = - TDM[26][0][1][:]
#print(TDM[26][0][1][:])
#print("{:e}".format(C[15][1][0][25][0][1]))

#print(TDM)
#print("")
#print(TDM[0][0][0][0])
#exit()

#aggregate ground state
Eg0 = 0.0
for i in range(0, A):
 Eg0 = Eg0 + E[i][0]			#sum all the ground states of the monomers
 for j in range(i, A):
  Eg0 = Eg0 + C[i][0][0][j][0][0]   	#sum all the couplings gg,gg between monomer i and all the other j-th,
     				#since C is upper triangular, there is not double counting
 
#######		FRENKEL-DAVYDOV 	#############
#number of aggregate states

N = (A * (s -1)) + 1

#Frankel matrix 
#F[A, state of A, B, state of B]
F = np.zeros(( N, N ))

#DIAGONAL ASSIST-VARIABLE
agg_s = 1
ENA = 0 #the monomer which excites
agg_s_end = 1
CA = 0
diff = 2
a = range(0, A)
p = 0
pp = 0
ppp = 0
w = 1


sss_start = 1
sss_end = 1

for i in range(0, N ):
 for j in range(0, N ):
  if i == j: #DIAGONAL TERMS 
   if i == 0:
    F[i][j] = Eg0
   else:
    F[i][j] = Eg0 - E[ENA][0] + E[ENA][agg_s]
    for M in range(0, A):
     F[i][j] +=  C[ENA][agg_s][agg_s][M][0][0] + C[M][0][0][ENA][agg_s][agg_s] - C[ENA][0][0][M][0][0] - C[M][0][0][ENA][0][0]
     #again we want double sums because the Cmn wants m higher than n
    agg_s += 1
    if agg_s >= s:
     agg_s = 1
     ENA += 1
#    if ENA >= A:
#     print('LOL')
  #end if diagonal
  elif (i == 0) and (j > 0): #FIRST ROW OFF DIAGONAL TERMS
   for M in range(0, A):
    F[i][j] +=  C[CA][0][agg_s_end][M][0][0] + C[M][0][0][CA][0][agg_s_end]
   agg_s_end += 1
   if agg_s_end >= s:
    agg_s_end = 1
    CA += 1
   if CA >= A:
#    print('LOL')
    CA = 0
    agg_s_end = 1
  #end if first row  
  elif (i > 0) and (j > i):
   F[i][j] = C[a[ppp]][sss_start][0][a[w]][0][sss_end]
#   if (a[ppp] == 15) and (a[w] == 25):
#     print(F[i][j]) 
#     print("sss_start = %i" %(sss_start))
#     print("sss_end = %i" %(sss_end))
#     print()
   sss_end += 1
   if sss_end >= s:
    w += 1
    sss_end = 1
    if (w >= A):
     sss_start += 1
     sss_end += 1
     w = ppp + 1
     if sss_start >= s:
      ppp += 1
      sss_start = 1
      sss_end = 1
      w = ppp + 1

#END FRENKEL MATRIX it is in Hartree 

#print(F) 
#print(F*27.211369917461)

#DEFINE THE ENERGY DIFFERENCE F - Eg0 in Hartree

F_del = (F) - (np.diag(np.ones(N)*Eg0))

#print(F_del*27.211369917461)

# DIAGONALISATION of F_del
eigenvalues, eigenvectors = LA.eigh(F_del,UPLO='U')
eigenvalues= eigenvalues - eigenvalues[0] # HERE THE ENERGY IS STILL IN HARTREE, eigenvalues[0] refers to the energy of the ground state that must be zero
					   #They are the exciton transitions from the ground state to the higher states 
#print(eigenvalues*27.211369917461)


#************** EXCITON'S SPECTRAL PROPERTIES SECTION ******************

#************* TRANSITION DIPOLE MOMENT

MON = 0
s_end = 1

#creation of the TDM matrices

TDMx = np.zeros(( N, N ))
TDMy = np.zeros(( N, N ))
TDMz = np.zeros(( N, N ))

#		X AXIS	

for i in range(0, N):
	for j in range(0,N):
		if i == j:
			TDMx[i][j] = float(0)
		elif i == 0:
			TDMx[i][j] = float(TDM[MON][0][s_end][0])
			s_end += 1
			if s_end >= s: 
				MON += 1
				s_end = 1
				if MON >= A:
					MON = 0
		elif j == 0:
			TDMx[i][j] = float(TDM[MON][0][s_end][0])
			s_end += 1
			if s_end >= s:
				MON += 1
				s_end = 1
				if MON >= A:
					MON = 0
#		elif i < s:
#			MON = 0
#			TDMx[i][j] = TDM[MON][i][j][0] 




#               Y AXIS 


for i in range(0, N):
        for j in range(0,N):
                if i == j:
                        TDMy[i][j] = 0.0
                elif i == 0:
                        TDMy[i][j] = float(TDM[MON][0][s_end][1])
                        s_end += 1
                        if s_end >= s:
                                MON += 1
                                s_end = 1
                                if MON >= A:
                                        MON = 0
                                
                elif j == 0:
                        TDMy[i][j] = float(TDM[MON][0][s_end][1])
                        s_end += 1
                        if s_end >= s:
                                MON += 1
                                s_end = 1
                                if MON >= A:
                                        MON = 0

#               z AXIS 


for i in range(0, N):
        for j in range(0,N):
                if i == j:
                        TDMz[i][j] = 0.0
                elif i == 0:
                        TDMz[i][j] = float(TDM[MON][0][s_end][2])
                        s_end += 1
                        if s_end >= s:
                                MON += 1
                                s_end = 1
                                if MON >= A:
                                        MON = 0
                elif j == 0:
                        TDMz[i][j] = float(TDM[MON][0][s_end][2])
                        s_end += 1
                        if s_end >= s:
                                MON += 1
                                s_end = 1
                                if MON >= A:
                                        MON = 0


#print(TDMz)
# at this point we must compute the TDM for the exciton transitions

inv_eigenvectors = LA.inv(eigenvectors)

#	DEFINE THE TDM X, Y, Z FOR EXCITON STATES

TDMx_times_eigvect = np.dot(TDMx, eigenvectors)
TDMy_times_eigvect = np.dot(TDMy, eigenvectors)
TDMz_times_eigvect = np.dot(TDMz, eigenvectors)

TDM_exc_x = np.dot(inv_eigenvectors,TDMx_times_eigvect)
TDM_exc_y = np.dot(inv_eigenvectors,TDMy_times_eigvect)
TDM_exc_z = np.dot(inv_eigenvectors,TDMz_times_eigvect)

# it would be nice to collect the exciton's TDM as done previously 
#	TDM_EXC[starting exciton state][final exciton state][x, y, or z]

TDM_EXC = np.zeros((N, N, 3))
for i in range(0,N): #starting exciton state
 for j in range(i+1, N): #final exciton state
  TDM_EXC[i][j][0] = float(TDM_exc_x[i][j])
  TDM_EXC[i][j][1] = float(TDM_exc_y[i][j])
  TDM_EXC[i][j][2] = float(TDM_exc_z[i][j])


#******************* OSCILLATOR STRENGTH

#define oscillator strength (osc) for each exciton transition
#in that case we limit to S_0 --> S_n transitions

osc = np.zeros(N-1)
S_n = 1
for i in range(0, N-1):
 osc[i] = (2/3)*float(eigenvalues[S_n])*((float(TDM_EXC[0][S_n][0]))**2 + (float(TDM_EXC[0][S_n][1]))**2 +(float(TDM_EXC[0][S_n][2]))**2)  
 S_n += 1


#***********************************************************************************

#***** PRINTOUT SECTION

#numerical names

#outputFileName = "frenkel_#.out"
#outputVersion = 1
#while os.path.isfile(outputFileName.replace("#", str(outputVersion))):
# outputVersion += 1
#outputFileName = outputFileName.replace("#", str(outputVersion))


#according to data

#moment=time.strftime("%d-%b-%H:%M",time.localtime())
moment=time.strftime("%Y-%m-%d--%H-%M",time.localtime())
outputFileName = "frenkel_#.out"
outputFileName = outputFileName.replace("#", moment)

#outputFileName = "frenkel_try.out" # comment if one of the above method are selected

out = open(outputFileName, "w")
out.write("=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=\n\n")
out.write("          FRENKEL-DAVIDOV EXCITON            \n\n")
out.write("=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=\n\n\n")
out.write("=> The Frenkel Matrix computed for your system considering "+str(s)+" states is\n\n")
L = []
for i in range(0,N):
 for j in range(0, N):
  #L += str(round(F_del[i][j]*27.211369917461,4))+'\t' 
  L += "{:.4e}".format(F_del[i][j]*27.211369917461)  + '\t' 
 L += "\n" 
out.writelines(L)
out.write("\n")
out.write("please, notice that the energies are in electronvolts.")
out.write("\n\n")
out.write("=> The diagonalisation process returns the following eigenvalues and eigenvectors \n\n")
out.write("Eigenvalues (eV) \n\n")
L = []
for i in range(0,N):
  L += str(round(eigenvalues[i]*27.211369917461,5)).ljust(10) 
out.writelines(L)
out.write("\n\n")
out.write("Eigenvectors \n\n")
L = []
for i in range(0,N):
 for j in range(0, N):
  L += str(round(eigenvectors[i][j],4)).ljust(10)
 L += "\n"
out.writelines(L)
out.write("\n\n\n")
out.write("=> For each exciton transition we can define the TDM and the oscillator strength \n\n")
L = ''
out.writelines('DeltaE (eV)\tmux(a.u.)\tmuy(a.u.)\tmuz(a.u.)\tOscillator Stregth\n\n')
for i in range(1, N):
 deltaE_eV = str(round(eigenvalues[i]*27.211369917461,5)).ljust(15)
 mux = str(round(TDM_EXC[0][i][0],5)).ljust(10)
 muy = str(round(TDM_EXC[0][i][1],5)).ljust(10)
 muz = str(round(TDM_EXC[0][i][2],5)).ljust(10)
 f = str(round(osc[i-1],5))
 L += deltaE_eV + '\t' + mux + '      ' + muy + '      ' + muz + '      ' + f + '\n'
out.writelines(L)
out.write("\n\n\n")





out.close()