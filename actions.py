#!/usr/bin/env python

from scipy.special import comb
from numpy import linalg as LA
import numpy as np
import re
import os

def collect(first_monomer: int, 
            last_monomer: int, 
            computed_roots: int, 
            wanted_roots: int, 
            source = os.getcwd(),
            key = ['Excitonic', 'couplings', 'section:']
            ):
    """
    Collect the results of Quantum Chemistry calculations from the OpenMolcas Software. It is specific for the building of Frenkel-Davydov Hamiltonian. 
    All the coupling outputs must be in the same directory and they have name "'coup_'+alpha+'_'+beta" where alpha and beta are the monomer indices. 
    The output files must have the name "'coupling_'+alfa+'_'+beta+'.out'".

    Parameters
    ----------
    first_monomer: int
        lowest index of the monomers

    last_monomer: int
        highest index of the monomers

    computed_roots: int
        number of computed roots in the CASSF calculation per monomer

    wanted_roots: int
        number of roots to include in the calculation per monomer

    Returns
    -------
    Returns a tuple with the following terms:

    C : np.ndarray
        triangular tensor of coupling terms with form:
        C = [Mol_i][Mol_i_state_start][Mol_i_state_end][Mol_j][Mol_j_state_start][Mol_j_state_end]
        Mol_i goes from 0 to the number of molecules.

    E : np.ndarray
        matrix of monomeric energies, with form:
        E = [Mol_i][Mol_i_state]
        Mol_i goes from 0 to the number of molecules.

    TDM : np.ndarray
        tensor of transition dipole moments, with form:
        TDM = [Mol_i][Mol_i_state_start][Mol_i_state_end][space components (x,y,z)]
        Mol_i goes from 0 to the number of molecules.

    """

    
    m = first_monomer #first monomer
    n = last_monomer #last monomer
    total_s = computed_roots #number of computed roots
    s = wanted_roots # number of states considered in FRENKEL e.g. 2 states for Qy band

    A = n - m + 1 # number of molecules
    number_coup = int(comb(total_s + total_s, 2))


    """
    variable that collects the coupling terms. Ordering:
      C = [Mol_i][Mol_i_state_start][Mol_i_state_end][Mol_j][Mol_j_state_start][Mol_j_state_end]
    Mol_i goes from 0 to A-1
    """
    C = np.zeros((A, s, s, A, s, s))

    """
    variable that collects the energy terms. Ordering:
        E = [Mol_i][Mol_i_state]
    Mol_i goes from 0 to A-1
    """
    E = np.zeros((A, s))

    """
    variable that collects the transition dipole moments (TDM) to compute spectroscopic properties. Ordering:
      TDM = [Mol_i][Mol_i_state_start][Mol_i_state_end][space components (x,y,z)]
    """
    TDM = np.zeros((A,total_s, total_s, 3))
    TOT_TDM = int(comb(total_s, 2))

    # The cycle starts here
    for i in range(m, n+1): 
        alfa = str(i)
        for j in range(i+1, n+1):
            beta = str(j)
            #define the destination directory for the files
            dest = os.path.join(source,'coup_'+alfa+'_'+beta)
            #define the path for the output file
            out = os.path.join(dest, 'coupling_'+alfa+'_'+beta+'.out')
            #go to the parsing keyword
            with open(out, "r") as F:
                output = F.readlines()


            #COUPLING TERMS COLLECTION

            #find the line number of the keyword 
            iline = 0
            for k in range(len(output),0, -1): # searches from the end of the output file
                try:
                    if output[k].split() == key:
                        iline = k
                        break
                except:
                    pass
            iline += 4 # adds four lines before the start of the couplings; only for new versions of Openmolcas
            for l in range(1, number_coup + 1):
                line = output[iline + l].split()
                indeces = re.split(r'\(|\)', line[0])
                A_start = int(indeces[7])
                A_end = int(indeces[3])
                B_start = int(indeces[5])
                B_end = int(indeces[1])
                if (A_start <= s) and (A_end <= s) and (B_start <= s) and (B_end <= s):
                    C[i - m][A_start - 1][A_end - 1][j - m][B_start - 1][B_end - 1] = float(line[2]) # upper triangular

            if i == m: #collect the energies and TDM only from the first monomer interaction term

                #ENERGY COLLECTION
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

                # TDM COLLECTION: note it is upper triangular
                jline = 0
                for q in range(iline, len(output)):  ###step for TDMs
                    try:
                        if (output[q].split()[0] == "Transition") and (output[q].split()[5] == "1:"): # Transition dipole vectors of system 1:
                            jline = q
                            break
                    except:
                        pass
                for q in range(jline + 2, len(output)):
                    if (output[q].split()[0] == "Transition") and (output[q].split()[5] == "2:"): # Transition dipole vectors of system 2:
                        jline = q
                        break
                    else:
                        TDM_comps = output[q].split()
                        TDM_s_start = int(TDM_comps[0]) - 1
                        TDM_s_end   = int(TDM_comps[1]) - 1
                        x_comp      = float(TDM_comps[2])
                        y_comp      = float(TDM_comps[3])
                        z_comp      = float(TDM_comps[4])
                        TDM[j-m][TDM_s_start][TDM_s_end][0] = x_comp
                        TDM[j-m][TDM_s_start][TDM_s_end][1] = y_comp
                        TDM[j-m][TDM_s_start][TDM_s_end][2] = z_comp
                if (j == m + 1):
                    for q in range(jline + 2, jline + 2 + TOT_TDM):
                        TDM_comps = output[q].split()
                        TDM_s_start = int(TDM_comps[0]) - 1
                        TDM_s_end   = int(TDM_comps[1]) - 1
                        x_comp      = float(TDM_comps[2])
                        y_comp      = float(TDM_comps[3])
                        z_comp      = float(TDM_comps[4])
                        TDM[i-m][TDM_s_start][TDM_s_end][0] = x_comp
                        TDM[i-m][TDM_s_start][TDM_s_end][1] = y_comp
                        TDM[i-m][TDM_s_start][TDM_s_end][2] = z_comp

    return C, E, TDM
 
def build(E: np.ndarray,
          C: np.ndarray,
          units: str = 'Hartree'
          ):
    """
    Build a Frenkel-Davydov Hamiltonian from the monomeric energies and coupling terms.

    Parameters
    ----------

    E : np.ndarray
        matrix of monomeric energies, with form:
        E = [Mol_i][Mol_i_state]
        Mol_i goes from 0 to the number of molecules.

    C : np.ndarray
        triangular tensor of coupling terms with form:
        C = [Mol_i][Mol_i_state_start][Mol_i_state_end][Mol_j][Mol_j_state_start][Mol_j_state_end]
        Mol_i goes from 0 to the number of molecules.
    
    units: str = 'Hartree'
        specify the unit of the energies. It can be 'Hartree' or 'eV'.

    Returns
    -------
    Returns a a matrix:

    F_Ha : np.ndarray
        Frenkel-Davydov Hamiltonian matrix in Hartree or eV.

    """
        
    A = E.shape[0] #number of monomers
    s = E.shape[1] # number of states considered in FRENKEL
    #define aggregate ground state
    Eg0 = sum(E[i][0] for i in range(A))
    Eg0 += sum(C[i][0][0][j][0][0] for i in range(A) for j in range(i, A))

    def compute_diagonal_energy(ENA, agg_s):
        energy = Eg0 - E[ENA][0] + E[ENA][agg_s]
        for M in range(A):
            energy += (
                C[ENA][agg_s][agg_s][M][0][0] + C[M][0][0][ENA][agg_s][agg_s]
                - C[ENA][0][0][M][0][0] - C[M][0][0][ENA][0][0]
            )
        return energy
    
    N = (A * (s -1)) + 1
    F_Ha = np.zeros(( N, N ))

    #DIAGONAL ASSIST-VARIABLE
    agg_s = 1
    ENA = 0 	#this is a variable that counts the excited monomer in the diagonal terms
    agg_s_end = 1
    CA = 0		#
    p = 0
    w = 1
    sss_start = 1
    sss_end = 1

    for i in range(0, N ):
        for j in range(0, N ):
            #Diagonal terms
            if i == j:
                if i == 0:
                    F_Ha[i][j] = Eg0
                else:
                    F_Ha[i][j] = compute_diagonal_energy(ENA, agg_s)
                    agg_s += 1 # each time i == j
                    if agg_s >= s:
                        agg_s = 1
                        ENA += 1
            #FIRST ROW OFF-DIAGONAL TERMS
            elif (i == 0) and (j > 0):
                F_Ha[i][j]  = sum(
                    C[CA][0][agg_s_end][M][0][0] + C[M][0][0][CA][0][agg_s_end] for M in range(A)
                    )
                agg_s_end += 1
                if agg_s_end >= s:
                    agg_s_end = 1
                    CA += 1
                    if CA >= A:
                        CA = 0
                        agg_s_end = 1
            #other OFF-DIAGONAL TERMS
            elif (i > 0) and (j > i):
                p = (i - 1) // (s - 1)
                sss_start = (i - 1) % (s - 1) + 1
                w = (j - 1) // (s - 1)
                sss_end = (j - 1) % (s - 1) + 1
                if p == w: # Intramolecular coupling within the same molecule
                    F_Ha[i][j] = (sum(C[p][sss_start][sss_end][i][0][0] for i in range(A)) 
                            +  sum(C[i][0][0][p][sss_start][sss_end] for i in range(A))
                            )
                else:  # Intermolecular coupling between different molecules
                    F_Ha[i][j] = C[p][sss_start][0][w][0][sss_end]
    
    if units == 'eV':
        F_Ha = F_Ha * 27.2114
    
    #create the hermitian matrix
    F_Ha = F_Ha - np.diag(np.diag(F_Ha)) + np.transpose(np.conj(F_Ha)) 
    F_Ha_del = (F_Ha) - (np.diag(np.ones(N) * np.min(F_Ha)))

    return F_Ha_del

def compute(Frenkel: np.ndarray,
            TDM: np.ndarray,
            printed_output: bool = False
            ):
    """
    Compute the physical properties of an excitonic system.

    Parameters
    ----------
    F : np.ndarray
        Frenkel-Davydov Hamiltonian matrix in Hartree or eV.
    
    TDM : np.ndarray
        tensor of transition dipole moments, with form:
        TDM = [Mol_i][Mol_i_state_start][Mol_i_state_end][space components (x,y,z)]
        Mol_i goes from 0 to the number of molecules.

    Returns
    -------
    Returns a tuple with the following terms:

    eigenvalues : np.ndarray
        array of eigenvalues of the Frenkel-Davydov Hamiltonian. The unit is the same as the input Hamiltonian.

    eigenvectors : np.ndarray
        array of eigenvectors of the Frenkel-Davydov Hamiltonian.
    
    f : np.ndarray
        array of oscillator strengths of the excitonic transitions only from the ground state to the excited states.
    
    TDM_EXC : np.ndarray
        tensor of transition dipole moments of the excitonic transitions, with form:
        TDM_EXC = np.zeros((N, N, 3))
        N is the number of states considered in the FRENKEL Hamiltonian.
    """

    
    A = TDM.shape[0]
    N = Frenkel.shape[0]
    s = int(1 + ((N - 1)/A))
    eigenvalues, eigenvectors = LA.eigh(Frenkel,UPLO='U')
    eigenvalues= eigenvalues - eigenvalues[0]
    
    MON = 0
    s_end = 1

    #creation of the TDM matrices

    TDMx = np.zeros(( N, N ))
    TDMy = np.zeros(( N, N ))
    TDMz = np.zeros(( N, N ))

    # X AXIS
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
    
    # Y AXIS
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
    # Z AXIS
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

    inv_eigenvectors = LA.inv(eigenvectors)

    TDM_exc_x = inv_eigenvectors @ TDMx @ eigenvectors
    TDM_exc_y = inv_eigenvectors @ TDMy @ eigenvectors
    TDM_exc_z = inv_eigenvectors @ TDMz @ eigenvectors

    TDM_EXC = np.zeros((N, N, 3))
    for i in range(0,N): #starting exciton state
        for j in range(i+1, N): #final exciton state
            TDM_EXC[i][j][0] = float(TDM_exc_x[i][j])
            TDM_EXC[i][j][1] = float(TDM_exc_y[i][j])
            TDM_EXC[i][j][2] = float(TDM_exc_z[i][j])

    f = np.zeros(N-1) #oscillator strength
    S_n = 1
    for i in range(0, N-1):
        f[i] = (2/3)*float(eigenvalues[S_n]) * ((float(TDM_EXC[0][S_n][0]))**2 + (float(TDM_EXC[0][S_n][1]))**2 +(float(TDM_EXC[0][S_n][2]))**2)  
        S_n += 1

    if printed_output == True:
        import time
        moment = time.strftime("%d_%b_%H:%M",time.localtime())
        outputFileName = "frenkel_#.out"
        outputFileName = outputFileName.replace("#", moment)

        with open(outputFileName, "w") as file:
            file.write("=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=\n\n")
            file.write("          FRENKEL-DAVIDOV EXCITON            \n\n")
            file.write("=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=\n\n\n")
            file.write("=> The Frenkel Matrix computed for your system considering "+str(s)+" states is\n\n")
            L = []
            for i in range(0,N):
                for j in range(0, N):
                    #L.append("{:.4e}".format(Frenkel[i][j])  + '\t')
                    L.append(f"{Frenkel[i][j]:<14.4e}")
                L .append("\n") 
            file.writelines(L)
            file.write("\n")
            file.write("please, notice that the energies have same units of the starting Frenkel matrix.")
            file.write("\n")
            file.write("=> The diagonalisation process returns the following eigenvalues and eigenvectors \n\n")
            file.write("Eigenvalues: \n")
            L = []
            for i in range(0,N):
                #L.append(str(round(eigenvalues[i],4)).ljust(10) )
                L.append(f"{eigenvalues[i]:<14.4e}")
            file.writelines(L)
            file.write("\n\n")
            file.write("Eigenvectors: \n")
            L = []
            for i in range(0,N):
                for j in range(0, N):
                    #L.append(str(round(eigenvectors[i][j],4)).ljust(10))
                    L.append(f"{eigenvectors[i][j]:<14.4e}")
                L.append("\n")
            file.writelines(L)
            file.write("\n\n\n")
            file.write("=> For each exciton transition we can define the TDM and the oscillator strength \n\n")
            file.write(f"{'DeltaE':<15}{'mux (a.u.)':<12}{'muy (a.u.)':<12}{'muz (a.u.)':<12}{'Oscillator Strength':<20}\n\n")
            L = []
            for i in range(1, N):
                deltaE_eV = f"{round(eigenvalues[i], 3):<15}"
                mux = f"{round(TDM_EXC[0][i][0], 5):<12}"
                muy = f"{round(TDM_EXC[0][i][1], 5):<12}"
                muz = f"{round(TDM_EXC[0][i][2], 5):<12}"
                osc = f"{round(f[i - 1], 5):<20}"
                # Append the formatted row to the lines list
                L.append(f"{deltaE_eV}{mux}{muy}{muz}{osc}\n")
            file.writelines(L)
            # Write the final newlines
            file.write("\n")
    
    return eigenvalues, eigenvectors, f, TDM_EXC

def get_spectrum(E: np.array, 
                 f: np.array,
                 sigma: float = 0.0272, #the sigma parameter must be adjusted according to the unit of measurement
                 ):
    """
    Compute a spectrum by applying a Gaussian function to the excitation energies.

    Parameters
    ----------

    E : np.ndarray
        array of excitation energies.
    
    f : np.ndarray
        array of oscillator strengths of the excitonic transitions only from the ground state to the excited states.
    
    sigma : float = 0.0272
        standard deviation of the Gaussian function.

    Returns
    -------
    Returns a tuple with the following terms:

    xaxis : np.ndarray
        x axis of the spectrum in units of energy. The unit is the same as the input energies.
    
    spectrum_out : np.ndarray
        y axis of the spectrum normalised to the maximum value. It is in arbitrary units.

    """
    if E[0] == 0:
        E = E[1:]
    
    E_min = np.min(E)
    E_max = np.max(E)
    xaxis = np.linspace(E_min-(3.5*sigma), 
                        E_max+(3.5*sigma), 1000)
    spectrum = 0
    for i in range(0,len(E)):
        spectrum += (f[i]/sigma) * np.exp(-((xaxis-E[i])**2)/(sigma**2))
    spectrum_out = spectrum/np.max(spectrum)

    return xaxis, spectrum_out
            










