#!/usr/bin/python
# -*- coding: utf-8 -*-

# Projected Site Occupation Embedding Theory
# This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import absolute_import

import os
import sys
import subprocess
import approximate_functionals as app_func
import numpy as np
from scipy.linalg import eigh
from scipy import optimize


def generate_Hamiltonian(L,N,t,v):
      """
      Function that generates the 
      non-interacting Hamiltonian 
      with potential v
      """
      H = np.zeros((L,L),dtype=float)
      for i in range(L-1):
       H[i,i+1] = H[i+1,i] = +t
       H[i,i] = v[i]
       H[L-1,0] = H[0,L-1] = +t
      if ((N/2) % 2 == 1): #periodic condition
       H[L-1,0] = H[0,L-1] = +t
      elif ((N/2) % 2 == 0): #antiperiodic condition
       H[L-1,0] = H[0,L-1] = -t
      H[L-1,L-1] = v[L-1]
      return H

def solve_Hamiltonian(L,N,H):
      E,C = eigh(H)
      # Fundamental energy
      E0 = E[0]
      # 1RDM
      onepdm = np.zeros((L,L),dtype=float)
      occ = [0]*L
      for i in range(L):
       for j in range(L):
        for k in range(int(N/2)):
         onepdm[i,j]+=C[i,k]*C[j,k]
       occ[i] = 2*onepdm[i,i] # Because the onepdm is for a given spin.
      return onepdm,occ

def write_mat(L,N,t,U,n_imp,mat):
    with open("mat.dat",'w') as f:
       f.write(str(L) + " " + str(N) + " " + str(t) + " " + str(U) + " " + str(n_imp) + "\n")
       for row in mat:
         for val in row:
           f.write('%20.15f' % val)
         f.write("\n")

def generate_potential(L,U,t,occ,beta,dbeta_dU,n_imp,approx):
      """
      Function that generates the embedding potential as well as
      the derivative of the Hxc per-site energy with respect to U
      """
      deHxc_dn = [0]*L
      dEHxcimp_dn = [0]*L
      dEHxcbath_dn = [0]*L
      for i in range(n_imp):
        if approx == "iBALDA":
          dEHxcimp_dn[i] = app_func.correlation_iBALDA(U,t,occ[i],beta,dbeta_dU)[3] + U*occ[i]*0.5
        elif approx == "2LBALDA":
          dEHxcimp_dn[i] = app_func.correlation_2LBALDA(U,t,occ[i])[3] + U*occ[i]*0.5
        elif approx == "SIAMBALDA":
          dEHxcimp_dn[i] = app_func.correlation_SIAMBALDA(U,t,occ[i])[3] + U*occ[i]*0.5
      for i in range(L):
        deHxc_dn[i] = app_func.correlation_BALDA(U,t,occ[i],beta,dbeta_dU)[3] + U*occ[i]*0.5
        dEHxcbath_dn[i] = deHxc_dn[i] - dEHxcimp_dn[i]
      return dEHxcbath_dn

def solve_dimer(U,t,deltav):
    occ = [0]
    dimp = [0]
    u = U/(2*t)
    nu = deltav/(2*t)
    w = np.sqrt(3.0*(1 + nu**2) + u**2)
    theta = (1/3.0)*np.arccos((u/(w**3))*(9*(nu**2 - 0.5) - u**2))
    Edimer = (4*t/3.0)*(u - w*np.sin(np.pi/6.0 + theta))
    dE_dt = 8*t*(U - Edimer)/(-3*Edimer**2 + 4*U*Edimer + (2*t)**2 - U**2 + deltav**2)
    dE_dU = ((2*t)**2 + 2*Edimer*(U - Edimer))/(-3*Edimer**2 + 4*U*Edimer + (2*t)**2 - U**2 + deltav**2)
    dE_ddeltav = (-2*deltav*Edimer)/(-3*Edimer**2 + 4*U*Edimer + (2*t)**2 - U**2 + deltav**2)
    gamma_12 = -0.5*dE_dt
    occ = [1 - dE_ddeltav]
    dimp = [0.5*dE_dU + (1 - occ[0])*(-1/2.0)]
    return occ, dimp

def persite_and_dblocc(occ,dimp,U,t,n_imp,beta,dbeta_dU,approx):
    persite = 0
    dblocc = 0
    for i in range(n_imp):
      (ec,dec_dU,dec_dt,dec_dn) = app_func.correlation_BALDA(U,t,occ[i],beta,dbeta_dU)
      if approx == "iBALDA":
        (Ecimp,dEcimp_dU,dEcimp_dt,dEcimp_dn) = app_func.correlation_iBALDA(U,t,occ[i],beta,dbeta_dU)
      elif approx == "2LBALDA":
        (Ecimp,dEcimp_dU,dEcimp_dt,dEcimp_dn) = app_func.correlation_2LBALDA(U,t,occ[i])
      elif approx == "SIAMBALDA":
        (Ecimp,dEcimp_dU,dEcimp_dt,dEcimp_dn) = app_func.correlation_SIAMBALDA(U,t,occ[i])
      decbath_dU = dec_dU - dEcimp_dU
      dblocc = dblocc + (dimp[i] + decbath_dU)/n_imp
      persite = persite + (-4*t*np.sin(np.pi*occ[i]/2.0)/np.pi + t*dec_dt + U*dimp[i] + U*decbath_dU)/n_imp
    return persite, dblocc

def write_FCIDUMP(U,h_emb,n_imp):
   with open("FCIDUMP","w") as f:
      f.write("&FCI NORB= {}, NELEC= {}, MS2= 0,\n".format(2*n_imp,2*n_imp))
      f.write(" ORBSYM="+"1,"*(2*n_imp)+"\n")
      f.write("ISYM=1,\n")
      f.write("&END\n")
      for i in range(n_imp):
         f.write('%15.10f %3d %3d %3d %3d\n' % (U,i+1,i+1,i+1,i+1))
      for i in range(n_imp,2*n_imp):
         f.write('%15.10f %3d %3d %3d %3d\n' % (0,i+1,i+1,i+1,i+1))
      for i in range(2*n_imp-1):
         for j in range(i+1,2*n_imp):
            f.write('%15.10f %3d %3d %3d %3d\n' % (h_emb[i,j],i+1,j+1,0,0))
      for i in range(2*n_imp):
         f.write('%15.10f %3d %3d %3d %3d\n' % (h_emb[i,i],i+1,i+1,0,0))

def write_dmrg_conf(n_imp,m,opt):
   with open("dmrg_"+str(opt)+".conf","w") as f:
      f.write("sym c1\n")
      f.write("orbitals FCIDUMP\n")
      f.write("nelec {}\n".format(2*n_imp))
      f.write("spin 0\n")
      f.write("irrep 1\n")
      f.write("hf_occ integral\n")
      f.write("schedule default\n")
      f.write("maxM {}\n".format(m))
      f.write("maxiter 200\n")
      f.write(str(opt))
  
def read_dmrg_onepdm_twopdm(n_imp):
    with open("spatial_onepdm.0.0.txt") as f:
     f.readline()
     onepdm_tmp = [[token for token in line.split()] for line in f.readlines()]
    with open("spatial_twopdm.0.0.txt") as f:
     f.readline()
     twopdm_tmp = [[token for token in line.split()] for line in f.readlines()]
    onepdm={}
    twopdm={}
    for row in range(len(onepdm_tmp)):
       onepdm[int(onepdm_tmp[row][0]),int(onepdm_tmp[row][1])] = float(onepdm_tmp[row][2])
    for row in range(len(twopdm_tmp)):
       twopdm[int(twopdm_tmp[row][0]),int(twopdm_tmp[row][1]),int(twopdm_tmp[row][2]),int(twopdm_tmp[row][3])] = float(twopdm_tmp[row][4])
    occ = [0]*n_imp
    dimp = [0]*n_imp
    for i in range(n_imp):
       occ[i] = onepdm[i,i]
       dimp[i] = twopdm[i,i,i,i]
    return occ,dimp

def self_consistent_loop(L,N,U,t,m,occ,beta,dbeta_dU,n_imp,approx,P,chem_pot,code_directory):
   occ_exact = [N/(1.0*L)]*L
   old_occ = [0]*L
   delta_occ = 0
   for i in range(L):
      old_occ[i] = occ[i]
   print("""   Step 1/5: Determine the one-body effective Hamiltonian.
            This is the SOET one without interaction 
            (or, equivalently, the noninteracting one with the embedding potential).\n
            You selected the following approximation : {}, with {} impurities\n""".format(approx,n_imp))
   dEHxcbath_dn = generate_potential(L,U,t,occ,beta,dbeta_dU,n_imp,approx)
   h_eff = np.asarray(generate_Hamiltonian(L,N,t,dEHxcbath_dn))

   print("   Step 2/5: Project the one-body effective Hamiltonian to create the one-body embedded Hamiltonian.\n")
   try:
      h_emb = np.dot(np.dot(P.T,h_eff),P)
   except:
      print("Failed to project the effective one-body Hamiltonian.")
      sys.exit(0)

   print("   Step 3/5: Add the coulomb interaction on the impurity sites and solve the full Embedded problem.\n")
   if n_imp == 1:
      deltav = h_emb[1,1] - h_emb[0,0]
      print("  Only one impurity is considered. The embedded problem results to a 2 electrons asymmetric Hubbard dimer.")
      print("  It is solved analytically with the following parameters:")
      print("      Coulomb repulsion                               : {}".format(U))
      print("      Projection of the hopping, t --> -h_emb[0,1]    : {}".format(-h_emb[0,1]))
      print("      Delta v = v_1 - v_0 = h_emb[1,1] - h_emb[0,0]   : {}\n".format(deltav))
      print("  The impurity problem is related to the fully-interacting one by Eimp(U,deltav) = E(U/2,deltav - U/2).\n")

      try:
         if chem_pot:
           mu = optimize.fminbound(lambda x: abs(solve_dimer(U/2.0,-h_emb[0,1],deltav - U/2.0 + x)[0][0] - N/(1.0*L)), -1000,1000)
           occ_imp, dimp = solve_dimer(U/2.0,-h_emb[0,1], deltav - U/2.0 + mu)
           print("Optimization of the chemical potential : mu = {}\n".format(mu))
         else:
           occ_imp, dimp = solve_dimer(U/2.0,-h_emb[0,1],deltav - U/2.0)
      except:
         print("Failed to solve the embedded Hubbard dimer.")
         sys.exit(0)

   else:
      print("  {} impurities are considered. The embedded problem has {} sites with {} interacting sites.".format(n_imp,2*n_imp,n_imp))
      print("  It is solved numerically by Density Matrix Renormalization Group (https://github.com/sanshar/Block).\n") 
      write_FCIDUMP(U,h_emb,n_imp)
      write_dmrg_conf(n_imp,m,"onepdm")
      write_dmrg_conf(n_imp,m,"twopdm")
      try:
         subprocess.check_call(code_directory+"dmrg_script.sh",shell=True)
      except:
         print("Could not run the dmrg_script.sh, check if it is in code_directory.")
         sys.exit(0)
      occ_imp,dimp = read_dmrg_onepdm_twopdm(n_imp)

   print("""   Step 4/5: Compute the SOET per-site energy and double occupation.
            Use the values of the occupation and the double occupation of the impurity
            together with the formula of SOET.\n""")
   persite, dblocc = persite_and_dblocc(occ_imp,dimp,U,t,n_imp,beta,dbeta_dU,approx)
   persite_nexact, dblocc_nexact = persite_and_dblocc(occ_exact,dimp,U,t,n_imp,beta,dbeta_dU,approx)

   print("   Step 5/5: Set all the occupation in the bath equal to the average of the impurity occupations.\n") 
   occ_average = 0
   for i in range(n_imp):
      occ[i] = occ_imp[i]
      delta_occ = delta_occ + abs(old_occ[i] - occ_imp[i])/n_imp
      occ_average = occ_average + occ_imp[i]/n_imp
   for i in range(n_imp,L):
      occ[i] = occ_average
   print("   ------------  RESULTS  ------------\n")
   print("   Occupation of the impurity site(s):")
   for i in range(n_imp):
      print("     site {}: {}".format(i,occ[i]))
   print("")
   print("   Per-site energy   : {}".format(persite))
   print("   Double occupation: {}\n".format(dblocc))
   print("   And using the exact uniform density:")
   print("   Per-site energy  : {}".format(persite_nexact))
   print("   Double occupation: {}\n".format(dblocc_nexact))
   print("   Delta occ = {}\n".format(delta_occ))
   return occ, delta_occ, persite, dblocc, persite_nexact, dblocc_nexact 




def run_psoet(L,N,U,t,n_imp,approx,code_directory,
              single_shot = False,
              chem_pot = False,
              m = 1000,
              MAXITER = 50,
              threshold = 0.0001):
   """
   Check if the initial values are valid:
   """
   if approx != "iBALDA" and approx != "2LBALDA" and approx != "SIAMBALDA":
     raise ValueError("approx must be either iBALDA, 2LBALDA or SIAMBALDA approximation.")
   if approx == "2LBALDA" or approx == "SIAMBALDA":
     if n_imp != 1:
       raise ValueError("The approximation is valid for a single impurity only.")
   if not isinstance(L,int):
     raise TypeError("The number of sites should be an integer.")
   if not isinstance(N,int):
     raise TypeError("The number of electrons should be an integer.")
   if not isinstance(n_imp,int):
     raise TypeError("The number of impurites should be an integer.")
   if (N % 2 != 0):
     raise ValueError("The number of electrons should be an even integer.")
   if N < 2*n_imp:
     raise ValueError("The number of electrons must be higher than twice the number of impurities.")
   if N > 2*L:
     raise ValueError("The number of electrons must be lower than twice the number of sites.")
   if n_imp > L/2.0:
     raise ValueError("The number of impurities must be lower than half the number of sites.\nIf n_imp = n_sites/2, the embedding is exact.")
   if chem_pot and n_imp > 1:
     raise ValueError("The optimization with a global chemical potential is not set for multiple impurities yet.")


   print("Starting PSOET for {} sites, {} eletrons, U/t = {} and {} impurity site(s)\n".format(L,N,U/t,n_imp))

   print("#"*60)
   print("            INITIAL STEPS BEFORE SELF-CONSISTENCE")
   print("#"*60 + "\n")
   print("""   Step 1/2: Perform a KS-SOFT calculation.
            Given that the system is uniform and that the potential is defined up to a constant,
            we set it to 0 and the calculation is done in only one iteration.
            This calculation gives us the 1RDM from which we will apply the Schmidt Decomposition 
            to obtain a new Schmidt basis (new bath states) and a Projector.\n""")
   try: 
      h_SOFT = generate_Hamiltonian(L,N,t,[0]*L)
      onepdm_SOFT, occ = solve_Hamiltonian(L,N,h_SOFT)
   except:
      print("Failed to do the KS-SOFT calculation.")
      sys.exit(0)

   print("""   Step 2/2: Perform the Schmidt decomposition of the KS-SOFT 1RDM to obtain the projector P.
            This projector will be used to project the effective Hamiltonian onto the Schmidt basis,
            corresponding to the basis of the new embedding problem.\n""")
   try:
     write_mat(L,N,U,t,n_imp,onepdm_SOFT)
     subprocess.check_call(code_directory+"schmidt_decomposition", shell=True)
   except:
     print("Failed to generate the projector P. Look at mat.dat and schmidt_decomposition.f90.")
     sys.exit(0)
   with open("projector.dat") as f:
     P = np.asarray([[float(token) for token in line.split()] for line in f.readlines()])


   # Compute beta and dbeta_dU for the approximate functionals.
   try:
      subprocess.check_call("echo " + str(U) + " " + str(t) + " | " + code_directory + "beta_and_derivatives",shell=True)
      with open("beta_dbetadU.dat","r") as f:
         line = f.read()
         beta = float(line.split()[0])
         dbeta_dU = float(line.split()[1])
   except:
      print("Failed to compute beta. Look at beta_and_derivatives.f90 and its compilation and execution.")
      sys.exit(0)


   print("#"*60)
   print("                     START ITERATIONS")
   print("#"*60 + "\n")
   iteration = 1
   delta_occ = 1
   if N/L == 1: # This option is necessary as otherwise, the occupation obtained from the KS calculation, which is subject to numerical error, oscillates around 1 and leads to convergence problem. This is known in the literature, and it is due to the Mott--Hubbard transition at half-filling (--> discontinuity in the potential).
      occ = [1.0*N/L]*L
   if single_shot is True:
      MAXITER = 1
   if chem_pot is True:
      name = "L{}_N{}_U{}_t{}_nimp{}_{}_chempot.out".format(str(L),str(N),str(U),str(t),str(n_imp),approx)
   else:
      name = "L{}_N{}_U{}_t{}_nimp{}_{}.out".format(str(L),str(N),str(U),str(t),str(n_imp),approx)
   with open(name,"w") as f:
      f.write("--------------------------------- SUMMARY OF THE RESULTS ---------------------------------\n")
      f.write("Exact uniform density : {}\n".format(1.0*N/L))
      f.write('%8s ' % ("Iter"))
      for i in range(n_imp):
         f.write('%12s ' % ("occ_" + str(i)))
      f.write('%12s %15s %15s %15s %15s\n' % ("Delta_occ","Per-site","P-s_nexact","Dbl-occ","D-o_nexact"))
      while delta_occ > threshold and iteration <= MAXITER:
         print("#"*30)
         print("        ITERATION " + str(iteration))    
         print("#"*30+"\n")
         occ, delta_occ, persite, dblocc, persite_nexact, dblocc_nexact = self_consistent_loop(L,N,U,t,m,occ,beta,dbeta_dU,n_imp,approx,P,chem_pot,code_directory)
         f.write('%8d ' % (iteration))
         for i in range(n_imp):
            f.write('%12.8f ' % occ[i])
         f.write('%12.8f %15.8f %15.8f %15.8f %15.8f\n' % (delta_occ,persite,persite_nexact,dblocc,dblocc_nexact))
         iteration += 1
