#!/usr/bin/python 
import sys
code_directory = "/b/home/quant/bsenjean/P-SOET/code/"
sys.path.append(code_directory)
import run_PSOET

n_site = 400
t      = 1.0
# approx : "iBALDA, "2LBALDA", "SIAMBALDA"
m = 1000
threshold = 0.0001
opt_pot = True

for approx in ["iBALDA","2LBALDA"]:
 for U in [4.0,8.0]:
  for n_imp in [1]:
   for n_elec in range(2*n_imp,402,2):

#optional arguments:
# m (--> number of renormalized states in DMRG calculation. Default : 1000)
# single_shot (True --> MAXITER = 1. Default = False)
# MAXITER (--> set the maximum number of iterations. Default : 50)
# threshold (--> set the threshold for the convergence of the impurity occupation(s). Default : 0.0001)
# opt_pot (True --> Optimization with the global chemical potential to recover the exact density. Implemented for a single impurity for now.)

      run_PSOET.run_psoet(n_site,n_elec,U,t,n_imp,approx,code_directory,
                    m = m,
                    opt_pot = opt_pot,
                    threshold = threshold)
