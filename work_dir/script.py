#!/usr/bin/python 
import sys
code_directory = "/b/home/quant/bsenjean/P-SOET/code/"
sys.path.append(code_directory)
import run_PSOET

n_site = 12
n_elec = 6
U      = 5.0
t      = 1.0
n_imp  = 1
# approx : "iBALDA, "2LBALDA", "SIAMBALDA"
approx = "2LBALDA"

#optional arguments:
# m (--> number of renormalized states in DMRG calculation. Default : 1000)
# single_shot (True --> MAXITER = 1. Default = False)
# MAXITER (--> set the maximum number of iterations. Default : 50)
# threshold (--> set the threshold for the convergence of the impurity occupation(s). Default : 0.0001)
# chem_pot (True --> Optimization with the global chemical potential to recover the exact density. Implemented for a single impurity for now.)

m = 1000
threshold = 0.0001
chem_pot = True
run_PSOET.run_psoet(n_site,n_elec,U,t,n_imp,approx,code_directory,
                    m = m,
                    chem_pot = chem_pot,
                    threshold = threshold)
