#!/usr/bin/python 
import sys
code_directory = "/b/home/quant/bsenjean/P-SOET_final/code/"
sys.path.append(code_directory)
import run_PSOET

n_site = 12
n_elec = 6
U      = 5.0
t      = 1.0
n_imp  = 2
# approx : "iBALDA, "2LBALDA", "SIAMBALDA"
approx = "iBALDA"

#optional arguments:
# m (--> number of renormalized states in DMRG calculation. Default : 1000)
# single_shot (True --> MAXITER = 1. Default = False)
# MAXITER (--> set the maximum number of iterations. Default : 50)
# threshold (--> set the threshold for the convergence of the impurity occupation(s). Default : 0.0001)

m = 1000
MAXITER = 50
threshold = 0.0001
run_PSOET.run_psoet(n_site,n_elec,U,t,n_imp,approx,code_directory,
                    m = m,
                    MAXITER=MAXITER,
                    threshold = threshold)
