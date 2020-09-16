#!/usr/bin/python 
import sys
import os

code_directory = os.getenv('PSOET_DIR')+"/code/"
work_directory = os.getcwd()
dmrg_code      = os.getenv('DMRG_code')

sys.path.append(code_directory)
import run_PSOET

n_site      = 12
n_elec      = 6
U           = 8.0
t           = 1.0
n_imp       = 2 # Solved analytically for n_imp = 1. DMRG impurity solver is called for n_imp > 1
approx      = ["iBALDA", "2LBALDA", "SIAMBALDA"][0]
m           = 1000 # default is 1000. Number of renormalized states in DMRG.
semi_opt    = False # default is False. This means only dEcimp_dn(n_0) will vary, while deHxc_dn and dEHximp_dn are calculated with fixed n=N/L.
chem_pot    = False # default is False. Optimized a chemical potential numerically to match the occupation between the impurity wavefunction and the KS one, like in DMET. Works for a single impurity only. 
single_shot = False # default is False. True --> MAXITER = 1
description = "" # default is "". Additional description added to the name of the output file. Note that semiopt, chempot and ss (single-shot) are already set by default if they are set to True.
MAXITER     = 2 # default is 50
threshold   = 0.0001 # default is 0.0001

run_PSOET.run_psoet(n_site,n_elec,U,t,n_imp,approx,code_directory,work_directory,dmrg_code, # this first line are mandatory arguments.
                    m = m,
                    semi_opt = semi_opt,
                    chem_pot = chem_pot,
                    single_shot = single_shot,
                    MAXITER = MAXITER,
                    threshold = threshold)
