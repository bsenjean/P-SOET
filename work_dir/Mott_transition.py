#!/usr/bin/python 

L = 400
t = 1.0
n_imp = 4
approx = "iBALDA"
description = "_singleshot"
for U in [4.0,8.0]:
  output = "L{}_U{}_t{}_nimp{}_{}{}_MIT_results.out".format(str(L),str(U),str(t),str(n_imp),approx,description)
  with open("OUTPUT/"+output,"w") as f:
     f.write("%15s %15s\n" % ("n_elec","chem_pot"))
  for i in range(-100*int(U)-1,100*int(U)+1,1):
    mu = i/100.0
    energy_min = 1000
    for N in range(2*n_imp,402,2): 
      name = "L{}_N{}_U{}_t{}_nimp{}_{}{}.out".format(str(L),str(N),str(U),str(t),str(n_imp),approx,description)
      with open("OUTPUT/"+name,"r") as f:
        for i, lines in enumerate(f):
          if i == 3:
            energy = float(lines.split()[4+n_imp-1])
            energy_mu_n = energy - mu*N/L
            if energy_mu_n < energy_min:
              energy_min = energy_mu_n
              Nmin = N
    with open("OUTPUT/"+output,"a") as f:
       f.write("%15d %15.6f\n" % (Nmin, mu))
