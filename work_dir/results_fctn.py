#!/usr/bin/python 
import os

L = 400
t = 1.0
n_imp = 4
approx = "iBALDA"
description ="_singleshot"
chem_pot = False
for U in [4.0,8.0]:
 output = "L{}_U{}_t{}_nimp{}_{}{}_fctn_results.out".format(str(L),str(U),str(t),str(n_imp),approx,description)
 with open("OUTPUT/"+output,"w") as f:
    if chem_pot is True:
      f.write("%15s %15s %15s %15s %15s %15s %15s\n" % ("N_elec","SS_energy","SS_dblocc","energy","dblocc","mu","deltav"))
    else:
      f.write("%15s %15s %15s %15s %15s\n" % ("N_elec","SS_energy","SS_dblocc","energy","dblocc"))
 for N in  range(2*n_imp,402,2):
    name = "L{}_N{}_U{}_t{}_nimp{}_{}{}.out".format(str(L),str(N),str(U),str(t),str(n_imp),approx,description)
    if os.stat("OUTPUT/"+name).st_size == 0:
       raise ValueError("File " + name + " is empty.")
    with open("OUTPUT/"+name,"r") as f:
       for i, lines in enumerate(f):
         if i == 3:
            single_shot = lines.split()
         else:
            pass
       conv = os.popen("tail -n 1 %s" % "OUTPUT/"+name).read().split()
    with open("OUTPUT/"+output,"a") as f:
       if chem_pot is True:
         f.write("%15d %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n" % (N,float(single_shot[4+n_imp-1]),float(single_shot[6+n_imp-1]),float(conv[3+n_imp-1]),float(conv[5+n_imp-1]),float(single_shot[7+n_imp-1]),float(single_shot[8+n_imp-1])))
       else:
         f.write("%15d %15.8f %15.8f %15.8f %15.8f\n" % (N,float(single_shot[4+n_imp-1]),float(single_shot[6+n_imp-1]),float(conv[3+n_imp-1]),float(conv[5+n_imp-1])))
