#!/usr/bin/python 
import os

L = 400
t = 1.0
n_imp = 1
approx = "2LBALDA"
description ="_singleshot"
Ulist = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,150.0,200.0,1000.0,10000.0]
chem_pot = False
for N in [400]:
 output = "L{}_N{}_t{}_nimp{}_{}{}_fctU_results.out".format(str(L),str(N),str(t),str(n_imp),approx,description)
 with open("OUTPUT/"+output,"w") as f:
    if chem_pot is True:
      f.write("%15s %15s %15s %15s %15s %15s %15s\n" % ("U/t","SS_energy","SS_dblocc","energy","dblocc","mu","deltav"))
    else:
      f.write("%15s %15s %15s %15s %15s\n" % ("U/t","SS_energy","SS_dblocc","energy","dblocc"))
 for U in Ulist:
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
         f.write("%15f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n" % (U,float(single_shot[4+n_imp-1]),float(single_shot[6+n_imp-1]),float(conv[3+n_imp-1]),float(conv[5+n_imp-1]),float(single_shot[7+n_imp-1]),float(single_shot[8+n_imp-1])))
       else:
         f.write("%15f %15.8f %15.8f %15.8f %15.8f\n" % (U,float(single_shot[4+n_imp-1]),float(single_shot[6+n_imp-1]),float(conv[3+n_imp-1]),float(conv[5+n_imp-1])))
