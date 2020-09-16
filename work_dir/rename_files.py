import subprocess

L=400
t=1.0
description="_singleshot"
Ulist = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,150.0,200.0,1000.0,10000.0]

for approx in ["iBALDA"]:
 for U in Ulist:
   for n_imp in [1,4]:
    for N in [400]:
     name1 = "OUTPUT/L{}_N{}_U{}_t{}_nimp{}_{}.out".format(str(L),str(N),str(U),str(t),str(n_imp),approx)
     name2 = "OUTPUT/L{}_N{}_U{}_t{}_nimp{}_{}{}.out".format(str(L),str(N),str(U),str(t),str(n_imp),approx,description)
     subprocess.call("mv -f " + name1 + " " + name2,shell=True)

