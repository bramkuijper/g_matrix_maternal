#!/usr/bin/env python

import math

rmu = [ 0, 0.5, -0.5, 0.9, -0.9 ]

romega = [ 0, 0.5, -0.5, 0.9, -0.9 ]

omega = [ 3, 9, 49 ]

mu = 0.0002

mu_m = [ 0, 0.01 ]
sdmu_m = 0.02

B = 2

exe = "./xg_maternal"

replicates = 10

ctr = 1

for replicate_i in range(0,replicates):
    for mu_mi in mu_m:
        for rmu_i in rmu:
            for omega_i in omega:
                for romega_i in romega:
                    omega12 = romega_i * omega_i

                    print("echo " + str(ctr))
                    ctr += 1
                    print(exe 
                            + " " + str(math.sqrt(0.05)) + " " + str(math.sqrt(0.05))
                            + " " + str(rmu_i) + " " + str(mu)
                            + " " + str(mu_mi) + " " + str(sdmu_m)
                            + " " + str(omega_i) + " " + str(omega12) + " " + str(omega12)
                            + " " + str(omega_i) + " " + str(B))
