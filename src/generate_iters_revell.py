#!/usr/bin/env python

import math

omega = [ 3, 9, 49 ]

sigma_r_omega = [ 0, 0.001, 0.01, 0.1 ]
sigma_r_mu = [ 0, 0.001, 0.01, 0.1 ]

mu = 0.0002

mu_m = [ 0, 0.01 ]
sdmu_m = 0.02

B = 2

exe = "./xg_m_revell"

replicates = 20

ctr = 1

for replicate_i in range(0,replicates):
    for mu_mi in mu_m:
        for omega_i in omega:
            for sigma_r_omega_i in sigma_r_omega:
                for sigma_r_mu_i in sigma_r_mu:
                    print("echo " + str(ctr))
                    ctr = ctr + 1
                    print(exe 
                        + " " + str(math.sqrt(0.05)) + " " + str(math.sqrt(0.05))
                        + " " + str(0) + " " + str(mu)
                        + " " + str(mu_mi) + " " + str(sdmu_m)
                        + " " + str(omega_i) 
                        + " " + str(omega_i) 
                        + " " + str(0) 
                        + " " + str(sigma_r_omega_i)
                        + " " + str(sigma_r_mu_i) 
                        )
