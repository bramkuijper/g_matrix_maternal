#!/usr/bin/env python

import math

rmu = [ 0, 0.5, -0.5, 0.9, -0.9 ]

romega = [ 0.02, 0.5, -0.5, 0.9, -0.9 ]

omega = [ 9, 49 ]

delta_t1 = [ 0, 0.0017 ]
delta_t2 = [ 0, 0.0017, -0.0017 ]

interval = [1, 100, 1000]

sigma_theta = [ 0, 0.05 ]

mu = 0.0002

mu_m = [ 0, 0.01 ]
sdmu_m = 0.02

B = 2

exe = "./xg_maternal_peakshift"

replicates = 10

ctr = 1

for replicate_i in range(0,replicates):
    for mu_mi in mu_m:
        for rmu_i in rmu:
            for omega_i in omega:
                for romega_i in romega:
                    for interval_i in interval:
                        if interval_i > 1:
                            delta_t1_sub = [ i * interval_i for i in delta_t1]
                            delta_t2_sub = [ i * interval_i for i in delta_t2]
                        else:
                            delta_t1_sub = delta_t1
                            delta_t2_sub = delta_t2

                        for delta_t1_i in delta_t1_sub:
                            for delta_t2_i in delta_t2_sub:

                                if delta_t1_i == 0 and delta_t2_i == 0:
                                    continue

                                for sigma_theta_i in sigma_theta:
                                    omega12 = romega_i * omega_i

                                    print("echo " + str(ctr))
                                    ctr += 1
                                    print(exe 
                                            + " " + str(math.sqrt(0.05)) + " " + str(math.sqrt(0.05))
                                            + " " + str(rmu_i) + " " + str(mu)
                                            + " " + str(mu_mi) + " " + str(sdmu_m)
                                            + " " + str(omega_i) + " " + str(omega12) + " " + str(omega12)
                                            + " " + str(omega_i) + " " + str(B)
                                            + " " + str(delta_t1_i) + " " + str(delta_t2_i)
                                            + " " + str(interval_i)
                                            + " " + str(sigma_theta_i) + " " + str(sigma_theta_i)
                                            )
