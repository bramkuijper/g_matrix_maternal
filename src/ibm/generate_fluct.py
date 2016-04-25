#!/usr/bin/env python

import math
import numpy as np

rmu = [ 0, 0.9, -0.9 ]

romega = [ 0, 0.9, -0.9 ]

omega = [ 0.7 ]

mu = 0.0002

mu_m = [ 0, 0.01 ]
sdmu_m = 0.02

sigma_e = math.sqrt(0.1)

B = 2

exe = "./xg_m_fluctuate_surv"
 
ampl1 = [ 1 ]
ampl2 = [ 1 ] 

stoch1 = [ 0 ]
stoch2 = [ 0 ]

freq1 = list(np.arange(0, math.pi + math.pi/20, math.pi/20))

shift = [ 0, 0.5 ]

replicates = 5

ctr = 1

for replicate_i in range(0,replicates):
    for mu_mi in mu_m:
        for rmu_i in rmu:
            for omega_i in omega:
                for romega_i in romega:
                    for ampl1_i in ampl1:
                        for ampl2_i in ampl2:
                            for stoch1_i in stoch1:
                                for stoch2_i in stoch2:
                                    for freq1_i in freq1:
                                        freq2_i = freq1_i
                                        for shift_i in shift:

                                            print("echo " + str(ctr))
                                            ctr += 1
                                            print(exe 
                                                    + " " + str(math.sqrt(0.05)) + " " + str(math.sqrt(0.05))
                                                    + " " + str(rmu_i) + " " + str(mu)
                                                    + " " + str(mu_mi) + " " + str(sdmu_m)
                                                    + " " + str(omega_i) + " " + str(omega_i) 
                                                    + " " + str(romega_i) + " " + str(ampl1_i)
                                                    + " " + str(ampl2_i) + " " + str(stoch1_i)
                                                    + " " + str(stoch2_i) + " " + str(freq1_i)
                                                    + " " + str(freq2_i) + " " + str(shift_i)
                                                    + " " + str(sigma_e)
                                            )
