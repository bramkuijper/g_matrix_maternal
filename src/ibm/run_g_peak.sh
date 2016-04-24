#!/usr/bin/env bash

#./xg_only  0.2236068  0.2236068 0.0  0.0002  49  44.1 44.1  49  2
#./xg_maternal  0.2236068  0.2236068 0.9  0.0002 0.01 0.02 49 44.1 44.1 49 2
#./xg_maternal_peakshift 0.22360679775 0.22360679775 0     0.0002 0.02 0.02    20 17 17 20 2       0.0017 -0.0017 1 0 0 &
./xg_maternal_peakshift 0.22360679775 0.22360679775 0     0.0002 0.02 0.02    20 17 17 20 2       0.0017 0.0017 1 0.01 0.01  &
./xg_maternal_peakshift 0.22360679775 0.22360679775 0     0.0002 0.02 0.02    20 17 17 20 2       0.2 0.2 100 0.0 0.0  &



#    a1 = atof(argv[1]);
#    a2 = atof(argv[2]);
#    rmu = atof(argv[3]);
#    mu = atof(argv[4]);
#    mu_m = atof(argv[5]);
#    sdmu_m = atof(argv[6]);
#    omega[0][0] = atof(argv[7]);
#    omega[0][1] = atof(argv[8]);
#    omega[1][0] = atof(argv[9]);
#    omega[1][1] = atof(argv[10]);
#    B = atof(argv[11]);
#    delta_t1 = atof(argv[12]);
#    delta_t2 = atof(argv[13]);
#    interval = atoi(argv[14]);
#    sigma_theta1 = atof(argv[15]);
#    sigma_theta2 = atof(argv[16]);
