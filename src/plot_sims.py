#!/usr/bin/env python3

from __future__ import print_function
import pandas as pd
import sys
import re, math, os
import numpy as np
import matplotlib.pyplot as plt


from matplotlib import rcParams
#rcParams['font.family'] = 'Myriad Pro'

filename = sys.argv[1]

f = open(filename);
fl = f.readlines();
f.close()

parline = -1

for idx, line in enumerate(fl):
    if re.match("^npop.*",line) != None:
        parline = idx - 1;
        break;



if parline > 0:
    dat = pd.read_csv(filename, nrows=parline-3, sep=";")
else:
    dat = pd.read_csv(filename, sep=";")

print(dat.describe())

# only take every tenth generation, otherwise too much data....
dat = dat[dat["generation"] % 10 == 0]

colnames = dat.columns.values

# number of graphs
num_rows = 5 

# when maternal effects are included extend the amount of rows
if "meanm11" in colnames:
    num_rows += 2

if "theta1" in colnames:
    num_rows += 1


plt.figure(figsize=(8,16))
plt.subplot(num_rows,1,1)
plt.grid(True)
plt.plot(dat["generation"], dat["meanz1"],"r",dat["generation"], dat["meanz2"],"b")
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'$z_{i}$')
plt.legend((r'$\bar{z}_{1}$',r'$\bar{z}_{2}$'), loc=2)

plt.subplot(num_rows,1,2)
plt.grid(True)
plt.plot(dat["generation"], dat["G11"],"r",dat["generation"], dat["G22"],"b",dat["generation"], dat["G12"],"g", dat["generation"], dat["G21"],"y")
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'genetic (co)variance $G_{ij}$')
plt.legend((r'$G_{11}$',r'$G_{22}$',r'$G_{12}$',r'$G_{21}$'),loc=2)

plt.subplot(num_rows,1,3)
plt.grid(True)
plt.plot(dat["generation"], dat["P11"],"c",dat["generation"], dat["P22"],"m",dat["generation"], dat["P12"],"y",dat["generation"], dat["P21"],"k")
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'phenotypic (co)variance $P_{ij}$')
plt.legend((r'$P_{11}$',r'$P_{22}$',r'$P_{12}$',r'$P_{21}$'),loc=2)

plt.subplot(num_rows,1,4)
plt.grid(True)
plt.plot(dat["generation"], dat["meanw"],"r")
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'survival $w$')

plt.subplot(num_rows,1,5)
plt.grid(True)
plt.plot(dat["generation"], dat["ev1"],"r",dat["generation"], dat["ev2"],"b")
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.legend((r'$\lambda_{1}$',r'$\lambda_{2}$'),loc=2)
plt.ylabel(r'eigenvalues $\lambda$')

row = 5

if "theta1" in colnames:
    row+=1
    plt.subplot(num_rows,1,row)
    plt.grid(True)
    plt.plot(dat["generation"], dat["theta1"],"r",dat["generation"], dat["theta2"],"b")
    plt.tick_params(axis='x',which='both',bottom='on',top='on')
    plt.legend((r'$\theta_{1}$',r'$\theta_{2}$'),loc=2)
    plt.ylabel(r'envt optima')

if "meanm11" in colnames:
    row += 1
    plt.subplot(num_rows,1,row)
    plt.grid(True)
    plt.plot(dat["generation"], dat["meanm11"],"r",dat["generation"], dat["meanm22"],"b", dat["generation"], dat["meanm12"], "m",dat["generation"], dat["meanm21"], "#ff8500")
    plt.tick_params(axis='x',which='both',bottom='on',top='on')
    plt.legend((r'$\bar{m}_{11}$',r'$\bar{m}_{22}$',r'$\bar{m}_{12}$',r'$\bar{m}_{21}$'),loc=2)
    plt.ylabel(r'mean maternal effect')

    row+=1
    plt.subplot(num_rows,1,row)
    plt.grid(True)
    plt.plot(dat["generation"], dat["Gm11"],"r",dat["generation"], dat["Gm22"],"b", dat["generation"], dat["Gm12"], "m", dat["generation"], dat["Gm21"], "#ff8500")
    plt.tick_params(axis='x',which='both',bottom='on',top='on')
    plt.legend((r'$\sigma_{{m}_{11}}^{2}$',r'$\sigma_{{m}_{22}}^{2}$',r'$\sigma_{{m}_{12}}^{2}$',r'$\sigma_{{m}_{21}}^{2}$'),loc=2)
    plt.ylabel(r'maternal effects variance')


graphname = "graph_" + os.path.basename(filename) + ".pdf"
plt.subplots_adjust(hspace=.3)
plt.savefig(graphname,format="pdf")

