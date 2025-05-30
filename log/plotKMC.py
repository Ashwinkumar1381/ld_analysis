print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
from ast import literal_eval as liteval
import math as m
import numpy as np
import file_utils

mpl.pyplot.rc('text', usetex=True)

runID = 3
dt = 5e-4
timeStart = 0e3
timeEnd = 30e3
plot_interval = 20

nr_list = {'0': {1: "0", 2: "1e-4", 3: "2e-4", 4: "5e-4", 5:"1e-3", 6:"2e-3", 7:"5e-3",
               8: "1e-2", 9: "2e-2", 10:"5e-2", 11:"0.1", 12:"0.2", 13:"0.5", 14:"1"},
           
           '1e-3': {15: "0", 16: "1e-4", 17: "2e-4", 18: "5e-4", 19:"1e-3", 20:"2e-3", 21:"5e-3",
               22: "1e-2", 23: "2e-2", 24:"5e-2", 25:"0.1", 26:"0.2", 27:"0.5", 28:"1"},
           
           '1e-2': {29: "0", 30: "1e-4", 31: "2e-4", 32: "5e-4", 33:"1e-3", 34:"2e-3", 35:"5e-3",
               36: "1e-2", 37: "2e-2", 38:"5e-2", 39:"0.1", 40:"0.2", 41:"0.5", 42:"1"},
           
           '1e-1': {43: "0", 44: "1e-4", 45: "2e-4", 46: "5e-4", 47:"1e-3", 48:"2e-3", 49:"5e-3",
               50: "1e-2", 51: "2e-2", 52:"5e-2", 53:"0.1", 54:"0.2", 55:"0.5", 56:"1"},
           
           '1e0': {57: "0", 58: "1e-4", 59: "2e-4", 60: "5e-4", 61:"1e-3", 62:"2e-3", 63:"5e-3",
               64: "1e-2", 65: "2e-2", 66:"5e-2", 67:"0.1", 68:"0.2", 69:"0.5", 70:"1"},
    
           '2': {17: "0", 10: "1e-4", 11: "5e-4", 12: "1e-3", 13: "5e-3",
                9: "1e-2", 14: "5e-2", 15: "0.1", 18: "0.5", 16: "1"},
           
           '5':  {7: "0", 6: "1e-4", 4: "5e-4", 3: "1e-3", 2: "5e-3",
                32: "1e-2", 1: "5e-2", 42: "0.1", 5: "0.5", 8: "1"}}

Pe = "1e0"

fig, ax = mpl.pyplot.subplots()

plot_list = [73]

for i in range(len(plot_list)):
    #fpath = "/media/ashwin/One Touch/ashwin_md/Feb2025/Pe_{pe}/Data{nr}/kmc.dat".format(pe = str(Pe), nr = plot_list[i])
    fpath = "//media/ashwin/One Touch/ashwin_md/Mar2025/Data{nr}/kmc.dat".format(nr = plot_list[i])
    figpath = "../../LD/LD-cpp/imgs/KMCplot2.png"

    [[step, numA, numB], _] = file_utils.readData(fpath, 1)
    
    step = [j*dt for j in step]
    step = step[0:len(step):plot_interval]
    numA = numA[0:len(numA):plot_interval]
    numB = numB[0:len(numB):plot_interval]
    
    y = [numB[j]/numA[j] for j in range(len(step))]

    ax.set(xlim = (timeStart, timeEnd))
    ax.set_xlabel(r"$t~(\sigma^2/D)$", fontsize = 14)
    ax.set_ylabel(r"$K = N_B/N_A$", fontsize = 14)
    #ax.plot(step, y, lw = 1, ls = "-", label = r"$\lambda =~$" + nr_list[Pe][plot_list[i]])
    ax.plot(step, y, lw = 1, ls = "-")

ax.axhline(y = 1, xmin = 0, xmax = 1, lw = 0.5, color = "k", ls = "-")
ax.axhline(y = 0, xmin = 0, xmax = 1, lw = 0.5, color = "k", ls = "-")
ax.legend(fontsize = 10)

#fig.savefig(figpath, dpi = 600, bbox_inches='tight')