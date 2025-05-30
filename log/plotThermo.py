print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
from ast import literal_eval as liteval
import math as m
import file_utils

mpl.pyplot.rc('text', usetex=True)

dt = 5e-4

fpath = os.getcwd().removesuffix("/ld_analysis/log") + "/LD/LD-cpp/Data42/thermo3.dat"

[step, pe, ke, etot, T] = file_utils.readData(fpath, 1)

fig, ax = mpl.pyplot.subplots()

plot_interval = 10
L = len(step)

step = step[0:L:plot_interval]
time = [i*dt for i in step]
ke = ke[0:L:plot_interval]
pe = pe[0:L:plot_interval]
etot = etot[0:L:plot_interval]

ax.semilogy(time, ke, lw= 0.8, label = "ke") 
ax.semilogy(time, pe, lw= 0.8, label = "pe")
ax.semilogy(time, etot, lw= 0.8, label = "etot")

ax.set_xlabel("Step")
ax.set_ylabel("Per atom energy")
ax.legend()

figpath = fpath.removesuffix("/LD-cpp/Data1/thermo.dat") + "/imgs"
#fig.savefig(figpath + "/nve.png", dpi = 600, bbox_inches = 'tight')