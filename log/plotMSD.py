print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
from ast import literal_eval as liteval
import math as m
import file_utils

mpl.pyplot.rc('text', usetex=True)
fig, ax = mpl.pyplot.subplots()

# Plotting LAMMPS output MSD data
dt = 5e-5
frameW = 2e5

dirpath = "//media/ashwin/One Touch/ashwin_md/Nov2024/msd_benchmarks"
contents = os.listdir(dirpath)

files = []
for fname in contents:
    if(fname.find("lmp_D_") >= 0):
        files.append(fname)

for file in files:
    fpath = os.path.join(dirpath, file, "msd.dat")
    
    [[frame, msd], _] = file_utils.readData(fpath, 1)
    time = [i*dt*frameW for i in frame]
    
    ax.loglog(time, msd)
    