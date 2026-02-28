"""
Extracts columns of thermo data from log.lammps and writes to thermo.dat

Author          : Ashwin Kumar M
Date created    : 16.10.24
Last modified   : 03.02.25
"""

print("\033[J\033[H",end='') # Clear screen

from file_utils import *

nRuns = 2

fname1 = "../../../TA-CH5041/Data8/log.lammps"
fobj1 = open(fname1, mode = "r", encoding = "utf-8")

fname2 = fname1.removesuffix("/log.lammps") + "/thermo.dat"
fobj2 = open(fname2, mode = "w")

run_nr = 0
for line in fobj1:
    array = line2array(line)
    
    if(run_nr == nRuns):
        if(array[0] == 'Loop'): break
        else:
            line = array2line(array) 
            fobj2.write(line + "\n")
    
    if(array[0] == 'Per'): 
        run_nr += 1

fobj1.close()
fobj2.close()