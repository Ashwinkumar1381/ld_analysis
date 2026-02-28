print("\033[J\033[H",end='') # Clear screen

import matplotlib.pyplot as plt
from ast import literal_eval as liteval
import os
import numpy as np
from file_utils import readData2D

plt.rc('text', usetex = True)

dirpath = "//media/ashwin/One Touch/ashwin_md/Apr2025/lane/lmp/Data19/"
fpath = os.path.join(dirpath, "rdf.dat")

rdf = readData2D(fpath, 4, 200, 200)

fig1, ax1 = plt.subplots()

z = rdf[3]
levels = np.linspace(z.min(), z.max(), 100)
cs = ax1.contourf(np.linspace(-10,10,200), np.linspace(-10,10,200), z,
                  levels = levels, cmap = 'bwr')
cbar = fig1.colorbar(cs)