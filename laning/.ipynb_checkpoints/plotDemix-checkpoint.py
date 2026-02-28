print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import file_utils
import numpy as np
from ast import literal_eval as liteval

mpl.pyplot.rc('text', usetex=True)

# Global params
markers = ["o", "s", "^", "D", "P", "p", "h"]
fs = ["none", "none", "full", "full", "full", "full", "full"]
colors = ["darkorange", "olive", "orangered","crimson", "#2544D0", "goldenrod", "k"]

Pe_list = [0, 1e-3, 1e-2, 0.1,1.0, 2.0, 5.0]

fpath = "../../LD/LD-cpp/demix.dat"
figpath = "../../LD/LD-cpp/imgs/demix3B.png"

data, Pe = file_utils.readData(fpath, len(Pe_list))

fig, ax = mpl.pyplot.subplots()
lines = []
points = []

for i in range(len(Pe_list)):
    lamb = data[2*i]
    line = ax.axhline(y = data[2*i + 1][0], xmin = 0, xmax = 0.08, linestyle = "--", color = colors[i], 
              label = r"$Pe = {pe}$".format(pe = Pe_list[i]))
    
    point, = ax.semilogx(lamb[1:], data[2*i + 1][1:], marker = markers[i], color = colors[i], 
                fillstyle = "none", lw = 0.5, label = r"$Pe = {pe}$".format(pe = Pe_list[i]))

    lines.append(line)
    points.append(point)

leg1 = ax.legend(handles = points, title = r"$(\bar{\lambda}\neq 0)$", fontsize = 12, 
                 bbox_to_anchor = (1, 1))
leg2 = mpl.pyplot.legend(handles = lines, title = r"$(\bar{\lambda}=0)$", fontsize = 12,
                 bbox_to_anchor = (0.78, 1))
ax.add_artist(leg1)

ax.set_xlabel(r"$\bar{\lambda}$", fontsize = 18)
ax.set_ylabel(r"$\delta$", fontsize = 18)

fig3, ax3 = mpl.pyplot.subplots()
lines = []
points = []

ax3.axhline(y = 0.0, xmin = 0, xmax = 1, linestyle = "--", color = "k", lw = 1)
ax3.axvline(x = 1.0, ymin = 0, ymax = 1, linestyle = "--", color = "k", lw = 1)
ax3.axvline(x = 2e-3, ymin = 0, ymax = 1, linestyle = "--", color = "k", lw = 1)

for i in range(1, len(Pe_list)):
    scaledlamb = [j/Pe_list[i] for j in data[2*i]]
    
    line = ax3.axhline(y = data[2*i + 1][0], xmin = 0, xmax = 0.05, linestyle = "--", color = colors[i-1])
    point, = ax3.semilogx(scaledlamb[1:], data[2*i + 1][1:], marker = markers[i], color = colors[i-1], 
                          linestyle = "-", lw = 0.8, label = r"$Pe_s = {pe}$".format(pe = Pe_list[i]))

ax3.set_xlabel(r"$\bar{\lambda} / Pe_s$", fontsize = 18)
ax3.set_ylabel(r"$\delta$", fontsize = 18)
ax3.legend(fontsize = 12)

"""
lamb_list = [0, 5e-4, 1e-3, 5e-3, 0.1, 1.0]

fig2, ax2 = mpl.pyplot.subplots()
data2 = [[] for i in range(len(lamb_list))]
lines2 = []
points2 = []

for i in range(len(Pe)):
    ctr = 0
    for j in range(len(data[2*i])):
        
        if(data[2*i][j] == lamb_list[ctr]):
            data2[ctr].append(data[2*i+1][j])
            ctr += 1
            
        if(ctr == len(lamb_list)): continue

for i in range(len(lamb_list)):
    
    point, = ax2.semilogx(Pe[1:], data2[i][1:], marker = markers[i], color = colors[i],
                fillstyle = "none", lw = 0.5, label = r"$\lambda = {lamb}$".format(lamb = lamb_list[i]))
    
    points2.append(point)

ax2.legend(fontsize = 12)
ax2.set_xlabel(r"$Pe$", fontsize = 15)
ax2.set_ylabel(r"$\delta$", fontsize = 15)
""" 

#fig.savefig(figpath, dpi = 500, bbox_inches = "tight")
#fig2.savefig(figpath, dpi = 500, bbox_inches = "tight")
#fig3.savefig(figpath, dpi = 500, bbox_inches = "tight")