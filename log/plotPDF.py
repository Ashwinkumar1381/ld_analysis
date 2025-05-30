print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import file_utils
import os
import numpy as np
from ast import literal_eval as liteval

mpl.pyplot.rc('text', usetex=True)

# Global params
markers = ["o", "o", "s", "^", "*", "D", "P"]
fs = ["full", "none", "full", "full", "full", "full", "full"]
colors = ["darkorange", "olive", "orangered","crimson", "#2544D0", "goldenrod", "k"]
#markers.remove("^")
#colors.remove("crimson")
lw = 0.5
ms = 4.0
fontsize = 18
fac = 1

option = "tavg_multi"
if(option == "frame"):
    # Params for frames-based data
    nr = 64
    Lx = 200
    Pe = "1.0"
    lamb = " 10^{-2}" 
    pdfSamples = [5e3, 30e3]
    nSamples = len(pdfSamples)
    sampleNum = [0, 1]
    
    fpath = "//media/ashwin/One Touch/ashwin_md/Feb2025/Pe_1e0/Data{nr}/pdf_frame.dat".format(nr = nr)
    figpath = "//media/ashwin/One Touch/ashwin_md/Feb2025/Pe_1e0/Data{nr}/pdf_frame_nr_{nr}.png".format(nr = nr)
    [pdfData, _] = file_utils.readData(fpath, nSamples)
    
    fig, ax = mpl.pyplot.subplots(1, nSamples, figsize = (12, 6))
    for i in sampleNum:
        x = pdfData[3*i + 0]
        pdfA = pdfData[3*i + 1]
        pdfB = pdfData[3*i + 2]
        
        x = [(j-Lx/(2*len(x))-Lx/2) for j in x]
        pdfA = [j*fac for j in pdfA]
        pdfB = [j*fac for j in pdfB]
        
        ax[i].plot(x, pdfA, "b-", marker = markers[i], lw = lw, ms = ms, label = r"$\rho_A$")
        ax[i].plot(x, pdfB, "r-", marker = markers[i], lw = lw, ms = ms, label = r"$\rho_B$")
    
    for i in range(nSamples):
        ax[i].set(xlim = (-Lx/2, Lx/2), ylim = (-0.01, 1.0))
        ax[i].set_xlabel(r"$x$", fontsize = fontsize)
        ax[i].set_ylabel(r"$\rho(x)$", fontsize = fontsize)
        ax[i].set_title(r"$t = {time}~\sigma^2/D$".format(time = int(pdfSamples[i])), 
                        fontsize = fontsize)
    
        ax[i].legend(fontsize = fontsize)

    #fig.suptitle(r"Simulation Parameters : $Pe_s = {Pe},~ \lambda = {lamb}$".format(Pe=Pe,lamb=lamb), fontsize = fontsize)
    #fig.savefig(figpath, dpi = 500, bbox_inches = "tight")

elif(option == "tavg"):
    Lx = 200
    nr = 264
    Pe = 1
    lamb = 1e-2
    
    fig, ax = mpl.pyplot.subplots(figsize = (6, 6))
    
    fpath = "//media/ashwin/One Touch/ashwin_md/Mar2025/Data{nr}/pdf_tavg.dat".format(nr = nr)
    figpath = "//media/ashwin/One Touch/ashwin_md/Mar2025/Data{nr}/pdf_tavg_nr{nr}".format(nr = nr)
    
    [[x, pdfA, pdfB], _] = file_utils.readData(fpath, 1)
    
    x = [(j - Lx/(2*len(x)) - Lx/2) for j in x]
    
    ax.plot(x, pdfA, "o-", lw = lw, ms = ms, color = "b", label = r"$\rho_A$")
    ax.plot(x, pdfB, "o-", lw = lw, ms = ms, color = "r", label = r"$\rho_B$")
    
    ax.set(xlim = (-Lx/2, Lx/2), ylim = (0, 1))
    ax.set_xlabel(r"$x (\sigma)$", fontsize = fontsize)
    ax.set_ylabel(r"$\rho (x)$", fontsize = fontsize)
    ax.set_title(r"$Pe = {pe}, ~\lambda = {lamb}$".format(pe = Pe, lamb = lamb), fontsize = fontsize)
    ax.legend(fontsize = fontsize)
    
    #fig.savefig(figpath, dpi = 500, bbox_inches = "tight")

elif(option == "tavg_multi"):
    Lx = 200
        
    fig, ax = mpl.pyplot.subplots(1, 2, figsize = (12, 6))
    
    plot_list = [80, 82, 84]
    m = [1e-2, 0.1, 1]
    
    dirpath = "//media/ashwin/One Touch/ashwin_md/Mar2025"
    figpath = "/home/ashwin/Desktop/ashwin_md/LD/LD-cpp/imgs/Mar2025"
    
    for i in range(len(plot_list)):
        fpath = os.path.join(dirpath, "Data{nr}".format(nr = plot_list[i]), "pdf_tavg.dat")
        
        [[x, pdfA, pdfB], _] = file_utils.readData(fpath, 1)
        
        x = [(j - Lx/(2*len(x)) - Lx/2) for j in x]
        
        ax[0].plot(x, pdfA, color = "b", marker = markers[i], lw = lw, ms = ms, 
                   label = "${m}$".format(m = m[i]))
        ax[1].plot(x, pdfB, color = "r", marker = markers[i], lw = lw, ms = ms,
                   label = "${m}$".format(m = m[i]))
    
    for i in range(len(ax)):
        ax[i].set(xlim = (-Lx/2, Lx/2), ylim = (0, 1))
        ax[i].set_xlabel(r"$x (\sigma)$", fontsize = fontsize)
        ax[i].legend(title = r"$m^*$", title_fontsize = fontsize - 4, 
                     fontsize = fontsize - 4, frameon = False)

    ax[0].set_ylabel(r"$\rho_A (x)$", fontsize = fontsize)
    ax[1].set_ylabel(r"$\rho_B (x)$", fontsize = fontsize)
    ax[0].set_title(r"$Pe_s = 10.0, ~\lambda = 10^{-3}$", fontsize = fontsize - 4, loc = "left")
    
    #fig.savefig(os.path.join(figpath, "pdf_tavg_Pe10.png"), dpi = 500, bbox_inches = "tight")
    
elif(option == "demix"):
    # Params for time-averaged data

    Lx = 200.0
    Ly = 50.0
    N = 5000
    
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
        
               '2e0': {17: "0", 10: "1e-4", 11: "5e-4", 12: "1e-3", 13: "5e-3",
                    9: "1e-2", 14: "5e-2", 15: "0.1", 18: "0.5", 16: "1"},
               
               '5e0':  {7: "0", 6: "1e-4", 4: "5e-4", 3: "1e-3", 2: "5e-3",
                    32: "1e-2", 1: "5e-2", 42: "0.1", 5: "0.5", 8: "1"}}

    Pe = "1e0"
    
    plot_list = [57] + list(np.arange(60, 71))
    demix = np.zeros(len(plot_list))
    
    for i in range(len(plot_list)):
        fpath = "/media/ashwin/One Touch/ashwin_md/Feb2025/Pe_{pe}/Data{nr}/pdf_tavg.dat".format(pe = str(Pe), 
                                                            nr = plot_list[i])
        figpath = "../../LD/LD-cpp/imgs/demix1.png"
        
        [[x, pdfA, pdfB], _] = file_utils.readData(fpath, 1)
        
        x = [(j - Lx/(2*len(x)) - Lx/2) for j in x]
        pdfA = [j*fac for j in pdfA]
        pdfB = [j*fac for j in pdfB]
        pdfALL = [(pdfA[j] + pdfB[j]) for j in range(len(x))]
        
        integ = [abs(pdfA[i] - pdfB[i]) for i in range(len(x))]
        demix[i] = Ly*np.trapz(integ, x)/N

        print(round(demix[i], 5))