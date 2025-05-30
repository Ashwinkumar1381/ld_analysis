print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import os
from ast import literal_eval as liteval
import math as m
import numpy as np
import file_utils

mpl.pyplot.rc('text', usetex=True)

# Global params
markers = ["o", "s", "^", "D", "*"]
colors = ["crimson", "#2544D0", "darkorange", "olive", "orangered"]
lw = 0.5
ms = 1
fs = "full"
fontsize = 12

option = "tavg"
if(option == "tavg"):
    fontsize = 15
    nr_list = [32]
    Pe_list = [5]
    lamb_list = [1e-2]
    Lx = 200
    Ly = 100
    binW = 4
    
    fig1, [(ax1, ax2), (ax3, ax4)] = mpl.pyplot.subplots(2,2)
    
    ax1.axvline(x = 0, ymin = 0, ymax = 1, color = "k", ls= "--", lw = lw+0.5, alpha = 0.6)
    ax2.axvline(x = 0, ymin = 0, ymax = 1, color = "k", ls= "--", lw = lw+0.5, alpha = 0.6)
    ax3.axvline(x = 0, ymin = 0, ymax = 1, color = "k", ls= "--", lw = lw+0.5, alpha = 0.6)
    ax4.axvline(x = 0, ymin = 0, ymax = 1, color = "k", ls= "--", lw = lw+0.5, alpha = 0.6)
    
    for i in range(len(nr_list)):
        
        fpath = "../../LD/LD-cpp/Data{nr}/pressure_tavg.dat".format(nr = nr_list[i])
        fpath2 = "../../LD/LD-cpp/Data{nr}/pdf_tavg.dat".format(nr = nr_list[i])
        figpath = "../../LD/LD-cpp/imgs/pressure_tavg_1_nr32.png"
        
        [x, Pin_xx, Pin_xy, Pin_yx, Pin_yy, Pkin_xx, Pkin_xy, Pkin_yx, Pkin_yy,
         Pswim] = file_utils.readData(fpath, 1)
        
        [x2, pdfA, pdfB] = file_utils.readData(fpath2, 1)
        
        x = [(j-Lx/(2*len(x))-Lx/2) for j in x]
        Pswim = [(Pe_list[0]*binW*i) for i in pdfB]
        Ptot_xx = [(Pin_xx[i]+Pkin_xx[i]+Pswim[i]) for i in range(len(x))]
        Ptot_yy = [(Pin_yy[i]+Pkin_yy[i]) for i in range(len(x))]
        Pdiff = [(Ptot_xx[i] - Ptot_yy[i]) for i in range(len(x))]
        
        ax1.plot(x, Pin_xx, marker = markers[0], ms = ms, color = colors[0], ls = "-", 
                lw = lw, fillstyle = fs, label = r"$\displaystyle\mathcal{P}^{(in)}_{xx}$")
        ax1.plot(x, Pin_yy, marker = markers[1], ms = ms, color = colors[0], ls = "-", 
                lw = lw, fillstyle = fs, label = r"$\displaystyle\mathcal{P}^{(in)}_{yy}$")
        ax2.plot(x, Pswim, marker = markers[0], ms = ms, color = colors[1], ls = "-", 
                lw = lw, fillstyle = fs, label = r"$\displaystyle\mathcal{P}^{(swim)}_{xx}$")
        ax3.plot(x, Ptot_xx, marker = markers[0], ms = ms, color = colors[2], ls = "-", 
                lw = lw, fillstyle = fs, label = r"$\displaystyle\mathcal{P}_{N}$")
        ax3.plot(x, Ptot_yy, marker = markers[1], ms = ms, color = colors[3], ls = "-", 
                lw = lw, fillstyle = fs, label = r"$\displaystyle\mathcal{P}_{T}$")
        ax4.plot(x, Pdiff, marker = markers[2], ms = ms, color = colors[4], ls = "-", 
                lw = lw, fillstyle = fs, label = r"$\displaystyle\mathcal{P}_{N}-\mathcal{P}_{T}$")
            
    ax1.set(xlim = (-Lx/2,Lx/2), ylim = (-2,max(Ptot_xx)+5), aspect = "auto")
    ax1.set_xlabel(r"$x/\sigma$", fontsize = fontsize)
    ax1.set_ylabel(r"$\mathcal{P}\sigma^3/k_BT$", fontsize = fontsize)
    ax1.legend(fontsize = fontsize - 5)
    
    ax2.set(xlim = (-Lx/2, Lx/2), ylim = (-0.2,max(Pswim)+2), aspect = "auto")
    ax2.set_xlabel(r"$x/\sigma$", fontsize = fontsize)
    ax2.set_ylabel(r"$\mathcal{P}\sigma^3/k_BT$", fontsize = fontsize)
    ax2.legend(fontsize = fontsize - 5)
    
    ax3.set(xlim = (-Lx/2, Lx/2), ylim = (-2,max(Ptot_xx)+5), aspect = "auto")
    ax3.set_xlabel(r"$x/\sigma$", fontsize = fontsize)
    ax3.set_ylabel(r"$\mathcal{P}\sigma^3/k_BT$", fontsize = fontsize)
    ax3.legend(fontsize = fontsize - 5)
    
    ax4.set(xlim = (-Lx/2, Lx/2), ylim = (-0.2,max(Pdiff)+2), aspect = "auto")
    ax4.set_xlabel(r"$x/\sigma$", fontsize = fontsize)
    ax4.set_ylabel(r"$\mathcal{P}\sigma^3/k_BT$", fontsize = fontsize)
    ax4.legend(fontsize = fontsize - 5)
    
    fig1.set_figheight(6)
    fig1.set_figwidth(9)
    fig1.tight_layout(pad = 1)
    #fig1.suptitle(r"$nr = {nr}, \lambda = {lamb}, Pe_s = {Pe}$".format(nr=nr_list[0],lamb=lamb_list[0],Pe=Pe_list[0]))
    
    ax1.set_facecolor("#ffffff")
    ax2.set_facecolor("#ffffff")
    ax3.set_facecolor("#ffffff")
    fig1.set_facecolor("#ffffff")
    #fig1.savefig(figpath, dpi = 600, bbox_inches = "tight")

if(option == "tension"):
    
    Lx = 200
    Ly = 100
    binW = 4
    
    fig2, [ax1, ax2] = mpl.pyplot.subplots(1,2)
    
    lamb_list = [0.001, 0.01, 0.1, 1]
    for i in range(2, len(lamb_list)+2):
        Pe_list = [1, 5, 10, 25, 50]
        nr_list = [i*10 + j for j in range(1, len(Pe_list)+1)]
    
        surface_tension = []
        for j in range(len(Pe_list)):
            fpath = "../../LD/LD-cpp/Data{nr}/pressure_tavg.dat".format(nr=nr_list[j])
            fpath2 = "../../LD/LD-cpp/Data{nr}/pdf_tavg.dat".format(nr = nr_list[j])
            figpath = "../../LD/LD-cpp/imgs/pressure_tavg_2.png"
        
            [x, Pin_xx, Pin_xy, Pin_yx, Pin_yy, Pkin_xx, Pkin_xy, Pkin_yx, Pkin_yy,
             Pswim] = file_utils.readData(fpath, 1)
            
            [x2, pdfA, pdfB] = file_utils.readData(fpath2, 1)
        
            x = [(k-Lx/(2*len(x))-Lx/2) for k in x]
            Pswim = [(Pe_list[j]*binW*k) for k in pdfB]
            Ptot_xx = [(Pin_xx[k]+Pkin_xx[k]+Pswim[k]) for k in range(len(x))]
            Ptot_yy = [(Pin_yy[k]+Pkin_yy[k]) for k in range(len(x))]
            Pdiff = [(Ptot_xx[k] - Ptot_yy[k]) for k in range(len(x))]
        
            surface_tension.append(np.trapz(Pdiff,x)/2)
            print("nr = {nr}, Pe = {Pe}, rate = {rate}, surface tension = {gamma}".format(nr = nr_list[j], 
                                Pe=Pe_list[j], rate=lamb_list[i-2], gamma = surface_tension[j]))
            
        ax1.loglog(Pe_list, surface_tension, marker = markers[0], ms = ms+2, color = colors[i-2], ls = "-", 
                lw = lw+0.1, fillstyle = fs, label = str(lamb_list[i-2]))
    
    ax1.loglog(())
    ax1.set(xlim = (0.9, 60), ylim = (10, 5e4))
    ax1.set_xlabel(r"$Pe_s = v_s\sigma/D$", fontsize = fontsize)
    ax1.set_ylabel(r"$\gamma \sigma^2/k_BT$", fontsize = fontsize)
    ax1.legend(fontsize = fontsize-4, title = r"$\bar{\lambda}$", loc = "upper right", ncol = 2)
    
    print("\n")
    
    Pe_list = [1, 5, 10, 25, 50]
    for i in range(1, len(Pe_list)+1):
        lamb_list = [0.001, 0.01, 0.1, 1]
        nr_list = [j*10 + i for j in range(2, len(lamb_list)+2)]
    
        surface_tension = []
        for j in range(len(nr_list)):
            fpath = "../../LD/LD-cpp/Data{nr}/pressure_tavg.dat".format(nr=nr_list[j])
            fpath2 = "../../LD/LD-cpp/Data{nr}/pdf_tavg.dat".format(nr = nr_list[j])
            figpath = "../../LD/LD-cpp/imgs/pressure_tavg_2.png"
        
            [x, Pin_xx, Pin_xy, Pin_yx, Pin_yy, Pkin_xx, Pkin_xy, Pkin_yx, Pkin_yy,
             Pswim] = file_utils.readData(fpath, 1)
            
            [x2, pdfA, pdfB] = file_utils.readData(fpath2, 1)
        
            x = [(k-Lx/(2*len(x))-Lx/2) for k in x]
            Pswim = [(Pe_list[i-1]*binW*k) for k in pdfB]
            Ptot_xx = [(Pin_xx[k]+Pkin_xx[k]+Pswim[k]) for k in range(len(x))]
            Ptot_yy = [(Pin_yy[k]+Pkin_yy[k]) for k in range(len(x))]
            Pdiff = [(Ptot_xx[k] - Ptot_yy[k]) for k in range(len(x))]
        
            surface_tension.append(np.trapz(Pdiff,x)/2)
        
            """
            Pdiff_avg = [0.5*(Pdiff[k]+Pdiff[len(x)-1-k]) for k in range(int(len(x)/2))]
            Pdiff_avg.reverse()
            x = x[int(len(x)/2):len(x)]
            surface_tension.append(np.trapz(Pdiff_avg,x))
            """
            
            print("nr = {nr}, Pe = {Pe}, rate = {rate}, surface tension = {gamma}".format(nr = nr_list[j], 
                                Pe=Pe_list[i-1], rate=lamb_list[j], gamma = surface_tension[j]))
            
        ax2.loglog(lamb_list, surface_tension, marker = markers[0], ms = ms+2, color = colors[i-1], ls = "-", 
                lw = lw+0.1, fillstyle = fs, label = str(Pe_list[i-1]))
    ax2.set(ylim = (10, 6e4))
    ax2.set_xlabel(r"$\bar{\lambda} = \lambda\sigma^2/D$", fontsize = fontsize)
    ax2.set_ylabel(r"$\gamma \sigma^2/k_BT$", fontsize = fontsize)
    ax2.legend(fontsize = fontsize-4, title = r"$Pe_s$", loc = "upper right", ncol = 3)
    
    fig2.set_figheight(5)
    fig2.set_figwidth(12)
    fig2.tight_layout(pad = 1.5)
    ax1.set_facecolor("#ffffff")
    ax2.set_facecolor("#ffffff")
    fig2.set_facecolor("#ffffff")
    #fig2.savefig(figpath, dpi = 600, bbox_inches = "tight")
