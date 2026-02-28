print("\033[J\033[H",end='') # Clear screen

import matplotlib as mpl
import file_utils

mpl.pyplot.rc('text', usetex=True)

# Global params
markers = ["o", "s", "^", "*", "D", "P"]
colors = ["darkorange", "olive", "orangered","crimson", "#2544D0", "goldenrod", ]
lw = 0.3
ms = 2.0
fontsize = 12
fac = 1

option = "tavg"
if(option == "frames"):
    # Params for frames-based data
    nr = 7
    Lx = 200
    Pe = "5.0"
    lamb = " 10^{-3}" 
    pdfSamples = [40e3, 42e3, 44e3, 47e3, 49e3]
    nSamples = len(pdfSamples)
    sampleNum = [0,1,2,3,4]
    
    fpath = "../../LD/LD-cpp/Data{nr}/pdf_frames.dat".format(nr = nr)
    figpath = "../../LD/LD-cpp/imgs/pdf_frames_{nr}.png".format(nr = nr)
    pdfData = file_utils.readData(fpath, nSamples)
    
    fig, [axA, axB] = mpl.pyplot.subplots(1, 2)
    for i in sampleNum:
        x = pdfData[3*i + 0]
        pdfA = pdfData[3*i + 1]
        pdfB = pdfData[3*i + 2]
        
        x = [(j-Lx/(2*len(x))-Lx/2) for j in x]
        pdfA = [j*fac for j in pdfA]
        pdfB = [j*fac for j in pdfB]
        
        axA.plot(x, pdfA, "b-", marker = markers[i], lw = lw, ms = ms, label = str(pdfSamples[i]))
        axB.plot(x, pdfB, "r-", marker = markers[i], lw = lw, ms = ms, label = str(pdfSamples[i]))
    
    axA.set(xlim = (-Lx/2, Lx/2), ylim = (-0.01, 1.0))
    axB.set(xlim = (-Lx/2, Lx/2), ylim = (-0.01, 1.0))
    axA.set_xlabel(r"$x$", fontsize = fontsize)
    axB.set_xlabel(r"$x$", fontsize = fontsize)
    axA.set_ylabel(r"$\rho_{A} (x), ~\rho_{B} (x)$", fontsize = fontsize)
    axA.legend(title = "time", fontsize = fontsize - 3)
    axB.legend(title = "time", fontsize = fontsize - 3)
    fig.suptitle(r"Simulation Parameters : $Pe_s = {Pe},~ \lambda = {lamb}$".format(Pe=Pe,lamb=lamb), 
                 fontsize = fontsize)
    
    #fig.savefig(figpath, dpi = 600, bbox_inches = "tight")
    
elif(option == "tavg"):
    # Params for time-averaged data
    nr_list = [7, 6, 5, 4, 3, 2, 32, 42, 8]
    lamb = ["0", "1e-4", "2.5e-4", "5e-4", "7.5e-4", "1e-3", "0.01", "0.1", "1"]
    width = [16,12,12,8,4,4,4,2,0]
    Lx = 200.0
    
    fig, ax = mpl.pyplot.subplots(1, 3)
    fig2, ax2 = mpl.pyplot.subplots(1, 2)
    plot_list = [7,6,5,4,3,2]
    plot_lamb = [lamb[nr_list.index(j)] for j in plot_list]
    x_range = [range(34,50,4), range(30,46,4), range(30,46,4), range(22,46,4)]
    
    plot_data = [
        {"nr":3, "lamb":7.5e-4, "width":12, "x_range":range(34,50,4)},
        {"nr":4, "lamb":5.0e-4, "width":12, "x_range":range(30,46,4)},
        {"nr":5, "lamb":2.5e-4, "width":12, "x_range":range(30,46,4)},
        {"nr":6, "lamb":1.0e-4, "width":20, "x_range":range(22,46,4)},
    ]
    
    for i in range(len(plot_list)):
        fpath = "../../LD/LD-cpp/Data{nr}/pdf_tavg.dat".format(nr = plot_list[i])
        figpath = "../../LD/LD-cpp/imgs/pdf_tavg3_3.png"
        
        pdfData = file_utils.readData(fpath, 1)
        
        x = pdfData[0]
        pdfA = pdfData[1]
        pdfB = pdfData[2]
        
        x = [(j - Lx/(2*len(x)) - Lx/2) for j in x]
        pdfA = [j*fac for j in pdfA]
        pdfB = [j*fac for j in pdfB]
        pdfALL = [(pdfA[j] + pdfB[j]) for j in range(len(x))]
        
        phiALL = []
        for j in range(len(x)):
            if(round(pdfALL[j],2) == 0): phiALL.append(0)
            else: phiALL.append((pdfA[j] - pdfB[j])/pdfALL[j])
                    
        ax[0].plot(x, pdfA, "b-", marker = markers[i], lw = lw, ms = ms, 
                      label = "{val}".format(val = lamb[nr_list.index(plot_list[i])]))
        ax[1].plot(x, pdfB, "r-", marker = markers[i], lw = lw, ms = ms, 
                      label = "{val}".format(val=lamb[nr_list.index(plot_list[i])]))
        ax[2].plot(x, pdfALL, "-", marker = markers[i], lw = lw, ms = ms, color = colors[i], 
                      label = "{val}".format(val=lamb[nr_list.index(plot_list[i])]))
        
        phi_avg = [0.5*(phiALL[j] + phiALL[len(phiALL)-1-j]) for j in range(int(len(phiALL)/2))]
        phi_avg.reverse()
        
        ax2[0].plot(x[int(len(x)/2):len(x)], phi_avg, "-", marker = markers[i], lw = lw, ms = ms, color = colors[i],
                      label = "{val}".format(val=lamb[nr_list.index(plot_list[i])]))
    
    """
        for j in range(len(phi_avg)):
            print(x[int(len(x)/2)+j], phi_avg[j])
        print()
    
        ax2[0].axvline(x=34,ymin=0,ymax=1,lw=0.8)
        ax2[0].axvline(x=38,ymin=0,ymax=1,lw=0.8)
    """
    ax2[1].plot(lamb, width, "o-", lw = lw + 0.5, ms = ms + 1)
    
    
    for i in range(3):
        ax[i].set_xlim((-Lx/2, Lx/2))
        ax[i].set_xlabel(r"$x / \sigma$", fontsize = fontsize)
        ax[i].legend(title = r"$\lambda$", fontsize = fontsize - 4)
    
    ax[0].set_ylim((-0.01, 1.3))
    ax[0].set_ylabel(r"$\rho_A(x)~\sigma^3$", fontsize = fontsize)
    
    ax[1].set_ylim((-0.01, 1.3))
    ax[1].set_ylabel(r"$\rho_B(x)~\sigma^3$", fontsize = fontsize)
    
    ax[2].set_ylim((-0.01, 1.3))
    ax[2].set_ylabel(r"$\rho(x)~\sigma^3$", fontsize = fontsize)
    
    ax2[0].set_xlim((-Lx/2*0, Lx/2))
    ax2[0].set_xlabel(r"$x / \sigma$", fontsize = fontsize)
    ax2[0].legend(title = r"$\lambda$", fontsize = fontsize - 3)
    ax2[0].set_ylim((-1.05, 1.05))
    ax2[0].set_ylabel(r"$\phi(x)$", fontsize = fontsize)
    
    ax2[1].set_xlabel(r"$\lambda$", fontsize = fontsize)
    ax2[1].set_ylabel(r"$w$", fontsize = fontsize)
    ax2[1].set_yticks(range(0,18,2))

    fig.set_figheight(3)
    fig.set_figwidth(9)
    fig.tight_layout(pad = 0.5)
    
    fig2.set_figheight(3)
    fig2.set_figwidth(9)
    fig2.tight_layout(pad = 0.5)
    
    #fig.savefig(figpath, dpi = 600, bbox_inches = "tight")  
    fig2.savefig(figpath, dpi = 600, bbox_inches = "tight")