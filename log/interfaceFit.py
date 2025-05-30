print("\033[J\033[H",end='') # Clear screen

import math as m
import scipy
import numpy as np
import matplotlib as mpl
import file_utils

mpl.pyplot.rc('text', usetex=True)

def getGibbsPlane(x, pdf):
    
    if(pdf[0] < pdf[len(pdf) - 1]):
        phiG = pdf[0]
        phiL = pdf[len(pdf) - 1]
    else:
        phiG = pdf[len(pdf) - 1]
        phiL = pdf[0]
        
    int1 = []
    int2 = []
    for pos in range(0,len(x)):
        x1 = x[0:pos+1]
        x2 = x[pos:len(x)]
        
        if(pdf[0] < pdf[len(pdf) - 1]):
            y1 = [(i - phiG) for i in pdf[0:pos+1]]
            y2 = [(phiL - i) for i in pdf[pos:len(x)]]
        else:
            y1 = [(phiL - i) for i in pdf[0:pos+1]]
            y2 = [(i - phiG) for i in pdf[pos:len(x)]]
            
        int1.append(np.trapz(y1,x1))
        int2.append(np.trapz(y2,x2))
            
    for i in range(len(x)):
        if(int1[i] >= int2[i]):
            xplus = x[i]
            xminus = x[i-1]
            break
    
    return(0.5*(xplus + xminus))

def phi(x, A, B, C):
    return(A + B*np.tanh(x/C))

def fitPhiProfile(x, pdf, xG):
    x = [(i-xG) for i in x]
    
    popt, pcov = scipy.optimize.curve_fit(phi, x, pdf)
    return(popt)

fac = 1
markers = ["o", "s", "^", "*", "D"]
lw = 0.5
ms = 3
fs = "none"
fontsize = 14

option = "tavg"
if(option == "frames"):
    # Params for frames-based data
    nr = 21
    Lx = 200
    Pe = "1.0"
    lamb = " 10^{-3}" 
    pdfSamples = [11e3, 12e3, 13e3, 14e3, 15e3]
    nSamples = 5
    sampleNum = [0, 4]
    
    fpath = "../../LD/LD-cpp/Data{nr}/pdf_frames.dat".format(nr = nr)
    figpath = "../../LD/LD-cpp/imgs/pdf_frames_all_{nr}.png".format(nr = nr)
    pdfData = file_utils.readData(fpath, nSamples)
    
    fig, ax = mpl.pyplot.subplots()
    for i in sampleNum:
        x = pdfData[3*i + 0]
        pdfA = pdfData[3*i + 1]
        pdfB = pdfData[3*i + 2]
        
        x = [(j-Lx/(2*len(x))-Lx/2) for j in x]
        pdfALL = []
        for j in range(len(x)):
            pdfALL.append((pdfA[j] + pdfB[j])*fac)
        
        ax.plot(x, pdfALL, "-o", lw = lw, ms = ms, label = str(pdfSamples[i]))

    ax.set(xlim = (-Lx/2, Lx/2), ylim = (-0.01, 1.0), aspect = Lx)
    ax.set_xlabel(r"$x ~ (\sigma)$", fontsize = fontsize - 4)
    ax.set_ylabel(r"$\phi ~ (x)$", fontsize = fontsize - 4)
    ax.legend(title = "time", fontsize = fontsize - 8)
    fig.suptitle(r"Simulation Parameters : $Pe_s = {Pe},~ \lambda = {lamb}$".format(Pe=Pe,lamb=lamb), 
                 fontsize = fontsize - 4)

if(option == "tavg"):
    # Params for time-averaged data
    nr_list = [32]
    Pe = [5]
    lamb = ["1"]
    Lx = 200.0
    
    fig, ax = mpl.pyplot.subplots()
    fig2, ax2 = mpl.pyplot.subplots()
    for i in range(len(nr_list)):
        fpath = "../../LD/LD-cpp/Data{nr}/pdf_tavg.dat".format(nr = nr_list[i])
        figpath = "../../LD/LD-cpp/imgs/pdf_fit1.png"
        pdfData = file_utils.readData(fpath, 1)
        
        x = pdfData[0]
        pdfA = pdfData[1]
        pdfB = pdfData[2]
        
        x = [(j-Lx/(2*len(x))-Lx/2) for j in x]
        
        pdfALL = []
        pdfALLavg = []
        for j in range(len(x)):
            pdfALL.append((pdfA[j] + pdfB[j])*fac)
        for j in range(int(len(x)/2)):
            pdfALLavg.append(0.5*(pdfALL[j] + pdfALL[len(x) - 1 - j]))
        
        mid = int(len(pdfALL)/2)
        ax.plot(x, pdfALL, "-o", lw = lw, ms = ms, fillstyle = fs, label = str(Pe[i]))
        ax2.plot(x[0:mid], pdfALLavg, "-o", lw = lw, ms = ms, fillstyle = fs, 
                 label = str(Pe[i]), color = "crimson")
    
        xG1 = getGibbsPlane(x[0:mid], pdfALL[0:mid])
        xG2 = getGibbsPlane(x[mid:len(x)], pdfALL[mid:len(x)])
        
        C1 = fitPhiProfile(x[0:mid], pdfALL[0:mid],xG1)
        C2 = fitPhiProfile(x[mid:len(x)], pdfALL[mid:len(x)], xG2)
        
        C = [0.5*(abs(C1[i])+abs(C2[i])) for i in range(len(C1))]
        print("Pe = {Pe}, rhoG = {rhoG}, rhoL = {rhoL}, D = {D}".format(Pe = Pe[i], 
                    rhoG = C[0]+C[1], rhoL = C[0]-C[1], D = m.sqrt(m.pi)*C[2]))
        
        pdfFit = [C[0]+C[1]*np.tanh((i-xG1)/C[2]) for i in x[0:mid]]
        ax2.plot(x[0:mid], pdfFit, "-o", lw = lw, ms = ms, fillstyle = fs)
        
        """
        SStot = 0
        SSres = 0
        pdfTestFit = pdfALL[0:mid]
        for i in range(len(pdfFit)):
            SSres += (pdfFit[i] - pdfTestFit[i])**2
            SStot += (pdfTestFit[i] - sum(pdfTestFit)/len(pdfTestFit))**2
        print(1-SSres/SStot)
        """
    
    ax.set(xlim = (-Lx/2, Lx/2), ylim = (-0.01, 1.8), aspect = Lx/2)
    ax.set_xlabel(r"$x~(\sigma)$", fontsize = fontsize - 4)
    ax.set_ylabel(r"$\rho(x) \sigma^3$", fontsize = fontsize - 4)
    ax.legend(title = r"$Pe_s$", fontsize = fontsize - 8)
    fig.suptitle(r"$\lambda = {lamb}$".format(lamb = lamb[0]))
    
    ax2.set(xlim = (-Lx/2,0), ylim = (-0.01, 1.3), aspect = Lx/2)
    ax2.set_xlabel(r"$x/\sigma$", fontsize = fontsize - 2)
    ax2.set_ylabel(r"$\rho(x) \sigma^3$", fontsize = fontsize - 2)
    #ax2.legend(title = r"$Pe_s$", fontsize = fontsize - 3)
    fig.set_figheight(5)
    fig.set_figwidth(7)
    
    fig2.savefig(figpath, dpi = 600, bbox_inches = "tight")