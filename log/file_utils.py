import os
from ast import literal_eval as liteval
import numpy as np

def readLogParam(entry, fpath):
    fobj = open(fpath, mode = "r", encoding = "utf-8")
    for line in fobj:
        line = line.removesuffix("\n").split(sep = " ")
        if(line[0] == entry):
            return(liteval(line[1]))
            break
    fobj.close()

def makeFilePath(nr, type1, type2, extras):
    fpath = os.getcwd().removesuffix("/ld_analysis/log") + \
            "/LD/{type1}/Data{nr}".format(type1 = type1, nr = nr)
    if(extras[0] == 'log'): fpath += "/log.dat"
    elif(extras[0] == 'dump'): 
        fpath += "/relaxa_{type2}/dump_{type2}.dat".format(type2 = type2)
    elif(extras[0] == 'thermo'): fpath += "/thermo.dat"
    elif(extras[0] == 'fig'): fpath = fpath.removesuffix("/{type1}/Data{nr}".format(type1 = type1, nr = nr)) +\
                                      "/imgs/msd_{type1}_{nr}.png".format(type1=type1, nr = extras[1])
    elif(extras[0] == 'msd' and type1 == 'lmp'): fpath += "/msd.dat"
    elif(extras[0] == 'msd' and type1 == 'cpp'): fpath += "/relaxa_{type2}/msd.dat".format(type2 = type2) 
    return(fpath)
    
def line2array(line):
    line = line.removesuffix("\n").split(sep = " ")
    line2 = []
    for element in line:
        if(element != ''): line2.append(element)
    if(line2 == []): line2.append('0')
    return(line2)

def array2line(array):
    line = ""
    for i in range(len(array)):
        line += str(array[i]) + " "
    return(line.removesuffix(" "))

def readData(fpath, nDataSets = 1):
    fobj = open(fpath, mode = 'r', encoding = "utf-8")
    line = fobj.readline().removesuffix("\n").split(sep = " ")
    nCols = len(line)
    fields = [line]
    
    if(nDataSets == 1):
        Cols = [[] for i in range(nCols)]
    
        for line in fobj:
            line = line.removesuffix("\n").split(sep = " ")
            if(line[0].isalpha() == False):
                for i in range(nCols):
                    if(line[i].isalpha() == True): 
                        Cols[i].append(line[i])
                    else:
                        Cols[i].append(liteval(line[i]))
            else:
                continue
        
    else:
        Cols = [[] for i in range(nCols*nDataSets)]
        
        dataCtr = 0
        for line in fobj:
            line = line.removesuffix("\n").split(sep = " ")
            if(line[0].isalpha() == False):
                for i in range(nCols):
                    if(line[i].isalpha() == True): 
                        Cols[nCols*dataCtr + i].append(line[i])
                    else:
                        Cols[nCols*dataCtr + i].append(liteval(line[i]))
            else:
                fields.append(line)
                dataCtr += 1
                continue
    
    fobj.close()
    return([Cols, fields])

def readData2D(fpath, nDatasets, nRows, nCols):
    
    fobj = open(fpath, mode = 'r', encoding = "utf-8")
    data = np.ndarray((nDatasets, nRows, nCols))
    
    for i in range(nDatasets):
        line = fobj.readline()
        
        for j in range(nRows):
            line = fobj.readline().removesuffix(" \n").split(" ")
            
            for k in range(len(line)):
                data[i][j][k] = liteval(line[k])
    
    fobj.close()
    return(data)
    
def copyList(a, lt = None, gt = None, leq = None, geq = None):
    temp = []
    
    if(lt != None):
        for i in a:
            if(i < lt):
                temp.append(i)

    elif(gt != None):
        for i in a:
            if(i > gt):
                temp.append(i)

    elif(leq != None):
        for i in a:
            if(i <= leq):
                temp.append(i)

    elif(geq != None):
        for i in a:
            if(i >= geq):
                temp.append(i)
    else:
        temp = a

    return(temp)