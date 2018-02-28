
import os
import sys
import math
import numpy as np
from bisect import bisect_left
import matplotlib.pyplot as plt
from scipy import interpolate

def calculatePhs(sigma, wavelength, inclination, thickness, threshold, ranges, stdevkev=50):
    EOfrag = [1420, 780, 1750, 960] # KeV, energy max fragments Alpha 94, Li 94, aplha 6, Li 6
    if threshold < 2:
        threshold = 2
        print("Threshold set to 2")
    EOfrag = [1420, 780, 1750, 960]

    #pf = '/Users/francescopiscitelli/Documents/MATLAB/000_data/SRIM_real_10B4C_layer_Linkoping';
    #fprintf('\n Linkoping coatings used as default for stopping power! \n')

    Efrag = np.loadtxt(fname=os.path.dirname(os.path.abspath(__file__)) + "/data/B10/10B4C224/Erim_fragments.txt", unpack=True)

    L = np.arange(0.001, 8, 0.05) #distance traveled in boron

    limitRange = ranges
    p=[[],[]]
    energ=[[],[],[],[]]
    defrag=[[],[],[],[]]

    dl = np.diff(Efrag[4])
    for d in range(0, len(defrag)):
        defrag[d] = np.abs(np.diff(Efrag[d])/dl)
    #Defrag = [Defrag ; Defrag(end, 1: 4)];
    pe = [[],[],[],[]]
    peT = [[],[],[],[]]
    for k in range(0,len(L)):
        ptemp = PHSprob(L[k],thickness,sigma[0])
        p[0].append(ptemp[0]) #BS
        p[1].append(ptemp[1]) #T
        for j in range (0,4):
            if L[k] > 0.05 and L[k] < ranges[j]: # um(SRIM lower limit)
                index = find_nearest(Efrag[4],L[k])
                energ[j].append(Efrag[j][index])
                pe[j].append(p[0][k]/defrag[j][index])
                peT[j].append(p[1][k] / defrag[j][index])
            elif L[k] <= 0.05:
                energ[j].append(EOfrag[j])
                pe[j].append(p[0][k]/defrag[j][101])
                peT[j].append(p[1][k]/defrag[j][101])
            else:
                energ[j].append(0)
                pe[j].append(0)
                peT[j].append(0)
    prob = [np.array(pe[0])*94/100,np.array(pe[1])*94/100,np.array(pe[2])*6/100,np.array(pe[3])*6/100]
    probT = [np.array(peT[0])*94/100,np.array(peT[1])*94/100,np.array(peT[2])*6/100,np.array(peT[3])*6/100]

    #prob = [0 0 0 0;prob];
    #probT = [0 0 0 0;probT];
    energ = [EOfrag,energ]


   # DE = np.diff(energ,1,1)
    #DE = [0 0 0 0;DE];
  #  effic = np.diag(prob*DE)



    return energ, pe, peT

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx


def PHSprob(l, thickness, sigma):
    p=[[],[]]
    if l >= 0 and l <= thickness:
        #BS
        p[0] = (1./2.)*(1./(l**2.))*( (1./sigma) - (( l+(1./sigma) )*np.exp([-l*sigma])) )[0]
        #TS
        p[1] = (1/2)*(1/(l**2))*(np.exp(-thickness*sigma))*( ((np.exp(l*sigma))*l) - (1/sigma)*((np.exp(l*sigma))-1) )
    elif l>thickness:
        #BS
        p[0] = (1 / 2) * (1 / (l ** 2)) * ((1 / sigma) - ((thickness + (1 / sigma)) * np.exp(-thickness * sigma)))
        #Transmission
        p[1] = (1 / 2) * (1 / (l ** 2)) * (thickness - (1 / sigma) * (1 - np.exp(-thickness * sigma)))
    else:
        p[0]=0
        p[1]=0

    return p
