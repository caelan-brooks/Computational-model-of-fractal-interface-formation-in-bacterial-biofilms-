from random import seed
from random import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import imageio
from PIL import Image
from scipy.optimize import leastsq
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap

# Adapted from fractal-dimension.py by Nicolas P. Rougier
# github url: https://gist.github.com/rougier/e5eafc276a4e54f516ed5559df4242c0

def fractal_dimension(Z, threshold=0.9):
    # for 2d image
    assert(len(Z.shape) == 2)

    def func(x):
        return 2*(x-3)**2+1
    def boxcount(Z, k):
        S = np.add.reduceat(
            np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0),
                               np.arange(0, Z.shape[1], k), axis=1)
        #print(S)

        # We count non-empty (0) and non-full/full boxes (k*k)
        return len(np.where((S > 0) & (S <= k*k))[0])
    def boxcount2(Z, k):
        xback=[]
        step=k
        xleft=Z.shape[0]
        while xleft>0:
            xleft=xleft-step
            if xleft<=0:
                xleft=0
            xback.append(xleft)
        yback=[]
        step=k
        yleft=Z.shape[1]
        while yleft>0:
            yleft=yleft-step
            if yleft<=0:
                yleft=0
            yback.append(yleft)
        S = np.add.reduceat(
            np.add.reduceat(Z, np.sort(xback), axis=0),
                               np.sort(yback), axis=1)
        #print(S)

        # We count non-empty (0) and non-full/full boxes (k*k)
        return len(np.where((S > 0) & (S <= k*k))[0])
    def boxcount3(Z, k):
        yback=[]
        step=k
        yleft=Z.shape[1]
        while yleft>0:
            yleft=yleft-step
            if yleft<=0:
                yleft=0
            yback.append(yleft)
        S = np.add.reduceat(
            np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0),
                               np.sort(yback), axis=1)
        #print(S)

        # We count non-empty (0) and non-full/full boxes (k*k)
        return len(np.where((S > 0) & (S <= k*k))[0])
    def boxcount4(Z, k):
        xback=[]
        step=k
        xleft=Z.shape[0]
        while xleft>0:
            xleft=xleft-step
            if xleft<=0:
                xleft=0
            xback.append(xleft)
        S = np.add.reduceat(
            np.add.reduceat(Z, np.sort(xback), axis=0),
                               np.arange(0, Z.shape[1], k), axis=1)
        #print(S)

        # We count non-empty (0) and non-full/full boxes (k*k)
        return len(np.where((S > 0) & (S <= k*k))[0])


    # Transform Z into a binary array
    #Z = (Z > threshold)

    # Minimal dimension of image
    p = min(Z.shape)

    # Greatest power of 2 less than or equal to p
    n = 2**np.floor(np.log(p)/np.log(2))

    # Extract the exponent
    n = int(np.log(n)/np.log(2))
    
    #backwards indices calculation
    

    # Build successive box sizes (from 2**n down to 2**1)
    sizes = np.unique(np.floor(2**np.linspace(0,n,num=3*n, endpoint=True)))
    sizecollect=[]

    # Actual box counting with decreasing size
    counts = []
    #print(sizes)
    for size in sizes:
        counts.append(boxcount(Z, int(size)))
        sizecollect.append(size)
        counts.append(boxcount2(Z,int(size)))
        sizecollect.append(size)
        counts.append(boxcount3(Z,int(size)))
        sizecollect.append(size)
        counts.append(boxcount4(Z,int(size)))
        sizecollect.append(size)
        

    # Fit the successive log(sizes) with log (counts)
    #print(sizecollect,counts)
    coeffs = np.polyfit(np.log(sizecollect), np.log(counts), 1)
    #print(sizecollect,counts)
    #plt.figure()
    #plt.plot(np.log(sizecollect),np.log(counts),"g*")
    #xaxisval=np.linspace(min(np.log(sizecollect)),max(np.log(sizecollect)),10)
    #ydot= coeffs[0]*xaxisval+coeffs[1]
    #plt.plot(xaxisval,ydot)
    #plt.show
    return -coeffs[0]
