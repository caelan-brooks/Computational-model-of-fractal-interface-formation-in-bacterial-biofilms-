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
def fixedinitialcondition_1(massratio=1):
    #massratio is the mass ratio of matrix producing cell to motile cell
    #fixed initial condition
    radius=7 #radius of proposed initial biofilm
    ylevel=radius
    xpo=[] #record of individual cell's horizontal position
    ypo=[] #record of individual cell's vertical position
    growthtime=[] #record of individual cell's growthtime
    cellmass=[] #record of individual cell's mass
    celltype=[] #record of individual cell's celltype -- matrix producing cell(type-1) or motile cell(type-0.5)
    norfacx=1  #growthtime for matrix producing cell
    norfacy=1  #growthtime for motile cell
    norfacmx=1*massratio #mass for matrix producing cell
    norfacmy=1 #mass for motile cell
    while ylevel>=-radius:
        xlevel=0
        xmotile=0
        if ylevel==7 or ylevel==-7:
            xmotile=0
            xlevel=2
        elif ylevel==6 or ylevel==-6:
            xmotile=0
            xlevel=4
        elif ylevel==5 or ylevel==-5:
            xmotile=2
            xlevel=5
        elif ylevel==4 or ylevel==-4:
            xmotile=3
            xlevel=6
        elif ylevel==3 or ylevel==-3:
            xmotile=4
            xlevel=6
        elif ylevel in range(-2,3):
            xmotile=5
            xlevel=7
        po=xlevel
        while po>=-xlevel:
            xpo.append(po)
            ypo.append(ylevel)
            if xmotile!=0:
                if po in range(-xmotile,xmotile+1):
                    celltype.append(0.5)
                    s=1
                    while s>0:
                        o=np.random.normal(norfacmy, 0.1)
                        if o>0:
                            cellmass.append(o)
                            break
                        else:
                            s=s+1
                    f=1
                    while f>0:
                        o=np.random.normal(norfacy, 0.1)
                        if o>0:
                            growthtime.append(o)
                            break
                        else:
                            f=f+1
                else:
                    celltype.append(1)
                    s=1
                    while s>0:
                        o=np.random.normal(norfacmx, 0.1)
                        if o>0:
                            cellmass.append(o)
                            break
                        else:
                            s=s+1
                    f=1
                    while f>0:
                        o=np.random.normal(norfacx, 0.1)
                        if o>0:
                            growthtime.append(o)
                            break
                        else:
                            f=f+1

            else:
                celltype.append(1)
                s=1
                while s>0:
                    o=np.random.normal(norfacmx, 0.1)
                    if o>0:
                        cellmass.append(o)
                        break
                    else:
                        s=s+1
                f=1
                while f>0:
                    o=np.random.normal(norfacx, 0.1)
                    if o>0:
                        growthtime.append(o)
                        break
                    else:
                        f=f+1
            po=po-1

        ylevel=ylevel-1
    initialmotilearea=sum(np.array(celltype)<0.7)
    initialmatrixarea=sum(np.array(celltype)>0.7)
    ##ploting initial condition
    #data = pd.DataFrame(data={'x':ypo, 'y':xpo, 'z':celltype})
    #data = data.pivot(index='x', columns='y', values='z')
    #fig, ax = plt.subplots(1, 1, figsize=[4, 4])
    #cMap = ListedColormap(['#00BF00','#BF00BF'])
    #sns.heatmap(data,cmap=cMap,linewidths=.5, ax=ax,cbar=False)
    #ax.axis(False)
    ##end of plot
    return xpo,ypo,growthtime,cellmass,celltype,initialmotilearea,initialmatrixarea
def fixedinitialcondition_2(massratio=1):
    #massratio is the mass ratio of matrix producing cell to motile cell
    #fixed initial condition
    radius=9 #radius of proposed initial biofilm
    ylevel=radius
    xpo=[] #record of individual cell's horizontal position
    ypo=[] #record of individual cell's vertical position
    growthtime=[] #record of individual cell's growthtime
    cellmass=[] #record of individual cell's mass
    celltype=[] #record of individual cell's celltype -- matrix producing cell or motile cell
    norfacx=1  #growthtime for matrix producing cell
    norfacy=1  #growthtime for motile cell
    norfacmx=1*massratio #mass for matrix producing cell
    norfacmy=1 #mass for motile cell
        
    while ylevel>=-radius:
        xlevel=0
        xmotile=0
        if ylevel==8 or ylevel==-8:
            xmotile=0
            xlevel=3
        elif ylevel==9 or ylevel==-9:
            xmotile=0
            xlevel=2
        elif ylevel==7 or ylevel==-7:
            xmotile=0
            xlevel=5
        elif ylevel==6 or ylevel==-6:
            xmotile=2
            xlevel=6
        elif ylevel==5 or ylevel==-5:
            xmotile=3
            xlevel=7
        elif ylevel==4 or ylevel==-4:
            xmotile=4
            xlevel=7
        elif ylevel==3 or ylevel==-3:
            xmotile=4
            xlevel=8
            
        elif ylevel in range(-2,3):
            xmotile=5
            xlevel=8
        po=xlevel
        while po>=-xlevel:
            xpo.append(po)
            ypo.append(ylevel)
            if xmotile!=0:
                if po in range(-xmotile,xmotile+1):
                    celltype.append(0.5)
                    s=1
                    while s>0:
                        o=np.random.normal(norfacmy, 0.1)
                        if o>0:
                            cellmass.append(o)
                            break
                        else:
                            s=s+1
                    f=1
                    while f>0:
                        o=np.random.normal(norfacy, 0.1)
                        if o>0:
                            growthtime.append(o)
                            break
                        else:
                            f=f+1
                else:
                    celltype.append(1)
                    s=1
                    while s>0:
                        o=np.random.normal(norfacmx, 0.1)
                        if o>0:
                            cellmass.append(o)
                            break
                        else:
                            s=s+1
                    f=1
                    while f>0:
                        o=np.random.normal(norfacx, 0.1)
                        if o>0:
                            growthtime.append(o)
                            break
                        else:
                            f=f+1

            else:
                celltype.append(1)
                s=1
                while s>0:
                    o=np.random.normal(norfacmx, 0.1)
                    if o>0:
                        cellmass.append(o)
                        break
                    else:
                        s=s+1
                f=1
                while f>0:
                    o=np.random.normal(norfacx, 0.1)
                    if o>0:
                        growthtime.append(o)
                        break
                    else:
                        f=f+1
            po=po-1

        ylevel=ylevel-1
    initialmotilearea=sum(np.array(celltype)<0.7)
    initialmatrixarea=sum(np.array(celltype)>0.7)
    ##ploting initial condition
    #data = pd.DataFrame(data={'x':ypo, 'y':xpo, 'z':celltype})
    #data = data.pivot(index='x', columns='y', values='z')
    #fig, ax = plt.subplots(1, 1, figsize=[4, 4])
    #cMap = ListedColormap(['#00BF00','#BF00BF'])
    #sns.heatmap(data,cmap=cMap,linewidths=.5, ax=ax,cbar=False)
    #ax.axis(False)
    ##end of plot
    return xpo,ypo,growthtime,cellmass,celltype,initialmotilearea,initialmatrixarea

def fixedinitialcondition_3(massratio=1):
    #massratio is the mass ratio of matrix producing cell to motile cell
    #fixed initial condition
    radius=6
    ylevel=radius
    xpo=[] #record of individual cell's horizontal position
    ypo=[] #record of individual cell's vertical position
    growthtime=[] #record of individual cell's growthtime
    cellmass=[] #record of individual cell's mass
    celltype=[] #record of individual cell's celltype -- matrix producing cell or motile cell
    norfacx=1  #growthtime for matrix producing cell
    norfacy=1  #growthtime for motile cell
    norfacmx=1*massratio #mass for matrix producing cell
    norfacmy=1 #mass for motile cell
    
    while ylevel>=-radius:
        xlevel=0
        xmotile=0
        
        if ylevel==6 or ylevel==-6:
            xmotile=0
            xlevel=2
        elif ylevel==5 or ylevel==-5:
            xmotile=0
            xlevel=3
        elif ylevel==4 or ylevel==-4:
            xmotile=1
            xlevel=4
        elif ylevel==3 or ylevel==-3:
            xmotile=2
            xlevel=5
        elif ylevel==2 or ylevel==-2:
            xmotile=3
            xlevel=6
        elif ylevel in range(-1,2):
            xmotile=4
            xlevel=6
        po=xlevel
        while po>=-xlevel:
            xpo.append(po)
            ypo.append(ylevel)
            if xmotile!=0:
                if po in range(-xmotile,xmotile+1):
                    celltype.append(0.5)
                    s=1
                    while s>0:
                        o=np.random.normal(norfacmy, 0.1)
                        if o>0:
                            cellmass.append(o)
                            break
                        else:
                            s=s+1
                    f=1
                    while f>0:
                        o=np.random.normal(norfacy, 0.1)
                        if o>0:
                            growthtime.append(o)
                            break
                        else:
                            f=f+1
                else:
                    celltype.append(1)
                    s=1
                    while s>0:
                        o=np.random.normal(norfacmx, 0.1)
                        if o>0:
                            cellmass.append(o)
                            break
                        else:
                            s=s+1
                    f=1
                    while f>0:
                        o=np.random.normal(norfacx, 0.1)
                        if o>0:
                            growthtime.append(o)
                            break
                        else:
                            f=f+1

            else:
                celltype.append(1)
                s=1
                while s>0:
                    o=np.random.normal(norfacmx, 0.1)
                    if o>0:
                        cellmass.append(o)
                        break
                    else:
                        s=s+1
                f=1
                while f>0:
                    o=np.random.normal(norfacx, 0.1)
                    if o>0:
                        growthtime.append(o)
                        break
                    else:
                        f=f+1
            po=po-1

        ylevel=ylevel-1
    initialmotilearea=sum(np.array(celltype)<0.7)
    initialmatrixarea=sum(np.array(celltype)>0.7)
    ##ploting initial condition
    #data = pd.DataFrame(data={'x':ypo, 'y':xpo, 'z':celltype})
    #data = data.pivot(index='x', columns='y', values='z')
    #fig, ax = plt.subplots(1, 1, figsize=[4, 4])
    #cMap = ListedColormap(['#00BF00','#BF00BF'])
    #sns.heatmap(data,cmap=cMap,linewidths=.5, ax=ax,cbar=False)
    #ax.axis(False)
    ##end of plot
    return xpo,ypo,growthtime,cellmass,celltype,initialmotilearea,initialmatrixarea



def fixedinitialcondition_4(massratio=1):
    #massratio is the mass ratio of matrix producing cell to motile cell
    #fixed initial condition
    radius=6
    ylevel=radius 
    xpo=[] #record of individual cell's horizontal position
    ypo=[] #record of individual cell's vertical position
    growthtime=[] #record of individual cell's growthtime
    cellmass=[] #record of individual cell's mass
    celltype=[] #record of individual cell's celltype -- matrix producing cell or motile cell
    norfacx=1  #growthtime for matrix producing cell
    norfacy=1  #growthtime for motile cell
    norfacmx=1*massratio #mass for matrix producing cell
    norfacmy=1 #mass for motile cell
    
    while ylevel>=-radius:
        xlevel=0
        xmotile=0
        
        if ylevel==6 or ylevel==-6:
            xmotile=0
            xlevel=6
        elif ylevel==5 or ylevel==-5:
            xmotile=0
            xlevel=6
        elif ylevel==4 or ylevel==-4:
            xmotile=4
            xlevel=6
        elif ylevel==3 or ylevel==-3:
            xmotile=4
            xlevel=6
        elif ylevel==2 or ylevel==-2:
            xmotile=4
            xlevel=6
        elif ylevel in range(-1,2):
            xmotile=4
            xlevel=6
        po=xlevel
        while po>=-xlevel:
            xpo.append(po)
            ypo.append(ylevel)
            if xmotile!=0:
                if po in range(-xmotile,xmotile+1):
                    celltype.append(0.5)
                    s=1
                    while s>0:
                        o=np.random.normal(norfacmy, 0.1)
                        if o>0:
                            cellmass.append(o)
                            break
                        else:
                            s=s+1
                    f=1
                    while f>0:
                        o=np.random.normal(norfacy, 0.1)
                        if o>0:
                            growthtime.append(o)
                            break
                        else:
                            f=f+1
                else:
                    celltype.append(1)
                    s=1
                    while s>0:
                        o=np.random.normal(norfacmx, 0.1)
                        if o>0:
                            cellmass.append(o)
                            break
                        else:
                            s=s+1
                    f=1
                    while f>0:
                        o=np.random.normal(norfacx, 0.1)
                        if o>0:
                            growthtime.append(o)
                            break
                        else:
                            f=f+1

            else:
                celltype.append(1)
                s=1
                while s>0:
                    o=np.random.normal(norfacmx, 0.1)
                    if o>0:
                        cellmass.append(o)
                        break
                    else:
                        s=s+1
                f=1
                while f>0:
                    o=np.random.normal(norfacx, 0.1)
                    if o>0:
                        growthtime.append(o)
                        break
                    else:
                        f=f+1
            po=po-1

        ylevel=ylevel-1
    initialmotilearea=sum(np.array(celltype)<0.7)
    initialmatrixarea=sum(np.array(celltype)>0.7)
    ##ploting initial condition
    #data = pd.DataFrame(data={'x':ypo, 'y':xpo, 'z':celltype})
    #data = data.pivot(index='x', columns='y', values='z')
    #fig, ax = plt.subplots(1, 1, figsize=[4, 4])
    #cMap = ListedColormap(['#00BF00','#BF00BF'])
    #sns.heatmap(data,cmap=cMap,linewidths=.5, ax=ax,cbar=False)
    #ax.axis(False)
    ##end of plot
    return xpo,ypo,growthtime,cellmass,celltype,initialmotilearea,initialmatrixarea

def fixedinitialcondition_5(massratio=1):
    #massratio is the mass ratio of matrix producing cell to motile cell
    #fixed initial condition
    radius=5
    ylevel=radius
    xpo=[] #record of individual cell's horizontal position
    ypo=[] #record of individual cell's vertical position
    growthtime=[] #record of individual cell's growthtime
    cellmass=[] #record of individual cell's mass
    celltype=[] #record of individual cell's celltype -- matrix producing cell or motile cell
    norfacx=1  #growthtime for matrix producing cell
    norfacy=1  #growthtime for motile cell
    norfacmx=1*massratio #mass for matrix producing cell
    norfacmy=1 #mass for motile cell
    
    while ylevel>=-radius:
        xlevel=0
        xmotile=0
        
        #if ylevel==6 or ylevel==-6:
         #   xmotile=0
          #  xlevel=6
        if ylevel==5 or ylevel==-5:
            xmotile=0
            xlevel=5
        elif ylevel==4 or ylevel==-4:
            xmotile=0
            xlevel=5
        elif ylevel==3 or ylevel==-3:
            xmotile=0
            xlevel=5
        elif ylevel==2 or ylevel==-2:
            xmotile=2
            xlevel=5
        elif ylevel in range(-1,2):
            xmotile=2
            xlevel=5
        po=xlevel
        while po>=-xlevel:
            xpo.append(po)
            ypo.append(ylevel)
            if xmotile!=0:
                if po in range(-xmotile,xmotile+1):
                    celltype.append(0.5)
                    s=1
                    while s>0:
                        o=np.random.normal(norfacmy, 0.1)
                        if o>0:
                            cellmass.append(o)
                            break
                        else:
                            s=s+1
                    f=1
                    while f>0:
                        o=np.random.normal(norfacy, 0.1)
                        if o>0:
                            growthtime.append(o)
                            break
                        else:
                            f=f+1
                else:
                    celltype.append(1)
                    s=1
                    while s>0:
                        o=np.random.normal(norfacmx, 0.1)
                        if o>0:
                            cellmass.append(o)
                            break
                        else:
                            s=s+1
                    f=1
                    while f>0:
                        o=np.random.normal(norfacx, 0.1)
                        if o>0:
                            growthtime.append(o)
                            break
                        else:
                            f=f+1

            else:
                celltype.append(1)
                s=1
                while s>0:
                    o=np.random.normal(norfacmx, 0.1)
                    if o>0:
                        cellmass.append(o)
                        break
                    else:
                        s=s+1
                f=1
                while f>0:
                    o=np.random.normal(norfacx, 0.1)
                    if o>0:
                        growthtime.append(o)
                        break
                    else:
                        f=f+1
            po=po-1

        ylevel=ylevel-1
    initialmotilearea=sum(np.array(celltype)<0.7)
    initialmatrixarea=sum(np.array(celltype)>0.7)
    ##ploting initial condition
    #data = pd.DataFrame(data={'x':ypo, 'y':xpo, 'z':celltype})
    #data = data.pivot(index='x', columns='y', values='z')
    #fig, ax = plt.subplots(1, 1, figsize=[4, 4])
    #cMap = ListedColormap(['#00BF00','#BF00BF'])
    #sns.heatmap(data,cmap=cMap,linewidths=.5, ax=ax,cbar=False)
    #ax.axis(False)
    ##end of plot
    return xpo,ypo,growthtime,cellmass,celltype,initialmotilearea,initialmatrixarea


def fixedinitialcondition_6(massratio=1):
    #massratio is the mass ratio of matrix producing cell to motile cell
    #fixed initial condition
    radius=5
    ylevel=radius
    xpo=[] #record of individual cell's horizontal position
    ypo=[] #record of individual cell's vertical position
    growthtime=[] #record of individual cell's growthtime
    cellmass=[] #record of individual cell's mass
    celltype=[] #record of individual cell's celltype -- matrix producing cell or motile cell
    norfacx=1  #growthtime for matrix producing cell
    norfacy=1  #growthtime for motile cell
    norfacmx=1*massratio #mass for matrix producing cell
    norfacmy=1 #mass for motile cell
    
    while ylevel>=-radius:
        xlevel=0
        xmotile=0
        
        #if ylevel==6 or ylevel==-6:
         #   xmotile=0
          #  xlevel=6
        if ylevel==5 or ylevel==-5:
            xmotile=0
            xlevel=5
        elif ylevel==4 or ylevel==-4:
            xmotile=0
            xlevel=5
        elif ylevel==3 or ylevel==-3:
            xmotile=2
            xlevel=5
        elif ylevel==2 or ylevel==-2:
            xmotile=3
            xlevel=5
        elif ylevel in range(-1,2):
            xmotile=3
            xlevel=5
        po=xlevel
        while po>=-xlevel:
            xpo.append(po)
            ypo.append(ylevel)
            if xmotile!=0:
                if po in range(-xmotile,xmotile+1):
                    celltype.append(0.5)
                    s=1
                    while s>0:
                        o=np.random.normal(norfacmy, 0.1)
                        if o>0:
                            cellmass.append(o)
                            break
                        else:
                            s=s+1
                    f=1
                    while f>0:
                        o=np.random.normal(norfacy, 0.1)
                        if o>0:
                            growthtime.append(o)
                            break
                        else:
                            f=f+1
                else:
                    celltype.append(1)
                    s=1
                    while s>0:
                        o=np.random.normal(norfacmx, 0.1)
                        if o>0:
                            cellmass.append(o)
                            break
                        else:
                            s=s+1
                    f=1
                    while f>0:
                        o=np.random.normal(norfacx, 0.1)
                        if o>0:
                            growthtime.append(o)
                            break
                        else:
                            f=f+1

            else:
                celltype.append(1)
                s=1
                while s>0:
                    o=np.random.normal(norfacmx, 0.1)
                    if o>0:
                        cellmass.append(o)
                        break
                    else:
                        s=s+1
                f=1
                while f>0:
                    o=np.random.normal(norfacx, 0.1)
                    if o>0:
                        growthtime.append(o)
                        break
                    else:
                        f=f+1
            po=po-1

        ylevel=ylevel-1
    initialmotilearea=sum(np.array(celltype)<0.7)
    initialmatrixarea=sum(np.array(celltype)>0.7)
    ##ploting initial condition
    #data = pd.DataFrame(data={'x':ypo, 'y':xpo, 'z':celltype})
    #data = data.pivot(index='x', columns='y', values='z')
    #fig, ax = plt.subplots(1, 1, figsize=[4, 4])
    #cMap = ListedColormap(['#00BF00','#BF00BF'])
    #sns.heatmap(data,cmap=cMap,linewidths=.5, ax=ax,cbar=False)
    #ax.axis(False)
    ##end of plot
    return xpo,ypo,growthtime,cellmass,celltype,initialmotilearea,initialmatrixarea

#checking different initial condition pattern

#xi,yi,growthtimei,cellmassi,celltypei,inimoi,inimai= fixedinitialcondition_6(1)
#data = pd.DataFrame(data={'x':yi, 'y':xi, 'z':celltypei})
#data = data.pivot(index='x', columns='y', values='z')
#fig, ax = plt.subplots(1, 1, figsize=[4, 4])

#cMap = ListedColormap(['#00BF00','#BF00BF'])
#sns.heatmap(data,cmap=cMap,linewidths=.5, ax=ax,cbar=False)
#ax.axis(False)
