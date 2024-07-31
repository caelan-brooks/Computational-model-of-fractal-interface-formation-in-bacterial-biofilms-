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
def noswitching(xpo=[1],ypo=[1],growthtime=[1],cellmass=[1],celltype=[1],massratio=1,shovecp=1):
    #initial condition imported from the input, together with mass ratio and shoving capacity
    #this function is supposed to iterate over time until the biofilm inner pattern freezes
    #output includes information about the biofilm at single cell level with estimated size and time when pattern freezes
    u=0 #u is used for graphing biofilm simulation at different time

    #parameter
    masstopush=shovecp
    norfacx=1  #growthtime for matrix producing cells
    norfacy=1  #growthtime for motile cells
    norfacmx=1*massratio #mass for matrix producing cells
    norfacmy=1 #mass for motile cell

    #looping over time
    recordpointt=1 #time point used to plot biofilm at specific time
  
    patternfreeze=0 #used to break the loop when pattern freezes
    tadvanced=0 #indicator of time
    celltypebefore=sum(np.where(np.array(celltype)<1)[0]) #sum of all motile cell, used to check whether the pattern freezes
    timepoint=1 #indicator of check point -- when time reaches certain timepoint, see if celltypebefore is the same
    areamotile=[] #record of amount of motile producing cells
    areamatrix=[] #record of amount of matrix cells
    areatime=[] #record of time when the amount of different cells is recorded
    areamotile.append(sum(np.array(celltype)<0.7))
    areamatrix.append(sum(np.array(celltype)>0.7))
    areatime.append(tadvanced)



    while patternfreeze ==0:
        minustime=min(growthtime)
        k=0
        #print(1)
        while k<len(growthtime):
            if growthtime[k]==minustime:
                
                #atedge=0
                cx=0
                cy=0
                atedge=[0,0,0,0] # np.array(cdown,cup,cright,cleft)
                while cx<len(xpo):
                    if xpo[cx]==xpo[k]:
                        if atedge[0]==0:
                            if ypo[cx]<ypo[k]:
                                #atedge=atedge+1000
                                atedge[0]=1
                        if atedge[1]==0:
                            if ypo[cx]>ypo[k]:
                                #atedge=atedge+100
                                atedge[1]=1
                    if ypo[cx]==ypo[k]:
                        if atedge[2]==0:
                            if xpo[cx]>xpo[k]:
                                atedge[2]=1
                        if atedge[3]==0:
                            if xpo[cx]<xpo[k]:
                                atedge[3]=1
                            
                    cx=cx+1
                    
                direction=0
                if sum(atedge)<4:
                    possibledirection=np.where(np.array(atedge)<0.5)[0]
                elif sum(atedge)==4:
                    possibledirection=np.where(np.array(atedge)>0.5)[0]
                    
                r=random()
                portion=len(possibledirection)
                chosen=0
                if r<1/portion:
                    chosen=0
                elif r<2/portion:
                    chosen=1
                elif r<3/portion:
                    chosen=2
                elif r<4/portion:
                    chosen=3
                direction=possibledirection[chosen]
                    
                if direction==0:
                    if r<1/portion:
                        j=0
                        csh=masstopush
                        massunit=[]
                        uplimi=ypo[k]
                        while j<len(ypo):
                            if xpo[k]==xpo[j]:
                                if ypo[j]<ypo[k]:
                                    massunit.append(j)
                            j=j+1
                        n=1

                        while n>0:
                            #print(2)
                            po=1
                            for m in massunit:
                                if ypo[m]==ypo[k]-n:
                                    csh=csh-cellmass[m]
                                    po=0
                                    if csh<0:
                                        po=1
                                        break
                                    break

                                else:
                                    po=1
                            n=n+1
                            if po==1:
                                break
                        n=n-1            
                        if csh>0:
                            #for m in massunit:
                             #   if n>0:
                              #      if ypo[m]==ypo[k]-n:
                               #         ypo[m]=ypo[m]-1
                                #        n=n-1
                            while n>0:
                                for m in massunit:
                                    if ypo[m]==ypo[k]-n:
                                        ypo[m]=ypo[m]-1
                                        break
                                n=n-1
                            #print("up")
                            ypo.append(ypo[k]-1)
                            xpo.append(xpo[k])
                            celltype.append(celltype[k])
                            if celltype[k]==1:
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
                                        growthtime.append(o+minustime)
                                        break
                                    else:
                                        f=f+1
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacx, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1
                            elif celltype[k]==0.5:
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
                                        growthtime.append(o+minustime)
                                        break
                                    else:
                                        f=f+1
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacy, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1
                        else:
                            if celltype[k]==0.5:
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacy, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1
                            elif celltype[k]==1:
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacx, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1

                if direction==1:
                    if r<2/portion:
                        j=0
                        csh=masstopush
                        massunit=[]
                        uplimi=ypo[k]
                        while j<len(ypo):
                            if xpo[k]==xpo[j]:
                                if ypo[j]>ypo[k]:
                                    massunit.append(j)
                            j=j+1
                        n=1

                        while n>0:
                            po=1
                            for m in massunit:
                                if ypo[m]==ypo[k]+n:
                                    csh=csh-cellmass[m]
                                    po=0
                                    if csh<0:
                                        po=1
                                    break

                                else:
                                    po=1
                            n=n+1
                            if po==1:
                                break
                        n=n-1            
                        if csh>0:
                            #for m in massunit:
                             #   if n>0:
                              #      if ypo[m]==ypo[k]+n:
                               #         ypo[m]=ypo[m]+1
                                #        n=n-1
                            while n>0:
                                for m in massunit:
                                    if ypo[m]==ypo[k]+n:
                                        ypo[m]=ypo[m]+1
                                        break
                                n=n-1
                            #print("down")
                            ypo.append(ypo[k]+1)
                            xpo.append(xpo[k])
                            celltype.append(celltype[k])
                            if celltype[k]==1:
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
                                        growthtime.append(o+minustime)
                                        break
                                    else:
                                        f=f+1
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacx, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1
                            if celltype[k]==0.5:
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
                                        growthtime.append(o+minustime)
                                        break
                                    else:
                                        f=f+1
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacy, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1
                        else:
                            if celltype[k]==0.5:
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacy, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1
                            elif celltype[k]==1:
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacx, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1

                if direction==2:
                    if r<3/portion:
                        j=0
                        csh=masstopush
                        massunit=[]
                        uplimi=xpo[k]
                        while j<len(xpo):
                            if ypo[k]==ypo[j]:
                                if xpo[j]>xpo[k]:
                                    massunit.append(j)
                            j=j+1
                        n=1

                        while n>0:
                            po=1
                            for m in massunit:
                                if xpo[m]==xpo[k]+n:
                                    csh=csh-cellmass[m]
                                    po=0
                                    if csh<0:
                                        po=1
                                    break

                                else:
                                    po=1
                            n=n+1
                            if po==1:
                                break

                        n=n-1            
                        if csh>0:
                            while n>0:
                                for m in massunit:
                                    if xpo[m]==xpo[k]+n:
                                        xpo[m]=xpo[m]+1
                                        break

                                n=n-1
                            #print("right")
                            xpo.append(xpo[k]+1)
                            ypo.append(ypo[k])
                            celltype.append(celltype[k])
                            if celltype[k]==1:
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
                                        growthtime.append(o+minustime)
                                        break
                                    else:
                                        f=f+1
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacx, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1
                            if celltype[k]==0.5:
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
                                        growthtime.append(o+minustime)
                                        break
                                    else:
                                        f=f+1
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacy, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1
                        else:
                            if celltype[k]==0.5:
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacy, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1
                            elif celltype[k]==1:
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacx, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1


                if direction==3:
                    if r<4/portion:
                        j=0
                        csh=masstopush
                        massunit=[]
                        uplimi=xpo[k]
                        while j<len(xpo):
                            if ypo[j]==ypo[k]:
                                if xpo[j]<xpo[k]:
                                    massunit.append(j)
                            j=j+1
                        n=1

                        #for m in massunit:
                         #   if xpo[m]==xpo[k]-n:
                          #      csh=csh-cellmass[m]
                           #     n=n+1
                            #else:
                             #   break
                        while n>0:
                            po=1
                            for m in massunit:
                                if xpo[m]==xpo[k]-n:
                                    csh=csh-cellmass[m]
                                    po=0
                                    if csh<0:
                                        po=1
                                    break

                                else:
                                    po=1
                            n=n+1
                            if po==1:
                                break

                        n=n-1
                        #print(n)
                        if csh>0:
                            while n>0:
                                po=1
                                for m in massunit:
                                    if xpo[m]==xpo[k]-n:
                                        xpo[m]=xpo[m]-1
                                        po=0
                                        break
                                    else:
                                        po=1
                                n=n-1
                            #print("left")
                            #print(xpo[k])
                            xpo.append(xpo[k]-1)
                            ypo.append(ypo[k])
                            celltype.append(celltype[k])
                            if celltype[k]==1:
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
                                        growthtime.append(o+minustime)
                                        break
                                    else:
                                        f=f+1
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacx, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1
                            if celltype[k]==0.5:
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
                                        growthtime.append(o+minustime)
                                        break
                                    else:
                                        f=f+1
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacy, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1
                        else:
                            if celltype[k]==0.5:
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacy, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1
                            elif celltype[k]==1:
                                f=1
                                while f>0:
                                    o=np.random.normal(norfacx, 0.1)
                                    if o>0:
                                        growthtime[k]=o+minustime
                                        break
                                    else:
                                        f=f+1









            k=k+1


        tadvanced=tadvanced+minustime
        growthtime=np.array(growthtime)-minustime
        growthtime=growthtime.tolist()

        aarea=len(celltype)
        aradius=(aarea/(math.pi))**0.5
        alocmo=0
        motilecellpo=np.where(np.array(celltype)<0.8)[0]
        while alocmo<len(motilecellpo):
            aloc=motilecellpo[alocmo]
            if celltype[aloc]==0.5:
                xloc=xpo[aloc]
                yloc=ypo[aloc]
                poploc=0
                popi5=0
                popi6=0
                popi7=0
                popi8=0
                popi1=0
                popi2=0
                popi3=0
                popi4=0
                for i in xpo:
                    if ypo[poploc]==yloc:
                        if i<xloc:
                            popi1=1
                            break
                    poploc=poploc+1
                poploc=0            
                for i in xpo:
                    if ypo[poploc]==yloc:
                        if i>xloc:
                            popi2=1
                            break
                    poploc=poploc+1
                poploc=0            
                for i in ypo:
                    if xpo[poploc]==xloc:
                        if i>yloc:
                            popi3=1
                            break
                    poploc=poploc+1

                poploc=0
                for i in ypo:
                    if xpo[poploc]==xloc:
                        if i<yloc:
                            popi4=1
                            break
                    poploc=poploc+1

                if popi1*popi2*popi3*popi4 ==0:
                    celltype.pop(aloc)
                    xpo.pop(aloc)
                    ypo.pop(aloc)
                    growthtime.pop(aloc)
                    cellmass.pop(aloc)
                    motilecellpo=np.where(np.array(celltype)<0.8)[0]
                    alocmo=-1

            elif celltype[aloc]==1:
                xloc=xpo[aloc]
                yloc=ypo[aloc]
                poploc=0
                popi5=0
                popi6=0
                popi7=0
                popi8=0
                popi1=0
                popi2=0
                popi3=0
                popi4=0

                for i in xpo:
                    if ypo[poploc]==yloc:
                        if i<xloc:
                            popi5=1
                            break
                    poploc=poploc+1
                poploc=0            
                for i in xpo:
                    if ypo[poploc]==yloc:
                        if i>xloc:
                            popi6=1
                            break
                    poploc=poploc+1
                poploc=0            
                for i in ypo:
                    if xpo[poploc]==xloc:
                        if i>yloc:
                            popi7=1
                            break
                    poploc=poploc+1

                poploc=0
                for i in ypo:
                    if xpo[poploc]==xloc:
                        if i<yloc:
                            popi8=1
                            break
                    poploc=poploc+1

            alocmo=alocmo+1

        areamotile.append(sum(np.array(celltype)<0.7))
        areamatrix.append(sum(np.array(celltype)>0.7))
        areatime.append(tadvanced)
        if tadvanced>=3*timepoint:
            #print("checking pattern")
            newcelltype=sum(np.where(np.array(celltype)<1)[0])
            if celltypebefore==newcelltype:
                patternfreeze=1
                frozenmotile=max(areamotile)
                frozenpointarray=np.where(np.array(areamotile)>=frozenmotile)[0]
                frozenpoint=frozenpointarray[0]
                frozenmotilearea=areamotile[frozenpoint]
                frozenmatrixarea=areamatrix[frozenpoint]
                frozenareatime=areatime[frozenpoint]
                
                #plot of amount of cells over time
                #plt.figure(u)
                #plt.plot(areatime,((np.array(areamatrix)+np.array(areamotile))/math.pi)**0.5,"r--")
                #plt.plot(areatime,(np.array(areamotile)/math.pi)**0.5,"b--")
                #plt.xlabel("time")
                #plt.ylabel("radius")
                #u=u+1
                
            else:
                patternfreeze=0
                celltypebefore=newcelltype
                timepoint=timepoint+1

        #plotting biofilm at time=recordpointt
        if tadvanced>recordpointt:
            if tadvanced-minustime<recordpointt:
                #print(tadvanced,3*timepoint)
                #plt.figure(u)
                #x=xpo
                #y=ypo
                #z=celltype

                #data = pd.DataFrame(data={'x':y, 'y':x, 'z':z})
                #data = data.pivot(index='x', columns='y', values='z')
                #sns.heatmap(data,linewidths=0.5, linecolor='white', cbar=True, square=True)
                #plt.show()
                #u=u+1
                recordpointt=recordpointt+1
    ##ploting biofilm at sometime after inner pattern freezes
    #u=u+1
    #plt.figure(u)
    #x=xpo
    #y=ypo
    #z=celltype

    #data = pd.DataFrame(data={'x':y, 'y':x, 'z':z})
    #data = data.pivot(index='x', columns='y', values='z')
    #sns.heatmap(data,linewidths=0.5, linecolor='white', cbar=True, square=True)
    #plt.show()
    ##end of plot
    ##a different way to plot
    #data = pd.DataFrame(data={'x':ypo, 'y':xpo, 'z':celltype})
    #data = data.pivot(index='x', columns='y', values='z')
    #fig, ax = plt.subplots(1, 1, figsize=[4, 4])
    #cMap = ListedColormap(['#00BF00','#BF00BF'])
    #sns.heatmap(data,cmap=cMap,linewidths=.5, ax=ax,cbar=False)
    #ax.axis(False)
    ##end of plot
    return xpo,ypo,growthtime,cellmass,celltype,frozenmotilearea,frozenmatrixarea,frozenareatime
