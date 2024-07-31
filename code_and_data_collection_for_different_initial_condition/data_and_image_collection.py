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
from initial_conditions import *
from biofilm_expansion_function import *
from radius_calculation_function import *
from binary_image_conversion_function import *
from fractal_dimension_calculation_function import *
from find_edgeormotile_cell_function import *

massratiocollect=[1,2,3,4,5,10]#,1.5,2,3,5,10,100]#[5,2,1]
shovingcapacity=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]#[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]#[5,10,15]#
fractaldimensioncollection=[]
totalsample=5
fractalwholesample=[]
motileradiustot=[]
matrixradiustot=[]
motileradiusgrowthtot=[]
biofilmradiusgrowthtot=[]
from matplotlib.colors import ListedColormap
initialtotalnumberofcell=0
initialratioofcell=0
for massr in massratiocollect:
    fractaldmass=[]
    fractalwholesamplemass=[]
    motileradiusmass=[]
    matrixradiusmass=[]
    motileradiusgrowthmass=[]
    biofilmradiusgrowthmass=[]
    for shovecp in shovingcapacity:
        samplefractal=[]
        totaltrial=totalsample
        motileradiussingle=[]
        matrixradiussingle=[]
        motileradiusgrowthsingle=[]
        biofilmradiusgrowthsingle=[]
        while totaltrial>0:
            totaltrial=totaltrial-1
            massratio=massr
            xpo,ypo,growthtime,cellmass,celltype,inimo,inima= fixedinitialcondition_1(massratio)
            xpo,ypo,growthtime,cellmass,celltype,moa,maa,ftime=noswitching(xpo,ypo,growthtime,cellmass,celltype,massratio,shovecp)
            #data = pd.DataFrame(data={'x':ypo, 'y':xpo, 'z':celltype})
            #data = data.pivot(index='x', columns='y', values='z')
            #fig, ax = plt.subplots(1, 1, figsize=[3, 3])
            #cMap = ListedColormap(['#00BF00','#BF00BF'])
            #im=plt.imshow(data,cmap=cMap)
            #ax.axis(False)
            #plt.title("whole bioflim "+str("massratio=")+str(massr)+str("shovingcapacity=")+str(shovecp),fontsize = 6)
            #plt.savefig("whole bioflim"+str("massratio=")+str(massr)+str("shovingcapacity=")+str(shovecp)+".pdf")
            motilex,motiley,motiletype=findmotile(xpo,ypo,celltype)
            edgex,edgey,edgetype,d2edge=whetheratedge(motilex,motiley,motiletype)
            #fig, ax = plt.subplots(1, 1, figsize=[3, 3])
            #data = pd.DataFrame(data={'x':edgey, 'y':edgex, 'z':edgetype})
            #data = data.pivot(index='x', columns='y', values='z')
            #ax.axis(False)
            #coloor=ListedColormap(['#00BF00'])
            #im = plt.imshow(data, cmap=coloor)
            #plt.title("motile edge "+str("massratio=")+str(massr)+str("shovingcapacity=")+str(shovecp),fontsize = 6)
            #plt.title("motile edge black"+str("massratio=")+str(massr)+str("shovingcapacity=")+str(shovecp),fontsize = 6)
            #plt.savefig("motile edge black"+str("massratio=")+str(massr)+str("shovingcapacity=")+str(shovecp)+".pdf")
            Z=binaryimage(edgex,edgey,edgetype)
            fractaldsingle=fractal_dimension(Z)
            
            motileradiussingle.append(findmotileradius(moa))
            matrixradiussingle.append(findbiofilmradius(moa,maa))
            motileradiusgrowthsingle.append(findmotileradiusgrowth(moa,inimo))
            biofilmradiusgrowthsingle.append(findbiofilmradiusgrowth(moa,inimo,maa,inima))
            samplefractal.append(fractaldsingle)
        fractaldsingletot=sum(samplefractal)/len(samplefractal)
        fractaldmass.append(fractaldsingletot)
        fractalwholesamplemass.append(samplefractal)
        print("massratio=",massr,"shovecap=",shovecp,"avg fd=",fractaldsingletot)
        motileradiusmass.append(motileradiussingle)
        matrixradiusmass.append(matrixradiussingle)
        motileradiusgrowthmass.append(motileradiusgrowthsingle)
        biofilmradiusgrowthmass.append(biofilmradiusgrowthsingle)
    fractaldimensioncollection.append(fractaldmass)
    fractalwholesample.append(fractalwholesamplemass)
    motileradiustot.append(motileradiusmass)
    matrixradiustot.append(matrixradiusmass)
    motileradiusgrowthtot.append(motileradiusgrowthmass)
    biofilmradiusgrowthtot.append(biofilmradiusgrowthmass)
initialtotalnumberofcell=inimo+inima
initialratioofcell=inimo/inima
print(initialtotalnumberofcell,initialratioofcell)

#plt.figure(1)
#bp=0
#colorbar=['b','g','r','c','m','y','k']
#while bp<len(massratiocollect):
#    plt.plot(shovingcapacity,fractaldimensioncollection[bp],color=colorbar[bp], marker='o', linestyle='dashed',label=massratiocollect[bp])
#    bp=bp+1
#plt.legend()
#plt.xlabel("shoving capacity")
#plt.ylabel("fractal dimension")
#plt.title("fractal dimension vs. shoving capacity with various mass ratio")

#data collection and plotting initial condition
#print("initial total number of cell+spore=",initialtotalnumberofcell)
#print("initial ratio of cell/spore",initialratioofcell)
#print("fractal dimension:",fractalwholesample)
#print("motile radius:",motileradiustot)
#print("biofilm radius:",matrixradiustot)
#print("average motile radius growth:",motileradiusgrowthtot)
#print("average biofilm radius growth:",biofilmradiusgrowthtot)
#xi,yi,growthtimei,cellmassi,celltypei,inimoi,inimai= fixedinitialcondition_6(1)
#data = pd.DataFrame(data={'x':yi, 'y':xi, 'z':celltypei})
#data = data.pivot(index='x', columns='y', values='z')
#fig, ax = plt.subplots(1, 1, figsize=[4, 4])

#cMap = ListedColormap(['#00BF00','#BF00BF'])
#sns.heatmap(data,cmap=cMap,linewidths=.5, ax=ax,cbar=False)
