from random import seed
from random import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
def whetheratedge(xpo,ypo,mtype):
    # finding position of cells on the edge of the biofilm
    edgex=[]
    edgey=[]
    edgetype=[]
    d2edge=[]
    loc=0
    while loc<len(xpo):
        e1=0
        e2=0
        e3=0
        e4=0
        notatedge=0
        n=0
        while n<len(ypo):
            if xpo[loc]==xpo[n]:
                if ypo[loc]==ypo[n]-1:
                    if e1==0:
                        notatedge=notatedge+1
                        e1=1
            if xpo[loc]==xpo[n]:
                if ypo[loc]==ypo[n]+1:
                    if e2==0:
                        notatedge=notatedge+1
                        e2=1
            if ypo[loc]==ypo[n]:
                if xpo[loc]==xpo[n]-1:
                    if e3==0:
                        notatedge=notatedge+1
                        e3=1
            if ypo[loc]==ypo[n]:
                if xpo[loc]==xpo[n]+1:
                    if e4==0:
                        notatedge=notatedge+1
                        e4=1
            n=n+1
        if notatedge<4:
            edgex.append(xpo[loc])
            edgey.append(ypo[loc])
            edgetype.append(mtype[loc])
            d2edge.append(np.array([xpo[loc],ypo[loc]]))
        loc=loc+1
        
    return edgex,edgey,edgetype,d2edge

def findmotile(xpo=[1],ypo=[1],celltype=[1]):
    #find positions of motile cell in the biofilm
    a=np.copy(celltype)
    b=np.where(np.array(a)<0.7)[0]
    motilex=[]
    motiley=[]
    motilepo=[]
    motiletype=[]
    for i in b:
        motilex.append(xpo[i])
        motiley.append(ypo[i])
        motilepo.append([xpo[i],ypo[i]])
        motiletype.append(celltype[i])
    #print(motilepo)
    #plt.imshow(xpo,ypo,celltype)
    return motilex,motiley,motiletype

