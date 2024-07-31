from random import seed
from random import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
def binaryimage(xpo,ypo,celltype):
    #return array for 2-D binary image plotting
    maxradius=max(max(xpo),max(ypo),-min(xpo),-min(ypo))
    xlength=max(xpo)-min(xpo)+1
    ylength=max(ypo)-min(ypo)+1
    if xlength>ylength:
        newxpo=np.array(xpo)-min(xpo)
        extra=np.floor((xlength-ylength)/2)
        newypo=np.array(ypo)-min(ypo)+extra
        length=xlength
    elif xlength<ylength:
        newypo=np.array(ypo)-min(ypo)
        extra=np.floor((ylength-xlength)/2)
        newxpo=np.array(xpo)-min(xpo)+extra
        length=ylength
    elif xlength==ylength:
        newxpo=np.array(xpo)-min(xpo)
        newypo=np.array(ypo)-min(ypo)
        length=xlength
    #print(newxpo,newypo,newlength)
    base=np.zeros((length,length))
    n=0
    while n<len(xpo):
        xloc=int(newxpo[n])
        yloc=int(newypo[n])
        base[yloc][xloc]=1
        n=n+1
    #print(base)
    return base
