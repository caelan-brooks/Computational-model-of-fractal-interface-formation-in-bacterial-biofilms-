import numpy as np
import matplotlib.pyplot as plt
import math

def findmotileradius(moarea=1):
    # moarea is the amount of motile cells
    radius=(moarea/math.pi)**0.5
    return radius
def findbiofilmradius(moarea=1,maarea=1):
    # moarea is the amount of motile cells, maarea is the amount of matrix producing cell
    radius=((moarea+maarea)/math.pi)**0.5
    return radius

def findmotileradiusgrowth(moarea=1,moini=1):
    # moarea is the amount of motile cells, moini is the initial amount of motile cells
    radius=(moarea/math.pi)**0.5-((moini/math.pi)**0.5)
    return radius
def findbiofilmradiusgrowth(moarea=1,moini=1,maarea=1,maini=1):
    # moarea is the amount of motile cells, moini is the initial amount of motile cells
    # moarea is the amount of matrix producing cells, moini is the initial amount of matrix producing cells
    radius=((moarea+maarea)/math.pi)**0.5-(((moini+maini)/math.pi)**0.5)
    #print(moarea,maarea,moini,maini,radius)
    return radius
