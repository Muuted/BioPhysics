import numpy as np
import matplotlib.pyplot as plt
import scipy as sp




def dCdt(C,D,x,y,dx,dy,dt):

    dcdxdx = (C[y][x + dx] - 2*C[y][x] + C[y][x-dx])/(dx^2)

    dcdydy = (C[y + dy][x] - 2*C[y][x] + C[y-dy][x])/(dy^2)
   
    dcdt = dcdxdx + dcdydy 
    
    return dcdt

    




"""
if __name__ == "__main__":
    pass
"""