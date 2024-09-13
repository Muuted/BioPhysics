import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


def cicle_boundary(x:int,y:int,boxlen:int,ref_matrix,refval):
    nx, ny = x+1 , y+1
    bx, by = x-1 , y-1

    #outer box boundary first
    if bx < 0:
        bx = x
    if nx == boxlen:
        nx = x
    if by < 0:
        by = y
    if ny == boxlen:
        ny = y

    # circle in the middle.
    if ref_matrix[y][nx] == refval:
        nx = x
    if ref_matrix[y][bx] == refval:
        bx = x
    if ref_matrix[ny][x] == refval:
        ny = y
    if ref_matrix[by][x] == refval:
        by = y
    
    return nx,x,bx , ny,y,by


def init_circle(boxlen:int, R:int ,dR:int, c_in:float, c_out:float):

    Grit = np.zeros(shape=(boxlen,boxlen))
    center = int(boxlen/2)

    for y in range(boxlen):
        for x in range(boxlen):

            radii = np.sqrt(
                (x-center)**2 + (y-center)**2
            )
            if radii < R :
                Grit[y][x] = c_in
            if radii > R + dR:
                Grit[y][x] = c_out

    return Grit

def circle_dCondt(C,pos,dx,dy,h) -> float:

    nx,x,bx,ny,y,by = pos
    dcdxdx = (C[y][nx] - 2*C[y][x] + C[y][bx])/(dx**2)

    dcdydy = (C[ny][x] - 2*C[y][x] + C[by][x])/(dy**2)
   
    dcdt = dcdxdx + dcdydy 
    
    return dcdt



def test_init_circle():
    sidelen = 100
    Grit = init_circle(boxlen=sidelen,R=(sidelen*0.6)/2,dR=2,c_in=2,c_out=10)

    plt.matshow(Grit)

    plt.show()


if __name__ == "__main__":
    test_init_circle()
    pass
