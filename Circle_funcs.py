import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

def replenish(Grid,boxsize:int,c_out):

    for x in range(boxsize):
        Grid[0][x] = c_out
        Grid[1][x] = c_out
        Grid[boxsize-2][x] = c_out
        Grid[boxsize-1][x] = c_out

    for y in range(boxsize):
        Grid[y][0] = c_out
        Grid[y][1] = c_out
        Grid[y][boxsize-2] = c_out
        Grid[y][boxsize-1] = c_out
    
    return Grid
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


def init_circle(boxlen:int,time:int
                , Radius:int ,dRadius:int
                ,c_in:float, c_out:float
                ,refval, center:int
                ):

    Grit = np.zeros(shape=(time,boxlen,boxlen))
    ref_Grit = np.zeros(shape=(boxlen,boxlen))
    
    for y in range(boxlen):
        for x in range(boxlen):
            radii = np.sqrt(    (x-center)**2 + (y-center)**2   )

            if radii < Radius :
                Grit[0][y][x] = c_in
            if radii > Radius + dRadius:
                Grit[0][y][x] = c_out

            if Radius <= radii <= Radius + dRadius:
                ref_Grit[y][x] = refval

    return Grit, ref_Grit

def circle_dCondt(C,pos,dx,dy) -> float:
    nx,x,bx,ny,y,by = pos

    dcdxdx = (C[y][nx] - 2*C[y][x] + C[y][bx])/(dx**2)

    dcdydy = (C[ny][x] - 2*C[y][x] + C[by][x])/(dy**2)
   
    dcdt = dcdxdx + dcdydy 
    
    return dcdt



def test_init_circle():
    sidelen = 100
    T_tot = 100
    wall_val = 100
    conc_list = []
    R = sidelen*0.7/2
    dR = 2
    dx,dy = 1, 1
    dt,D = 0.1, 0.5
    center = int(sidelen/2)

    Grid,ref_Grid = init_circle(boxlen=sidelen,time=T_tot
                                ,R=(sidelen*0.6)/2,dR=2
                                ,c_in=0,c_out=10
                                ,refval=100,center=center
                                )

    Grid[0][center][center] = 1000

    for t in np.arange(0,T_tot-1):
        t1, t2 = t%2, (t+1)%2 # to change the matricies back and forth
        for x in range(0,sidelen):
            for y in range(0,sidelen):
                radii = np.sqrt(
                    (x-center)**2 + (y-center)**2
                )

                if ref_Grid[y][x] != wall_val:
                    pos = cicle_boundary(x=x,y=y,boxlen=sidelen
                                            ,ref_matrix=ref_Grid
                                            ,refval=wall_val)

                    dcdt = circle_dCondt(C=Grid[t],pos=pos,dx=dx,dy=dy)

                    Grid[t+1][y][x] = Grid[t][y][x] + dt*D*dcdt
                
        if t%(int(T_tot*0.1)) == 0:
            conc_list.append(
                np.sum(Grid[t])
                    )


    plt.matshow(ref_Grid)
    plt.matshow(Grid[T_tot-1])

    plt.figure()
    plt.plot(conc_list)
    plt.show()


if __name__ == "__main__":
    #test_init_circle()
    pass
