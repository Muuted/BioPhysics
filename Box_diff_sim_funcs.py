import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

def check_volume(grid)-> float:
    Ysize,Xsize = np.shape(grid)
    
    tot_conc = 0
    for x in range(Xsize):
        for y in range(Ysize):
            tot_conc += grid[y][x]
    
    return tot_conc

def boundary_conditions(x,y,Xsize,Ysize,holesize,open):
    nx, ny = x + 1, y + 1
    bx,by = x - 1, y - 1

    if bx < 0 or bx == Ysize+1:
        bx = x

    if nx == Ysize+1 or nx == Xsize:
        nx = x

    if by < 0:
        by = y
    if ny == Ysize:
        ny = y
    
    if open == True:
        if x == Ysize+1 and  by <= (Ysize - holesize)/2 :
            by = y
        if x == Ysize+1 and  ny >= (Ysize + holesize)/2 :
            ny = y

        if (Ysize-holesize)/2 < y < (Ysize+holesize)/2:
            bx = x - 1
            nx = x + 1
            
            if bx < 0:
                bx = x
            if nx == Xsize:
                nx = x

    return nx,x,bx , ny,y,by

def dCondt(C,pos,dx,dy,h) -> float:

    nx,x,bx,ny,y,by = pos
    dcdxdx = (C[y][nx] - 2*C[y][x] + C[y][bx])/(dx**2)

    dcdydy = (C[ny][x] - 2*C[y][x] + C[by][x])/(dy**2)
   
    dcdt = dcdxdx + dcdydy 
    
    return dcdt

    
def test():
    #Testing diff equation setup.
    #print(Grid_diff[0])
    N = 10
    test_grid = np.zeros(shape=(2,N,N))
    t = 0
    T= 100
    test_grid[0][int(N/2)][int(N/2)] = 1000

    D,dt = 0.1,0.1 # there is something about numerical stability
    dx = 1         # where dt/(dx)^2 < 0.5 or something like that.
    dy = dx

    plt.matshow(test_grid[0])


    for t in range(T):
        if t%10==0:
            print(f"time ellapsed={t}, of total time={T}")


        t1, t2 = t%2, (t+1)%2
        
        for x in range(1,N-1):
            for y in range(1,N-1):
                change_val = dt*D*dCondt(C=test_grid[t1]
                                    ,x=x,y=y,dx=dx,dy=dy
                                    )
                test_grid[t2][y][x] = test_grid[t1][y][x] + change_val


    plt.matshow(test_grid[0])
    plt.show()





if __name__ == "__main__":
    pass
