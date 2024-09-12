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



def dCondt(C,x,y,dx,dy,h) -> float:

    dcdxdx = (C[y][x + dx] - 2*C[y][x] + C[y][x-dx])/(dx**2)

    dcdydy = (C[y + dx][x] - 2*C[y][x] + C[y-dx][x])/(dy**2)
   
    dcdt = dcdxdx + dcdydy 
    
    return dcdt


def RungeKutta(C,x,y,dx,dy,h):

    #first part Runge Kutta
    dx,dy = 2*dx, 2*dy
    
    dcdxdx = (C[y][x+dx] - 2*C[y][x] + C[y][x-dy] )/(dx**2)

    dcdydy = (C[y+dy][x] - 2*C[y][x] + C[y-dx][x])/(dy**2)

   
    Cxy = C[y][x] + ( dx*dcdxdx/2 + dy*dcdydy/2)/2
    
    #k1x,ky1 = dcdxdx, dcdydy 

    # 2nd part Runge kutta
    dx, dy = dx/2 , dy/2
    dcdxdx = (C[y][x+dx] - 2*Cxy + C[y][x-dy] )/(dx**2)

    dcdydy = (C[y+dy][x] - 2*Cxy + C[y-dx][x])/(dy**2)

    


    
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
