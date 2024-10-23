import numpy as np
import matplotlib.pyplot as plt


def circle_dCondt(C,pos,const,D) -> float:
    nx,x,bx,ny,y,by = pos
    t,dt,dx,dy,R,dR,r = const
    
    dcdxdx = ((C[t][y][nx] - C[t][y][x]) + (-C[t][y][x]+ C[t][y][bx]))/(dx**2)

    dcdydy = (C[t][ny][x] - 2*C[t][y][x] + C[t][by][x])/(dy**2)
   
    dcdt = dcdxdx + dcdydy 
    
    C[t+1][y][x] = C[t][y][x] + dt*D*dcdt

def circle_dAdt(A_free,A_bound,C,C_bound,pos,const,D):
    nx,x,bx,ny,y,by = pos
    t,dt,dx,dy,k1,k2,R,dR,r = const

    dA_freedxdx = D*(A_free[t][y][nx] - 2*A_free[t][y][x] + A_free[t][y][bx])/(dx**2)
    dA_freedydy = D*(A_free[t][ny][x] - 2*A_free[t][y][x] + A_free[t][by][x])/(dy**2)

    dA_freedt = (dA_freedxdx + dA_freedydy) - k1*A_free[t][y][x]*C[t][y][x] + k2*A_bound[t][y][x]
    
    A_free[t+1][y][x] = A_free[t][y][x] + dt*dA_freedt
    

    dAbounddt = k1*A_free[t][y][x]*C[t][y][x] - k2*A_bound[t][y][x]
    A_bound[t+1][y][x] = A_bound[t][y][x] + dt*dAbounddt

    #C_bound[t+1][y][x] += 4*dAbounddt
    #C[t+1][y][x] += -4*dAbounddt

def circle_dAbounddt(A_free,A_bound,C,pos,k1,k2) -> float:
    nx,x,bx,ny,y,by = pos
    a = C[y][x]
    #a = 1
    dAbounddt = k1*A_free[y][x]*a - k2*A_bound[y][x]
    
    return dAbounddt

def stabil_condi(dt,dx,dy,D_list):
    # Von Neumann stability condition
    const1 = (1/dx)**2 + (1/dy)**2
    for D in D_list:
        const2 = 2*D*const1
        for _ in range(0,20):
            if dt > 1/const2:
                dt *= 0.9
        print(f"D={D} and dt={dt}")

    print(f"final dt={dt}")
    return dt

def open_close_membrane(Grid
                        ,Radius:int,dRadius:int
                        ,offsets,holesize:int
                        ,open_val:int ,wall_val:int
                        ,open_wall:bool
                        ):
    
    boxlen = np.shape(Grid)[0]
    x0, y0 = offsets
    xlowlim, xuplim  = int( x0 + Radius), int( x0 + Radius + 2*dRadius)
    ylowlim, yuplim  = int( y0 - holesize)  ,int( y0 + holesize)

    for x in range(xlowlim ,xuplim):
        for y in range(ylowlim ,yuplim):
            r = np.sqrt( (x-x0)**2 + (y-y0)**2 )
            #if Radius -1 < r <= Radius + dRadius +1:
            if open_wall == True:
                if Grid[y][x] == wall_val:
                    Grid[y][x] = open_val

            if open_wall == False:
                if Grid[y][x] == open_val:
                    Grid[y][x] = wall_val


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
def cicle_boundary(x:int,y:int,boxlen:int,ref_matrix, refval):
    nx, ny = x+1 , y+1
    bx, by = x-1 , y-1

    # circle in the middle.
    if ref_matrix[y][nx] == refval:
        nx = x
    if ref_matrix[y][bx] == refval:
        bx = x
    if ref_matrix[ny][x] == refval:
        ny = y
    if ref_matrix[by][x] == refval:
        by = y
    
    return nx,x,bx ,ny,y,by  


def init_conc(ref_grid, time:int
              ,c_in:float ,c_out:float
              ,inside_val:int ,outside_val:int
              ):
    ref_shape = np.shape(ref_grid)[0] # assuming NxN box shape.
    Grid = np.zeros(shape=(time,ref_shape,ref_shape),dtype=float)
    
    for x in range(ref_shape):
        for y in range(ref_shape):

            if ref_grid[y][x] == inside_val:
                Grid[0][y][x] = c_in
            
            if ref_grid[y][x] == outside_val:
                Grid[0][y][x] = c_out

            #if x == 0 or y == 0 or x == ref_shape-1 or y == ref_shape-1:
             #   Grid[0][y][x] = c_out
    return Grid

def init_ref_circle(boxlen:int
                ,Radius:int ,dRadius:int
                ,offsets:list
                ,inside_val:int,outside_val:int,wall_val:int
                ):
    
    x0,y0 = offsets
    ref_Grit = np.zeros(shape=(boxlen,boxlen),dtype=float)
    
    for y in range(boxlen):
        for x in range(boxlen):
            radii = np.sqrt(    (x-x0)**2 + (y-y0)**2   )

            if Radius <= radii <= Radius + dRadius:
                ref_Grit[y][x] = wall_val
            
            if radii < Radius:
                ref_Grit[y][x] = inside_val
            
            if Radius + dRadius < radii:
                ref_Grit[y][x] = outside_val
            
            if x == 0 or y == 0 or x ==boxlen-1 or y == boxlen-1:
                # this one comes last, as to override the value from the
                # if statement -> Radius + dRadius < radii:
                ref_Grit[y][x] = wall_val 

    print("shape of refgrid =",np.shape(ref_Grit))
    return ref_Grit



def Annexin_stablilization(
        k1:float,k2:float
        ,A_fo:float
        ,c_in:float
        ,realtime:float
        ,dt:float
        )->list:
    
    K = k1*c_in/k2
    C = A_fo/(1+K)
    b = k1*c_in + k2

    A_f = []
    A_b = []

    for t in np.arange(0,realtime,dt):
        #A_f.append( C*( 1 + np.exp(-b*t) ))
        A_f.append(
            C + (A_fo - C)*np.exp(-b*t)
        )
        A_b.append(
            ((k1*c_in*A_fo)/b)*( 1   -  np.exp(-b*t) )
            )

    return A_f, A_b

def points_inside_cell(refgrid,inside_val):

    ylen,xlen = np.shape(refgrid)
    point_count = 0
    for x in range(xlen):
        for y in range(ylen):
           if refgrid[y][x] == inside_val:
               point_count += 1
    
    return point_count
            


if __name__ == "__main__":
    
    pass
