import numpy as np
import matplotlib.pyplot as plt


def Ring_sum(boxsize:int, offsets
             #,xstart:int,ystart:int
             ,Radius:int,dR:int
             ):
    xstart = offsets[0] + Radius
    xend = xstart - 2*Radius
    ystart = offsets[1] - Radius
    yend = offsets[1] + Radius
    for x in range(xstart,xend,-1):
        for y in range(ystart,yend):
            pass

    

def stabil_condi(dt,dx,dy,D):
    # Von Neumann stability condition
    const1 = (1/dx)**2 + (1/dy)**2
    const2 = 1/(2*D*const1)

    for _ in range(0,20):
        if dt > const2:
            dt *= 0.1
        
    return dt


def open_close_membrane(Grid
                        ,Radius:int,dRadius:int
                        ,offsets,holesize:int
                        ,open_val:int ,wall_val:int
                        ,open_wall:bool
                        ):
    
    boxlen = np.shape(Grid)[0]
    x0, y0 = offsets
    xlowlim, xuplim  = int( x0 + Radius - holesize), int( x0 + Radius + holesize)
    ylowlim, yuplim  = int( y0 - holesize)  ,int( y0 + holesize)

    for x in range(xlowlim ,xuplim):
        for y in range(ylowlim ,yuplim):
            if open_wall == True:
                if Grid[y][x] == wall_val:
                    Grid[y][x] = open_val

            if open_wall == False:
                if Grid[y][x] == open_val:
                    Grid[y][x] = wall_val
    return Grid

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


def init_conc(ref_grid, time:int
              ,c_in:float ,c_out:float
              ,inside_val:int ,outside_val:int
              ):
    ref_shape = np.shape(ref_grid)[0] # assuming NxN box shape.
    Grid = np.zeros(shape=(time,ref_shape,ref_shape))
    
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
    ref_Grit = np.zeros(shape=(boxlen,boxlen))
    
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

    return ref_Grit

def circle_dCondt(C,pos,dx,dy) -> float:
    nx,x,bx,ny,y,by = pos

    dcdxdx = (C[y][nx] - 2*C[y][x] + C[y][bx])/(dx**2)

    dcdydy = (C[ny][x] - 2*C[y][x] + C[by][x])/(dy**2)
   
    dcdt = dcdxdx + dcdydy 
    
    return dcdt




if __name__ == "__main__":
    
    pass
