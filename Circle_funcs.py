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

"""def circle_dAbounddt(A_free,A_bound,C,pos,k1,k2) -> float:
    nx,x,bx,ny,y,by = pos
    a = C[y][x]
    #a = 1
    dAbounddt = k1*A_free[y][x]*a - k2*A_bound[y][x]
    
    return dAbounddt"""

def stabil_condi(dt,dx,dy,D_list):
    # Von Neumann stability condition    
    N = 3000
    for D in D_list:
        for _ in range(0,N):
            if dt > (dx*dx)/(4*D):
                dt *= 9.5e-1
            if dt < (dx*dx)/(4*D):
                break
            if _ == N-1:
                print("Not enough iterations for dt")
    
    dt *= 9.5e-1
    #print(f"final dt={dt}")
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

def open_close_membrane2(
        Grid
        ,offsets
        ,Xholesize:int
        ,Yholesize:int
        ,open_val:int ,wall_val:int
        ,open_wall_bool:bool
                        ):
    x0,y0 = offsets
    xlow,xhigh= int(x0 - Xholesize) , int(x0 + Xholesize)
    ylow,yhigh = int(y0 - Yholesize) , int(y0 + Yholesize)
    
    for y in range(ylow,yhigh):
        for x in range(xlow,xhigh):
            
            if open_wall_bool == True:
                if Grid[y][x] == wall_val:
                    Grid[y][x] = open_val

            if open_wall_bool == False:
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
    ref_shapey,ref_shapex = np.shape(ref_grid) # assuming NxN box shape.
    Grid = np.zeros(shape=(time,ref_shapey,ref_shapex),dtype=float)
    
    for x in range(ref_shapex):
        for y in range(ref_shapey):

            if ref_grid[y][x] == inside_val:
                Grid[0][y][x] = c_in
            
            if ref_grid[y][x] == outside_val:
                Grid[0][y][x] = c_out

            if x == 0 or x == ref_shapex - 1 :
                Grid[0][y][x] = c_out
            if y == 0 or y == ref_shapey - 1 :
                Grid[0][y][x] = c_out
            
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
            
            if x == 0 or x == boxlen - 1:
                ref_Grit[y][x] = wall_val
            if y == 0 or y == boxlen - 1:
                ref_Grit[y][x] = wall_val
            

    print("shape of refgrid =",np.shape(ref_Grit))
    return ref_Grit



def Annexin_stablilization(
        k1:float,k2:float
        ,A_tot:float
        ,c_in:float
        ,realtime:float
        ,dt:float
        )->list:
    
    #K = k1*c_in/k2
    #C = A_tot/(1+K)
    b = k1*c_in + k2

    A_f = []
    A_b = []

    for t in np.arange(0,realtime,dt):
        A_f.append(
            (A_tot/b)*(k2 + k1*c_in*np.exp(-b*t ))
        )
        
        A_b.append(
            ((k1*c_in*A_tot)/b)*( 1 - np.exp(-b*t) )
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
            

def main_circle_sim(
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,open_hole=True
        ,ref_bakteria = ""
        ,data_path=""
                    ) -> list:
    
    # Creation of our Grids
    if ref_bakteria == "":
        ref_structure = init_ref_circle(
            boxlen=len_size
            ,Radius=R,dRadius=dR ,offsets=[x0,y0]
            ,inside_val=inside_val
            ,outside_val=outside_val
            ,wall_val=wall_val
                                        )
    else:
        py_ref_struct,outside_val,inside_val,wall_vall,xoffset,yoffset = make_ref_structure(
            path=data_path
            ,ref_name=ref_bakteria
        )


    Free_Ca = init_conc(
        ref_grid=ref_structure
        ,time=T_tot
        ,c_in=c_in
        ,c_out=c_out
        ,inside_val=inside_val
        ,outside_val=outside_val
                        )
    
    Free_annexin = init_conc(
        ref_grid=ref_structure
        ,time=T_tot
        ,c_in= c_in_annexin
        ,c_out=0
        ,inside_val=inside_val
        ,outside_val=outside_val
                            )

    Bound_annexin = init_conc(
        ref_grid=ref_structure
        ,time=T_tot
        ,c_in= bound_annexin_start
        ,c_out=0
        ,inside_val=inside_val
        ,outside_val=outside_val
                            )
    
    Bound_Ca = np.zeros(shape=(T_tot,len_size,len_size))

    open_close_membrane(
        Grid=ref_structure
        ,Radius=R, dRadius=dR
        ,offsets=[x0,y0],holesize=holesize
        ,open_val=open_val
        ,wall_val=wall_val
        ,open_wall=open_hole
                        )
    i = 0 # for showing that the theory steady state concentration  
          # matches the simulated one.
    for t in np.arange(0,T_tot-1): 
        if t%(int(T_tot/10)) == 0:
            #print(f"time={t} of {T_tot}")
            print(f" Simulation Progress :  {int((t/T_tot)*100)} %")   
        t1, t2 = t%2, (t+1)%2
        
        for x in range(0,len_size):
            for y in range(0,len_size):
                #Bound_Ca[t+1][y][x] += Bound_Ca[t][y][x]
                #if ref_structure[y][x] != wall_val:
                if ref_structure[y][x] == wall_val or ref_structure[y][x] == outside_val:
                    Free_Ca[t+1][y][x] = Free_Ca[t][y][x]
                    Free_annexin[t+1][y][x] = Free_annexin[t][y][x] 
                    Bound_annexin[t+1][y][x] = Bound_annexin[t][y][x]

                if ref_structure[y][x] == inside_val or ref_structure[y][x] == open_val:
                    radii = np.sqrt( (x-x0)**2 + (y-y0)**2 )
                    pos = cicle_boundary(x=x,y=y,boxlen=len_size
                                            ,ref_matrix=ref_structure
                                            ,refval=wall_val
                                        )

                    circle_dCondt(
                        C=Free_Ca,pos=pos
                        ,const=[t,dt,dx,dy,R,dR,radii]
                        ,D=D_Ca_cyto
                    )

                    circle_dAdt(
                        A_free= Free_annexin
                        ,A_bound=Bound_annexin
                        ,C=Free_Ca
                        ,C_bound=Bound_Ca
                        ,pos=pos 
                        ,const=[t,dt,dx,dy,k1,k2,R,dR,radii]
                        ,D=D_Annexin_cyto
                                )

                    #if c_in < Free_Ca[t+1][y][x] <= c_in:
                        #Free_Ca[t+1][y][x] = c_in

                    #if Free_Ca[t+1][y][x] >= c_pump + c_in:
                    if Free_Ca[t+1][y][x] > c_in:
                        if Free_Ca[t+1][y][x] - c_pump <= c_in:
                            Free_Ca[t+1][y][x] = c_in
                        else:
                            Free_Ca[t+1][y][x] += -c_pump # the pumping mechanism
                                                #, for only inside the cell


                    if t == 0:
                        i += 1 #count points in cell.
        if t >= close_time and open_hole==True:
            if ref_bakteria == "":
                open_close_membrane(
                    Grid=ref_structure
                    ,Radius=R, dRadius=dR
                    ,offsets=[x0,y0],holesize=holesize
                    ,open_val=open_val
                    ,wall_val=wall_val
                    ,open_wall=False
                                    )
                open_hole = False
                print(f"wall closure time t={t}")
            else:
                print(" \n \n \n Need to code the other ref structure -- closing program \n \n \n")
                exit()
                
    return ref_structure,Free_Ca,Free_annexin,Bound_annexin,Bound_Ca
    

if __name__ == "__main__":
    pass
