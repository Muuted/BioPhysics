import matplotlib.pyplot as plt
import numpy as np
from Constants import *
from Circle_funcs import *


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

    Grid,ref_Grid = init_ref_circle(boxlen=sidelen,time=T_tot
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


def test_draw_rings():
    len_size=40
    Radius = int(len_size*0.5/2)
    dR = 2
    dRadius = int(2*Radius/10)
    x0,y0 = int(len_size/2), int(len_size/2)

    xstart = x0 + Radius + 2*dR
    xend = x0 - Radius - 2*dR

    Grid = np.zeros(shape=(len_size,len_size))
    sum = 0
    Z = 0
    for R in range(0,int(2*Radius + x0 ),dRadius):
        sum += 10    
        Z += 1 # normalizing constant
        for x in range(xstart,xend,-1):#(xstart,xend,-1):
            for y in range(len_size):
                if x < len_size-1 and y < len_size:# and x < xstart:
                    r1 = np.sqrt((x-xstart)**2 + (y-y0)**2)
                    r2 = np.sqrt((x-x0)**2 + (y-y0)**2)
                    if R-dRadius < r1 <= R:# and r < Radius:
                        if r2 < Radius:
                            Grid[y][x] = sum

        
        print(R)
    
    plt.matshow(Grid)
    plt.show()


def test_annexin_diffusion():
    k1,k2 = 0.1,0.1
    T_tot = 3000
    free_Ca = np.zeros(shape=(T_tot,len_size,len_size))
    free_annexin = np.zeros(shape=(T_tot,len_size,len_size))
    bound_annexin = np.zeros(shape=(T_tot,len_size,len_size))
    ref_grid = np.zeros(shape=(len_size,len_size))

    free_annexin[0][20][31] = 100

    ref_grid = init_ref_circle(boxlen=len_size
                               ,Radius=R,dRadius=dR,offsets=[x0,y0]
                               ,inside_val=inside_val,outside_val=outside_val,wall_val=wall_val
                               )
    open_close_membrane(Grid=ref_grid,Radius=R,dRadius=dR,offsets=[x0,y0]
                        ,holesize=holesize,open_val=open_val,wall_val=wall_val
                        ,open_wall=True
                        )

    free_Ca = init_conc(ref_grid=ref_grid
              ,time=T_tot,c_in=c_in,c_out=c_out
              ,inside_val=inside_val,outside_val=outside_val
              )
    for t in range(T_tot-1):
        if t%(T_tot*0.1)== 0:
            print(f"t={t} out of {T_tot}")
        for y in range(len_size):
            for x in range(len_size):
                if ref_grid[y][x] == wall_val:
                    free_Ca[t+1][y][x] = free_Ca[t][y][x]
                    free_annexin[t+1][y][x] = free_annexin[t][y][x] 
                    bound_annexin[t+1][y][x] = bound_annexin[t][y][x]

                if ref_grid[y][x] != wall_val:
                    radii = np.sqrt( (x-x0)**2 + (y-y0)**2 )

                    pos = cicle_boundary(x=x,y=y,boxlen=len_size
                                        ,ref_matrix=ref_grid
                                        ,refval=wall_val
                                        )
                    
                    circle_dAdt(A_free=free_annexin
                                ,A_bound=bound_annexin
                                ,C=free_Ca
                                ,pos=pos
                                ,const=[t,dt,dx,dy,k1,k2,R,dR,radii]
                                ,D_list=[D_Annexin_cyto,D_Annexin_water]
                                )
                    circle_dCondt(C=free_Ca,pos=pos
                                    ,const=[t,dt,dx,dy,R,dR,radii]
                                    ,D_list=[D_Ca_cyto,D_Ca_water]
                                )
    
    sumAfree, sumAbound, sumtot= sum_annexin(A_free=free_annexin,A_bound=bound_annexin)

    plt.figure()
    plt.plot(sumtot)
    plt.title("total sum of Annexins")
    print(min(sumtot),max(sumtot))

    plt.figure()
    plt.plot(sumAfree)
    plt.title("Free Annexins")
    print(min(sumAfree),max(sumAfree))

    plt.figure()
    plt.plot(sumAbound)
    plt.title("bound Annexins")
    print(min(sumAbound),max(sumAbound))

    #plt.show()
    
    plt.matshow(ref_grid)
    plt.title("reference grid")

    plt.matshow(free_annexin[0])
    plt.title("annexin at t=0")
    plt.matshow(free_annexin[T_tot-1])
    plt.title("Annexin at end" + "\n"+"Does anything look weird on the boundaries?")
    plt.show()
    

if __name__ =="__main__":
    #test_draw_rings()
    test_annexin_diffusion()
    