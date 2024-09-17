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

if __name__ =="__main__":
    test_draw_rings()