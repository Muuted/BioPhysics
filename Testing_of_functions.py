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

if __name__ =="__main__":
    test_init_circle()