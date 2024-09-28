import matplotlib.pyplot as plt
import numpy as np
from Constants import *
from Circle_funcs import *
from Data_extraction_funcs import *
import time as tm


def main_circle_sim():
    open_hole = False
    add_Ca = True
    conc_list = []
    conc_time_list = []
    
    # Creation of our Grids
    ref_structure = init_ref_circle(
        boxlen=len_size
        ,Radius=R,dRadius=dR ,offsets=[x0,y0]
        ,inside_val=inside_val
        ,outside_val=outside_val
        ,wall_val=wall_val
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
        if t%(T_tot/10) == 0:
            print(f"time={t} of {T_tot}")   
        t1, t2 = t%2, (t+1)%2
        t1, t2 = t, t+1
        for x in range(0,len_size):
            for y in range(0,len_size):
                #Bound_Ca[t+1][y][x] += Bound_Ca[t][y][x]
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
                    
                    if Free_Ca[t+1][y][x] > c_in:
                        Free_Ca[t+1][y][x] += -c_pump # the pumping mechanism
                                                    #, for only inside the cell

                    if t == 0:
                        i += 1 #count points in cell.
        if t >= close_time and open_hole==True:
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
                
        conc_list.append(np.sum(Free_Ca[t]))
        conc_time_list.append(t)

    data_list = Ring_sum(
        ref_Grid=Free_Ca[T_tot-1]
        ,offsets=[x0,y0]
        ,Radius=R   ,dRadius=dR
        ,num_of_rings=3
    )
    Ring_sums = data_list[0]
    Ring_radius = data_list[1]
    conc_removed_grid = data_list[2]
    visual_rings = data_list[3]

    sumfree,sumbound,sumtot = sum_annexin(
        A_free=Free_annexin
        ,A_bound=Bound_annexin
    )

    plt.figure()
    plt.plot(sumtot,label="Total annexins")
    plt.title("total sum of Annexins")
    print(min(sumtot),max(sumtot))

    plt.figure()
    plt.hlines(y=sumfree[len(sumfree)-1],xmin=0,xmax=len(sumfree),label=f"{sumfree[len(sumfree)-1]}")
    plt.plot(sumfree,'k',label="free annnexins")
    plt.title("Free Annexins")
    plt.legend()

    plt.figure()
    plt.hlines(y=sumbound[len(sumbound)-1],xmin=0,xmax=len(sumbound),label=f"A_b_eval={sumbound[len(sumbound)-1]}")
    plt.hlines(y=A_b_init*i,xmin=0,xmax=len(sumbound),label=f"A_b_init={A_b_init*i}")
    plt.plot(sumbound,'k',label="Bound Annexins")
    plt.title("bound Annexins")
    plt.legend()

    plt.figure()
    plt.plot(conc_time_list,conc_list)
    plt.xlabel("timestep")
    plt.title("Ca concentration over time")

    plt.show()
    
    """ 
    plt.figure()
    plt.plot(Ring_radius,Ring_sums,'-.')
    plt.xlabel("number of pixels from opening")
    plt.title("Ring sums try, normalized by Area")

    plt.matshow(Free_annexin[0])
    plt.title("free annexin start")
    plt.matshow(Free_annexin[T_tot-1])
    plt.title("free annexin end")

    plt.matshow(Bound_annexin[T_tot-1])
    plt.title("bound annexin")
    #plt.show()



    #plt.figure()
    #plt.plot(Free_Ca[90][y0][:],label="t=90,y=y0")
    #plt.legend()

    plt.matshow(Free_Ca[T_tot-1])
    plt.title("free ca end")

    plt.show()
    """

if __name__ == "__main__":
    main_circle_sim()
    
