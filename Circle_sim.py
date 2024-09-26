import matplotlib.pyplot as plt
import numpy as np
from Constants import *
from Circle_funcs import *

def main_circle_sim():
    open_hole = True
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
        ,c_in=c_in,c_out=c_out
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

    Bound_annexin = np.zeros(shape=(T_tot,len_size,len_size))
    Bound_ca = np.zeros(shape=(T_tot,len_size,len_size))

    #plt.matshow(ref_structure)
    #plt.title("ref closed")

    #ref_structure = 
    open_close_membrane(
        Grid=ref_structure
        ,Radius=R, dRadius=dR
        ,offsets=[x0,y0],holesize=holesize
        ,open_val=open_val
        ,wall_val=wall_val
        ,open_wall=True
    )

    #plt.matshow(ref_structure)
    #plt.title("ref open")
    #plt.show()

    ref_tot_conc_Annexin = np.sum(Free_Ca[0])

    for t in np.arange(0,T_tot-1): 
        if t%(T_tot/10) == 0:
            print(f"time={t} of {T_tot}")   
        t1, t2 = t%2, (t+1)%2
        t1, t2 = t, t+1
        for x in range(0,len_size):
            for y in range(0,len_size):
                if ref_structure[y][x] == wall_val:
                    Free_Ca[t+1][y][x] = Free_Ca[t][y][x]
                    Free_annexin[t+1][y][x] = Free_annexin[t][y][x] 
                    Bound_annexin[t+1][y][x] = Bound_annexin[t][y][x]

                if ref_structure[y][x] != wall_val:
                    radii = np.sqrt( (x-x0)**2 + (y-y0)**2 )
                    pos = cicle_boundary(x=x,y=y,boxlen=len_size
                                            ,ref_matrix=ref_structure
                                            ,refval=wall_val
                                        )

                    circle_dCondt(C=Free_Ca,pos=pos
                                        ,const=[t,dt,dx,dy,R,dR,radii]
                                        ,D_list=[D_Ca_cyto,D_Ca_water]
                                        )

                    circle_dAdt(
                        A_free= Free_annexin
                        ,A_bound=Bound_annexin
                        ,C=Free_Ca
                        ,pos=pos 
                        ,const=[t,dt,dx,dy,k1,k2,R,dR,radii]
                        ,D_list=[D_Annexin_cyto,D_Annexin_water]
                    )
                    
                    if radii < R and Free_Ca[t+1][y][x] > c_pump:
                        Free_Ca[t+1][y][x] += -c_pump # the pumping mechanism
                                                #, for only inside the cell

                if add_Ca == True:
                    if x < 3 or y < 3 or x > len_size-4 or y > len_size-4:
                        Free_Ca[t+1][y][x] = c_out 
                        #We assume infinite supply, i.e the value
                        #is allways the same the walls

                if t >= close_time and open_hole==True:
                    #ref_structure = 
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
                
        #if t%(int(T_tot*0.1)) == 0:
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

    timevec = np.linspace(0,T_tot,len(sumfree))


    plt.figure()
    plt.plot(sumtot)
    plt.title("total sum of Annexins")
    print(min(sumtot),max(sumtot))

    plt.figure()
    plt.plot(sumfree)
    plt.title("Free Annexins")
    #plt.ylim(min(sumfree)*1.1,max(sumfree)*1.1)
    #plt.ylim(3.05e-5 - max(sumbound),3.05e-5 + max(sumbound))
    print(min(sumfree),max(sumfree))

    plt.figure()
    plt.plot(sumbound)
    plt.title("bound Annexins")
    print(min(sumbound),max(sumbound))
    #plt.show()


    plt.figure()
    plt.plot(Ring_radius,Ring_sums,'-.')
    plt.xlabel("number of pixels from opening")
    plt.title("Ring sums try, normalized by num of points")

    plt.matshow(conc_removed_grid)
    plt.title("concentration outside inner cell removed")

    plt.matshow(visual_rings)
    plt.title("Visualizing the rings")

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


    plt.figure()
    plt.plot(conc_time_list,conc_list)
    plt.xlabel("timestep")
    plt.title("Ca concentration over time")


    plt.show()

if __name__ == "__main__":
    main_circle_sim()
