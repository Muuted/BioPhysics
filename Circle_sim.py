import matplotlib.pyplot as plt
import numpy as np
from Constants import *
from Circle_funcs import *
from Data_extraction_funcs import *
from Constants import constants
import time as tm


def main_circle_sim(
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,open_hole=True
        ,ref_bakteria = ""
        ,data_path=""
                    ) -> list:
    #open_hole = True
    
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

                    if c_in < Free_Ca[t+1][y][x] <= c_pump + c_in:
                        Free_Ca[t+1][y][x] += c_in
                    if Free_Ca[t+1][y][x] >= c_pump + c_in:
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
                
    return ref_structure,Free_Ca,Free_annexin,Bound_annexin,Bound_Ca
    





    

if __name__ == "__main__":
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val,Real_sim_time, real_close_time = constants()
    
    data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"
    ref_struct_name = "ref_struct__filenum4.txt"

    #T_tot = 100

    Sim_data_list = main_circle_sim(
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,open_hole=True
                    )
    ref_structure,Free_Ca,Free_annexin,Bound_annexin,Bound_Ca = Sim_data_list

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


    
    cmap_type = "hot"
    cmap_type = "viridis"

    fig, ax = plt.subplots(3,1)
    
    ax[0].plot(sumtot,label="Total annexins")
    if len(sumtot) >= close_time:
        ax[0].vlines(x=np.ceil(close_time),ymin=min(sumtot),ymax=max(sumtot),colors='k',label=f"holeclose t={close_time}")
    ax[0].legend()
    ax[0].set_title(f"total sum of Annexins \n dt = {dt}")
    
    
    #ax[1].hlines(y=sumfree[len(sumfree)-1],xmin=0,xmax=len(sumfree),label=f"{sumfree[len(sumfree)-1]}")
    ax[1].plot(sumfree,'k',label="free annnexins")
    ax[1].set_title(f"Free Annexins, holeclosed={close_time}")
    ax[1].legend()
    
    
    #ax[2].hlines(y=sumbound[len(sumbound)-1],xmin=0,xmax=len(sumbound),label=f"A_b_eval={sumbound[len(sumbound)-1]}")
    #ax[2].hlines(y=A_b_init*i,xmin=0,xmax=len(sumbound),label=f"A_b_init={A_b_init*i}")
    ax[2].plot(sumbound,'k',label="Bound Annexins")
    ax[2].set_title(f"bound Annexins, holeclosed={close_time}")
    ax[2].legend()
    
    plt.matshow(ref_structure,cmap=cmap_type)
    plt.title("Reference structure")
    plt.colorbar()

    """ 
    plt.figure()
    plt.plot(Ring_radius,Ring_sums,'-.')
    plt.xlabel("number of pixels from opening")
    plt.title("Ring sums try, normalized by Area")
    """

    fig, ax = plt.subplots(1,2)
    pos0 = ax[0].matshow(Free_annexin[T_tot-1],cmap=cmap_type)
    ax[0].set_title("free annexin end")
    #fig.colorbar(pos0,ax=ax[0])

    pos1 = ax[1].matshow(Bound_annexin[T_tot-1], cmap =cmap_type)
    ax[1].set_title("bound annexin")
    #fig.colorbar(pos1, ax=ax[1])


    removed_conc_free_Ca = remove_conc_outside(
        ref_grid=ref_structure
        ,grid=Free_Ca[T_tot-1]
        ,outside_val=outside_val
        )

    fig, ax = plt.subplots(1,2)
    
    pos3 = ax[0].matshow(Free_Ca[T_tot-1],cmap=cmap_type)
    ax[0].set_title("free ca end")
    #fig.colorbar(pos3,ax=ax[0])
    
    pos4 = ax[1].matshow(Bound_Ca[T_tot-1],cmap=cmap_type)
    ax[1].set_title("bound Ca end")
    #fig.colorbar(pos4,ax=ax[1])



    conc_over_time_Free_Ca = sum_in_cell(ref_Grid=ref_structure
                                 ,Matrix_Free=Free_Ca
                                 ,Matrix_Bound=Bound_Ca
                                 ,inside_val=inside_val
                                 )
    
    plt.figure()
    plt.plot(conc_over_time_Free_Ca/max(conc_over_time_Free_Ca))
    plt.title("Concentration free Ca over time inside the cell")

    

    
    
    plt.show()
    
    
