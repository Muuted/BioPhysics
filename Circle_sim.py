import matplotlib.pyplot as plt
import numpy as np
from Constants import *
from Circle_funcs import *
from Data_extraction_funcs import *
from Constants import constants
import time as tm
import os
import pandas as pd
from Real_data import main_ring_summing

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
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val,Real_sim_time, real_close_time = constants()
    
    data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"
    ref_struct_name = "ref_struct__filenum4.txt"

    save_data = True
    c_pump  = 0 

    Sim_data_list = main_circle_sim(
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,open_hole=True
                    )
    ref_structure,Free_Ca,Free_annexin,Bound_annexin,Bound_Ca = Sim_data_list

    """data_list = Ring_sum(
        ref_Grid=Free_Ca[T_tot-1]
        ,offsets=[x0,y0]
        ,Radius=R   ,dRadius=dR
        ,num_of_rings=3
    )
    
    Ring_sums = data_list[0]
    Ring_radius = data_list[1]
    conc_removed_grid = data_list[2]
    visual_rings = data_list[3]"""

    sumfree,sumbound,sumtot = sum_annexin(
        A_free=Free_annexin
        ,A_bound=Bound_annexin
    )

    time_vec = np.linspace(0,Real_sim_time,len(sumfree))
    
    print(f"len(sumfree)={len(sumfree)}")
    print(f"len(sumbound)={len(sumbound)}")
    print(f"len(sumtot)={len(sumtot)}")
    
    fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"

    cmap_type = "hot"
    #cmap_type = "viridis"
    cmap_type = "RdBu"

    fig_ann_time, ax = plt.subplots(3,1)
    fig_ann_time.canvas.manager.window.showMaximized()

    real_time_close_time = int(close_time*dt) 
    ax[0].plot(time_vec,sumtot,label="Total annexins")
    if len(sumtot) >= real_time_close_time:
        ax[0].vlines(x=np.ceil(real_time_close_time),ymin=min(sumtot),ymax=max(sumtot),colors='k',label=f"holeclose t={close_time}")
    ax[0].legend()
    ax[0].set_title(f"total sum of Annexins \n dt = {dt}")
    
    
    #ax[1].hlines(y=sumfree[len(sumfree)-1],xmin=0,xmax=len(sumfree),label=f"{sumfree[len(sumfree)-1]}")
    ax[1].plot(time_vec,sumfree,'k',label="free annnexins")
    ax[1].set_title(f"Free Annexins, holeclosed={close_time}")
    ax[1].legend()
    
    
    #ax[2].hlines(y=sumbound[len(sumbound)-1],xmin=0,xmax=len(sumbound),label=f"A_b_eval={sumbound[len(sumbound)-1]}")
    #ax[2].hlines(y=A_b_init*i,xmin=0,xmax=len(sumbound),label=f"A_b_init={A_b_init*i}")
    ax[2].plot(time_vec,sumbound,'k',label="Bound Annexins")
    ax[2].set_title(f"bound Annexins, holeclosed={close_time}")
    ax[2].legend()
    
    
    
    

    fig_ref_structure,axs = plt.subplots()

    fig_ref_structure.canvas.manager.window.showMaximized()
    pos = axs.matshow(ref_structure,cmap=cmap_type)
    axs.set_title("Reference structure")
    fig_ref_structure.colorbar(pos,ax=axs,shrink=0.7)

    """ 
    plt.figure()
    plt.plot(Ring_radius,Ring_sums,'-.')
    plt.xlabel("number of pixels from opening")
    plt.title("Ring sums try, normalized by Area")
    """

    fig_mat_ann, ax = plt.subplots(1,2)
    fig_mat_ann.canvas.manager.window.showMaximized()

    pos0 = ax[0].matshow(Free_annexin[T_tot-1],cmap=cmap_type)
    ax[0].set_title("free annexin end")
    fig_mat_ann.colorbar(pos0,ax=ax[0],shrink=0.7)

    pos1 = ax[1].matshow(Bound_annexin[T_tot-1], cmap =cmap_type)
    ax[1].set_title("bound annexin")
    fig_mat_ann.colorbar(pos1, ax=ax[1],shrink=0.7)
    

    removed_conc_free_Ca = remove_conc_outside(
        ref_grid=ref_structure
        ,grid=Free_Ca[T_tot-1]
        ,outside_val=outside_val
        )


    fig_mat_Ca, ax = plt.subplots(1,2)

    fig_mat_Ca.canvas.manager.window.showMaximized()
    pos3 = ax[0].matshow(Free_Ca[T_tot-1],cmap=cmap_type)
    ax[0].set_title("free ca end")
    fig_mat_Ca.colorbar(pos3,ax=ax[0],shrink=0.7)
    
    pos4 = ax[1].matshow(Bound_Ca[T_tot-1],cmap=cmap_type)
    ax[1].set_title("bound Ca end")
    fig_mat_Ca.colorbar(pos4,ax=ax[1],shrink=0.7)
    
    


    conc_over_time_Free_Ca = sum_in_cell(ref_Grid=ref_structure
                                 ,Matrix_Free=Free_Ca
                                 ,Matrix_Bound=Bound_Ca
                                 ,inside_val=inside_val
                                 )
    
    fig_conc_Ca_time = plt.figure()
    fig_conc_Ca_time.canvas.manager.window.showMaximized()
    plt.plot(time_vec,conc_over_time_Free_Ca/max(conc_over_time_Free_Ca))
    plt.title("Concentration free Ca over time inside the cell")
    plt.xlabel("time [seconds]")
    plt.ylabel("[Ca]")  
    
    
    if save_data == True:
        df = pd.DataFrame({
            'Free Calcium': [Free_Ca],
            'Bound Calcium': [Bound_Ca],
            'Free Annexins': [Free_annexin],
            'Bound Annexins': [Bound_annexin],
            'Reference Matrix': [ref_structure],
            'hole position X': x0 + R,
            'hole position Y': y0,
            'Sim time (s)': Real_sim_time,
            'time steps': T_tot,
            'k1': k1,
            'k2' :k2,
            'dt': dt,
            'free annexin start': c_in_annexin,
            'bound annexin start': bound_annexin_start,
            'D_Ca': D_Ca_cyto,
            'D_annexin': D_Annexin_cyto,
            'hole closure time': close_time,
            'c_pump': c_pump,
            'hole size': holesize,
            'dx':dx,
            'dy':dy,
            'side length':len_size,
            'inside val': inside_val,
            'outside val': outside_val,
            'wall val': wall_val
                        })

        print(df.info())

        fig_folder_path =  fig_save_path + f"simtime={Real_sim_time}\\"
        if not os.path.exists(fig_folder_path):
            os.makedirs(fig_folder_path)

        df_name = f"Simulation_data_simtime={Real_sim_time}.pkl"
        df.to_pickle(fig_folder_path + df_name)


        plt.show(block=False)
        plt.pause(10)


        fig_name = f"Annexins over time Realsimtime={Real_sim_time}s"
        fig_ann_time.savefig(fig_folder_path + fig_name)

        fig_name = f"ref_structure Realsimtime={Real_sim_time}s"
        fig_ref_structure.savefig(fig_folder_path + fig_name)

        fig_name = f"mat Annexin Realsimtime={Real_sim_time}s"
        fig_mat_ann.savefig(fig_folder_path + fig_name)

        fig_name = f"Mat Ca Realsimtime={Real_sim_time}s"
        fig_mat_Ca.savefig(fig_folder_path + fig_name)

        fig_name = f"[Ca] over time Realsimtime={Real_sim_time}s"
        fig_conc_Ca_time.savefig(fig_folder_path + fig_name)

        plt.close("all")

        
        main_ring_summing(
        fig_save_path = fig_save_path
        ,fig_folder_path = fig_folder_path
        ,df_name = df_name
        ,hole_pos=""#[74,40]
                        )
        