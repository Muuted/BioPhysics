import matplotlib.pyplot as plt
import numpy as np
from Constants import *
from Circle_funcs import *
from Data_extraction_funcs import *
from Constants import constants
import time as tm
import os
import pandas as pd
from Ring_sum_file import main_ring_summing
from Make_movie import Make_video2
from Compare_data import main_compare

if __name__ == "__main__":
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val,Real_sim_time, real_close_time = constants()
    
    data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"
    ref_struct_name = "ref_struct__filenum4.txt"

    save_data = True
    

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

    
    ax[0].plot(time_vec,sumtot,label="Total annexins")
    if len(sumtot) >= np.ceil(close_time*dt):
        ax[0].vlines(x=np.ceil(close_time*dt)+1,ymin=min(sumtot),ymax=max(sumtot),colors='k',label=f"holeclose t={close_time}")
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
            
        main_compare(
        Real_sim_time=Real_sim_time
        )
        fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
        fig_folder_path =  fig_save_path + f"simtime={Real_sim_time}\\"   
        video_save_path = fig_folder_path + f"video_folder\\" 
        Make_video2(
            output_path=fig_folder_path
            ,input_path=video_save_path
            ,video_name= "movie.avi"
            ,fps= 4
        )
        