import os
import matplotlib.pyplot as plt
import numpy as np
import time as tm
from Constants import constants
import pandas as pd
from Data_extraction_funcs import *

from Make_movie import *

def main_compare(
        Real_sim_time
        ,fig_save_path
        ,fig_folder_path
        ,video_save_path
        ,df_name
        ):
    Real_time_steps_data = 235
    #Real_sim_time = 120

    """   The simulation data loaded  """

    #fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
    
    #fig_folder_path =  fig_save_path + f"Cell structure simtime={Real_sim_time}\\"
    #fig_name = f"Cell structure Simulation_data_simtime={Real_sim_time}.pkl"
    
    
    #fig_folder_path =  fig_save_path + f"simtime={Real_sim_time}\\"   
    #video_save_path = fig_folder_path + f"video_folder\\"
    #df_name = f"Simulation_data_simtime={Real_sim_time}.pkl"

    df_sim = pd.read_pickle(fig_folder_path + df_name)
    print(df_sim.info())

    
    Free_Ca = df_sim['Free Calcium'][0]
    Free_annexin = df_sim['Free Annexins'][0]
    Bound_annexin = df_sim['Bound Annexins'][0]

    k1 = df_sim['k1'][0]
    k2 = df_sim['k2'][0]
    
    sim_T_tot, Ysize, Xsize = np.shape(Free_Ca)
    
    sim_dt = df_sim['dt'][0]
    
    #sim_T_tot = int(df_sim['time steps'][0])
    
    hole_closure_time = df_sim['hole closure time'][0]

    for i in range(5):
        a = df_sim['Ring sum list Ca'][i]
        if isinstance(a , np.ndarray):
            sim_ring_data_Ca = df_sim['Ring sum list Ca'][i]
            sim_ring_data_Ann = df_sim['Ring sum list Annexin'][i]
            print(f"i = {i}")
            break
    
    
       
    print(f"shape of simulated ring sums = {np.shape(sim_ring_data_Ca)}")

    """   Experimental data loaded  """

    data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"
    Ca_data = "Ring_sum_data_Ca_filenum4.txt"
    Annexin_data = "Ring_sum_data_Annexins_filenum4.txt"

    data_Ca = pd.read_csv(data_path + Ca_data)
    data_Ann= pd.read_csv(data_path + Annexin_data)
    
    exp_data_shape_t, exp_data_shapeX = np.shape(data_Ca)

    real_data_Ca = np.zeros(shape=(exp_data_shape_t,exp_data_shapeX))
    real_data_Ann = np.zeros(shape=(exp_data_shape_t,exp_data_shapeX))

    for t in range(exp_data_shape_t):
        for R in range(exp_data_shapeX):
            real_data_Ca[t][R] = data_Ca.loc[t][R]
            real_data_Ann[t][R] = data_Ann.loc[t][R]
    """
    First thing is that we need to make a list for the sim data, that has as many time rows as
    our experimental data set.
    """


    time_check_vec = np.linspace(0,sim_T_tot,exp_data_shape_t)

    animate_Ca = np.zeros(shape=(exp_data_shape_t,exp_data_shapeX))
    animate_Ann = np.zeros(shape=(exp_data_shape_t,exp_data_shapeX))


    
    vec = np.linspace(0,exp_data_shapeX,exp_data_shapeX)    

    max_ca_sim = np.max(np.max(sim_ring_data_Ca))
    max_ann_sim = np.max(np.max(sim_ring_data_Ann))

    max_ca_data =np.max( np.max(real_data_Ca))
    max_ann_data = np.max(np.max(real_data_Ann))

        
    scaling_factor_ca = max_ca_data/max_ca_sim
    scaling_factor_ann = max_ann_data/max_ann_sim
    
    for t in range(exp_data_shape_t):
        for R in range(exp_data_shapeX):
            real_data_Ca[t][R] *= 1/scaling_factor_ca
            real_data_Ann[t][R] *= 1/scaling_factor_ann


    if np.shape(sim_ring_data_Ca)[0] !=  np.shape(real_data_Ca)[0]:
        for i in range(exp_data_shape_t-1):
            t = int(time_check_vec[i])
            animate_Ca[i] = sim_ring_data_Ca[t]
            animate_Ann[i] = sim_ring_data_Ann[t]
    if np.shape(sim_ring_data_Ca)[0] ==  np.shape(real_data_Ca)[0]:
        animate_Ca = sim_ring_data_Ca.copy
        animate_Ann = sim_ring_data_Ann.copy
    else:
        print(f" we have a problem, the shapes doesnt match in Compare_data-py")
        print(f"shape of simulated ring sums = {np.shape(sim_ring_data_Ca)}")
        print(f"shape of animate ring sums = {np.shape(animate_Ca)}")
        print(f"shape of data ring sums = {np.shape(real_data_Ca)}")
        exit()

    
    fig, ax = plt.subplots(2,2)
    cmap_type = "gray" #"RdBu"
    #cmap_type = "RdBu"

    vmin_val_Ca , vmax_val_Ca = 5e-8 , 4e-4
    vmin_val_Ann, vmax_val_Ann = 1e-7 , 8e-6
    normalize = True
    data_opening_frame = 15
    end_i = exp_data_shape_t - 1 - data_opening_frame

    for i in range(end_i):
        fig.canvas.manager.window.showMaximized()
        if np.shape(sim_ring_data_Ca)[0] !=  np.shape(animate_Ca):
            t = int(time_check_vec[i])
        if np.shape(sim_ring_data_Ca)[0] !=  np.shape(animate_Ca):
            t = i
        t_matrix = int(time_check_vec[i])
        t_show = round(time_check_vec[i]*sim_dt,3)
        T_final = round(time_check_vec[len(time_check_vec)-1]*sim_dt,3)
        
        j = i + data_opening_frame  # getting the opening of the hole
                                    # in the experiment to match the simulation


        ax[0,0].plot(vec,sim_ring_data_Ca[t],label="simulation")
        ax[0,0].plot(vec,real_data_Ca[j],label="Experiment")
        ax[0,0].set_title(f"concentration Calcium rings \n "
                          +f"t={t_show}s of {T_final}"
                          +f"and k1={k1:.1e} , k2={k2:.1e}"
                          )
        ax[0,0].set_xlabel(f"Ring")
        ax[0,0].set_ylabel(r"[Ca]       "
                           , rotation='horizontal'
                           ,fontsize=15,y=0.45
                           )
        ax[0,0].set_ylim(0,max_ca_sim*1.1)
        ax[0,0].legend()

        ax[1,0].plot(vec,sim_ring_data_Ann[t],label="simulation")
        ax[1,0].plot(vec,real_data_Ann[j],label="Experiment")
        ax[1,0].set_title(
            f" Total concentration of Annexin rings \n"
            +f" time ={t_show}s of {T_final}"
            )
        ax[1,0].set_xlabel(f"Ring")
        ax[1,0].set_ylabel(r"[Ann]         "
                           , rotation='horizontal'
                           ,fontsize=15,y=0.45
                           )
        ax[1,0].set_ylim(0,max_ann_sim*1.1)
        ax[1,0].legend()


        pos0 = ax[0,1].matshow(Free_Ca[t_matrix],cmap=cmap_type
                               ,vmin=vmin_val_Ca,vmax=vmax_val_Ca
                               )
        ax[0,1].set_title("Concentration Ca")

        ToT_ann =Bound_annexin[t_matrix] + Free_annexin[t_matrix]
        pos0 = ax[1,1].matshow(ToT_ann ,cmap=cmap_type
                               ,vmin=vmin_val_Ann,vmax=vmax_val_Ann
                               )
        ax[1,1].set_title("total concentration annexins")
        
        if i == 0:
            fig.colorbar(pos0,ax=ax[0,1],shrink=0.7)
            fig.colorbar(pos0,ax=ax[1,1],shrink=0.7)
            #plt.show()
            #fig, ax = plt.subplots(2,2)
            #fig.colorbar(pos0,ax=ax[0,1],shrink=0.7)
            #fig.colorbar(pos0,ax=ax[1,1],shrink=0.7)
        
        plt.draw()
        

        if t < 20 :
            plt.pause(0.05)
        else:
            plt.pause(0.05)
        
        if not os.path.exists(video_save_path):
            os.makedirs(video_save_path)
        
        fig.savefig(video_save_path + f"{i}")
        plt.pause(0.01)

        ax[0,0].clear()
        ax[0,1].clear()
        ax[1,0].clear()
        ax[1,1].clear()

    
        
            


if __name__ == "__main__":
    Real_sim_time = 60
    data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"
    ref_struct_name_cell = "ref_struct_from_Ca_filenum4.txt"

    fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
    fig_folder_path =  fig_save_path + f"Cell structure simtime={Real_sim_time}\\"
    video_save_path = fig_folder_path + f"video_folder\\"     
    fig_name_df = f"Cell structure Simulation_data_simtime={Real_sim_time}.pkl"
    time1 = tm.time()
    main_compare(
        Real_sim_time=Real_sim_time
        ,fig_save_path=fig_save_path
        ,fig_folder_path=fig_folder_path
        ,video_save_path=video_save_path
        ,df_name=fig_name_df
        )
    
    #fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
    #fig_folder_path =  fig_save_path + f"simtime={Real_sim_time}\\"   
    #video_save_path = fig_folder_path + f"video_folder\\" 
    Make_video2(
            output_path=fig_folder_path
            ,input_path=video_save_path
            ,video_name= "movie.avi"
            ,fps= 8
        )
    time2 = tm.time()
    print(f"it took {(time2-time1)/60} s")
    print("\n \n \n ----------- DONE ------------- \n \n ")