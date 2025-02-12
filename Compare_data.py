import os
import matplotlib.pyplot as plt
import matplotlib
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
        ,data_path
        ,Ca_data
        ,Annexin_data
        ,opening_frame
        ):
    Real_time_steps_data = 235

    """   The simulation data loaded  """
    df_sim = pd.read_pickle(fig_folder_path + df_name)
    print(df_sim.info())

    
    Free_Ca = df_sim['Free Calcium'][0]
    Free_annexin = df_sim['Free Annexins'][0]
    Bound_annexin = df_sim['Bound Annexins'][0]

    k1 = df_sim['k1'][0]
    k2 = df_sim['k2'][0]
    
    sim_T_tot, Ysize, Xsize = np.shape(Free_Ca)
    
    sim_dt = df_sim['dt'][0]
    
       
    hole_closure_time = df_sim['hole closure time'][0]

    for i in range(5):
        a = df_sim['Ring sum list Ca'][i]
        if isinstance(a , np.ndarray):
            sim_ring_data_Ca = df_sim['Ring sum list Ca'][i]
            sim_ring_data_Ann = df_sim['Ring sum list Annexin'][i]
            print(f"i = {i}")
            break


    Tsize, ringsize = np.shape(sim_ring_data_Ann)
    background_noise = [sim_ring_data_Ann[0][i] for i in range(ringsize)]
    for t in range(Tsize):
        for ring in range(ringsize):
            pass #sim_ring_data_Ann[t][ring] += - background_noise[ring]
    
       
    print(f"minimum is = {min([min(sim_ring_data_Ann[i]) for i in range(Tsize) ])}")

    print(f"shape of simulated ring sums = {np.shape(sim_ring_data_Ca)}")

    """   Experimental data loaded  """

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

    """
    corr_ca, corr_Ann  = [],[]
    scale_ca = 0.83#1 - sim_ring_data_Ca[0][0]/real_data_Ca[0][0]
    scale_ann = 0.83#1 - sim_ring_data_Ann[0][0]/real_data_Ann[0][0]

    print(f"scale_ca={scale_ca} \n scale_ann={scale_ann}")
    for t in range(np.shape(sim_ring_data_Ca)[0]):
        for ring in range(np.shape(sim_ring_data_Ca)[1]):
            if t == 0:
                corr_ca.append(sim_ring_data_Ca[t][ring]*scale_ca)
                corr_Ann.append(sim_ring_data_Ann[t][ring]*scale_ann)
            sim_ring_data_Ca[t][ring] -= corr_ca[ring]
            sim_ring_data_Ann[t][ring] -= corr_Ann[ring]
    
    """

    vec = np.linspace(
        start=0
        ,stop=exp_data_shapeX
        ,num=exp_data_shapeX)    

    max_ca_sim = np.max(np.max(sim_ring_data_Ca))
    max_ann_sim = np.max(np.max(sim_ring_data_Ann))

    max_ca_data = np.max( np.max(real_data_Ca))
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
    cmap_type = "gist_ncar"
    #cmap_type = "nipy_spectral"
    #cmap_type = "gnuplot"
    
    cmap_type = matplotlib.cm.get_cmap('gist_ncar')
    #cmap_type.set_under('black')
    
    vmin_val_Ca , vmax_val_Ca = 1e-7 , 6e-4
    vmin_val_Ann, vmax_val_Ann = 2e-6 , 8e-6
    normalize = True
    
    end_i = exp_data_shape_t - 1 - opening_frame - 1

    min_sim_ann = min([min(sim_ring_data_Ann[t]) for t in range(Tsize-1)])
    
    
    for i in range(end_i):
        fig.canvas.manager.window.showMaximized()
        if np.shape(sim_ring_data_Ca)[0] !=  np.shape(animate_Ca):
            t = int(time_check_vec[i])
        if np.shape(sim_ring_data_Ca)[0] !=  np.shape(animate_Ca):
            t = i
        t_matrix = int(time_check_vec[i])
        t_show = round(time_check_vec[i]*sim_dt,3)
        T_final = round(time_check_vec[len(time_check_vec)-1]*sim_dt,3)
        
        j = i + opening_frame - 1  # getting the opening of the hole
                                    # in the experiment to match the simulation

        fig.suptitle(
            f"t={t_show}s of {T_final} for all the plots"
            ,fontsize=25
            )
        

        ax[0,0].plot(vec,sim_ring_data_Ca[t]
                     ,label="simulation (Molar)"
                     ,linestyle="-"
                     ,marker = "o"
                     )
        ax[0,0].plot(vec,real_data_Ca[j]
                     ,label="Experiment"
                     ,linestyle="--"
                     ,marker = "*"
                     )
        ax[0,0].set_title(f"     Calcium rings average value "
                          #+"\n "+f"t={t_show}s of {T_final}"
                          #+f"and k1={k1:.1e} , k2={k2:.1e}"
                          ,fontsize = 15
                          )
        ax[0,0].set_xlabel(f"Ring"
                           ,fontsize = 15
                           ,x=1
                           )
        ax[0,0].set_ylabel(r"[Ca]       "
                           ,rotation='horizontal'
                           ,fontsize=15,y=0.45
                           )
        ax[0,0].set_ylim(np.min(sim_ring_data_Ca),max_ca_sim*1.1)
        ax[0,0].legend()

        ax[1,0].plot(vec,sim_ring_data_Ann[t]
                     ,label="simulation (Molar)"
                     ,linestyle="-"
                     ,marker = "o"
                     )
        ax[1,0].plot(vec,real_data_Ann[j]
                     ,label="Experiment"
                     ,linestyle="--"
                     ,marker = "*"
                    )
        ax[1,0].set_title(
            f"Total Annexins rings average value "
            #+"\n"+f"time={t_show}s of {T_final}"
            ,fontsize = 15
            ,y=0.98
            )
        ax[1,0].set_xlabel(f"Ring"
                           ,fontsize = 15
                           ,x=1
                           )
        ax[1,0].set_ylabel(r"[Ann]         "
                           , rotation='horizontal'
                           ,fontsize=15
                           ,y=0.45
                           )
        ax[1,0].set_ylim(min_sim_ann,max_ann_sim*1.1)
        ax[1,0].legend()


        pos0 = ax[0,1].matshow(Free_Ca[t_matrix]
                               ,cmap=cmap_type
                               ,vmin = vmin_val_Ca , vmax = vmax_val_Ca
                               )       
        
        ax[0,1].set_title("Simulation [Ca] distribution (Molar) "
                          ,fontsize = 15
                          )

        ToT_ann =Bound_annexin[t_matrix] + Free_annexin[t_matrix]
        pos1 = ax[1,1].matshow(ToT_ann ,cmap=cmap_type
                               ,vmin = vmin_val_Ann , vmax = vmax_val_Ann
                               )
        ax[1,1].set_title("total [annexins] distribution in (Molar)"
                          ,fontsize = 15
                          )
        

        if i == 0:
            cbformat0 = matplotlib.ticker.ScalarFormatter()
            cbformat0.set_powerlimits((0,0)) 
            
            fig.colorbar(
                pos0,ax=ax[0,1]#,shrink=1
                ,format = cbformat0 #"%.1e"
                ,label="Molar"
                #,size=15
                    )
            
            
            cbformat1 = matplotlib.ticker.ScalarFormatter()   # create the formatter
            cbformat1.set_powerlimits((0,0))

            fig.colorbar(
                pos1,ax=ax[1,1]#,shrink=1
                ,format = "%.1e"#cbformat1
                ,label="Molar"
                #,size=15
                         )
            #plt.show()
            #fig, ax = plt.subplots(2,2)
            #fig.colorbar(pos0,ax=ax[0,1],shrink=0.7)
            #fig.colorbar(pos0,ax=ax[1,1],shrink=0.7)
        
        plt.draw()
        

        if t < 20 :
            plt.pause(0.01)
        else:
            plt.pause(0.01)
        
        if not os.path.exists(video_save_path):
            os.makedirs(video_save_path)
        
        fig.savefig(video_save_path + f"{i}")
        plt.pause(0.01)

        ax[0,0].clear()
        ax[0,1].clear()
        ax[1,0].clear()
        ax[1,1].clear()

    
        
            


if __name__ == "__main__":
    const_list = constants()
    c_in ,c_out, D_Ca_cyto, T_tot, len_size, dx, dy, k1, k2 = const_list[0:9]
    c_in_annexin ,bound_annexin_start ,D_Annexin_cyto = const_list[9:12]
    dt ,close_time, c_pump, holesize ,dR ,R ,x0 ,y0 = const_list[12:20]
    wall_val ,inside_val ,outside_val ,open_val = const_list[20:24]
    Real_sim_time, real_close_time = const_list[24:26]
    ref_struct_name_cell ,fig_save_path = const_list[26:28]
    fig_folder_path ,video_save_path ,fig_name_df, data_path = const_list[28:32]
    Ca_data_experiment ,Annexin_data_experiment = const_list[32:34]
    frame_open ,frame_close = const_list[34:36]
    
    time1 = tm.time()
    main_compare(
        Real_sim_time=Real_sim_time
        ,fig_save_path=fig_save_path
        ,fig_folder_path=fig_folder_path
        ,video_save_path=video_save_path
        ,df_name=fig_name_df
        ,data_path = data_path
        ,Ca_data = Ca_data_experiment
        ,Annexin_data = Annexin_data_experiment
        ,opening_frame=frame_open
        )
    
    
    Make_video2(
            output_path=fig_folder_path
            ,input_path=video_save_path
            ,video_name= "movie.avi"
            ,fps= 8
        )
    time2 = tm.time()
    print(f"it took {(time2-time1)/60} s")
    print("\n \n \n ----------- DONE ------------- \n \n ")