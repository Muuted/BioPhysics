
import matplotlib.pyplot as plt
import numpy as np

from Constants import constants
import pandas as pd
from Data_extraction_funcs import *

def main_compare():
    Real_time_steps_data = 235
    Real_sim_time = 120

    """   The simulation data loaded  """

    fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
    fig_folder_path =  fig_save_path + f"Cell structure simtime={Real_sim_time}\\"
    fig_folder_path =  fig_save_path + f"simtime={Real_sim_time}\\"

    #fig_name = f"Cell structure Simulation_data_simtime={Real_sim_time}.pkl"
    fig_name = f"Simulation_data_simtime={Real_sim_time}.pkl"

    df_sim = pd.read_pickle(fig_folder_path + fig_name)
    print(df_sim.info())

    
    Free_Ca = df_sim['Free Calcium'][0]
    Free_annexin = df_sim['Free Annexins'][0]
    Bound_annexin = df_sim['Bound Annexins'][0]

    for i in range(5):
        a = df_sim['Ring sum list Ca'][i]
        if isinstance(a , np.ndarray):
            sim_ring_data_Ca = df_sim['Ring sum list Ca'][i]
            sim_ring_data_Ann = df_sim['Ring sum list Annexin'][i]
            print(f"i = {i}")
            break
    
    sim_dt = df_sim['dt'][0]
    sim_real_time = df_sim['Sim time (s)'][0]
    sim_T_tot = int(df_sim['time steps'][0])
    
        


    """   Experimental data loaded  """

    data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"
    Ca_data = "Ring_sum_data_Ca_filenum4.txt"
    Annexin_data = "Ring_sum_data_Annexins_filenum4.txt"
    real_data_Ca = pd.read_csv(data_path + Ca_data)
    real_data_Ann= pd.read_csv(data_path + Annexin_data)
    


    """
    First thing is that we need to make a list for the sim data, that has as many time rows as
    our experimental data set.
    """
    exp_data_shape_t, exp_data_shapeX = np.shape(real_data_Ca)

    sim_shorten_list_Ca = np.zeros(shape=(exp_data_shape_t,exp_data_shapeX))
    sim_shorten_list_Ann = np.zeros(shape=(exp_data_shape_t,exp_data_shapeX))

    time_check_vec = np.linspace(0,sim_T_tot,exp_data_shape_t)

    
    i = 0
    for t in range(sim_T_tot):
        if t >= time_check_vec[i]:
            for x in range(exp_data_shapeX):
                sim_shorten_list_Ca[i][x] = sim_ring_data_Ca[t][x]
                sim_shorten_list_Ann[i][x] = sim_ring_data_Ann[t][x]
            i += 1
    
    
    vec = np.linspace(0,exp_data_shapeX,exp_data_shapeX)    

        
    fig, ax = plt.subplots(2,2)
    cmap_type = "RdBu"

    vmin_val_Ca , vmax_val_Ca = 1e-7 , 1e-5
    vmin_val_Ann, vmax_val_Ann = 1e-6 , 6e-5
    normalize = True
    data_opening_frame = 20
    end_i = exp_data_shape_t - 1 - data_opening_frame

    for i in range(end_i):
        fig.canvas.manager.window.showMaximized()
        t = int(time_check_vec[i])
        t_show = round(t*sim_dt,3)
        T_final = round(time_check_vec[len(time_check_vec)-1]*sim_dt,3)
        
        j = i + data_opening_frame  # getting the opening of the hole
                                    # in the experiment to match the simulation

        if normalize == True:
            norm_sim_Ca = max(sim_shorten_list_Ca[i])
            norm_sim_Ca = max(sim_ring_data_Ca[t])
            norm_data_Ca = max(real_data_Ca.loc[j])

            norm_sim_ann = max(sim_shorten_list_Ann[i])
            norm_sim_ann = max(sim_ring_data_Ann[t])
            norm_data_ann = max(real_data_Ann.loc[j])
        if normalize == False:
            norm_sim_Ca, norm_data_Ca = 1 ,1
            norm_sim_ann, norm_data_ann = 1 ,1

        ax[0,0].plot(vec,sim_ring_data_Ca[t]/norm_sim_Ca,label="simulation")
        ax[0,0].plot(vec,real_data_Ca.loc[j]/norm_data_Ca,label="Experiment")
        ax[0,0].set_title(f"Calcium rings t={t_show}s of {T_final}")
        ax[0,0].set_xlabel(f"Ring")
        ax[0,0].set_ylabel(r" $ \frac{ [Ca] }{ max([Ca]) } $ ", rotation='horizontal')
        ax[0,0].legend()

        ax[1,0].plot(vec,sim_ring_data_Ann[t]/norm_sim_ann,label="simulation")
        ax[1,0].plot(vec,real_data_Ann.loc[j]/norm_data_ann,label="Experiment")
        ax[1,0].set_title(f"Annexin rings t={t_show}s of {T_final}")
        ax[1,0].set_xlabel(f"Ring")
        ax[1,0].set_ylabel(r" $ \frac{ [Ann] }{ max([Ann]) } $ ", rotation='horizontal')
        ax[1,0].legend()


        pos0 = ax[0,1].matshow(Free_Ca[t],cmap=cmap_type,vmin=vmin_val_Ca,vmax=vmax_val_Ca)
        ax[0,1].set_title("free Ca")

        #ToT_ann =Bound_annexin[t] + Free_annexin[t]
        pos0 = ax[1,1].matshow(Free_annexin[t] ,cmap=cmap_type,vmin=vmin_val_Ann,vmax=vmax_val_Ann)
        ax[1,1].set_title("Free annexin")
        
        if i == 0:
            fig.colorbar(pos0,ax=ax[0,1],shrink=0.7)
            fig.colorbar(pos0,ax=ax[1,1],shrink=0.7)
            plt.show()
            fig, ax = plt.subplots(2,2)
            fig.colorbar(pos0,ax=ax[0,1],shrink=0.7)
            fig.colorbar(pos0,ax=ax[1,1],shrink=0.7)
        
        plt.draw()
        plt.pause(0.01)
        ax[0,0].clear()
        ax[0,1].clear()
        ax[1,0].clear()
        ax[1,1].clear()

        
            


if __name__ == "__main__":
    main_compare()
    print("\n \n \n ----------- DONE ------------- \n \n ")