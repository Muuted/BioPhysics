
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
    fig_folder_path =  fig_save_path + f"simtime={Real_sim_time}\\"
    fig_name = f"Simulation_data_simtime={Real_sim_time}.pkl"

    df_sim = pd.read_pickle(fig_folder_path+fig_name)
    print(df_sim.info())

    
    Free_Ca = df_sim['Free Calcium'][0]
    Free_annexin = df_sim['Free Annexins'][0]
    Bound_annexin = df_sim['Bound Annexins'][0]
    sim_ring_data_Ca = df_sim['Ring sum list Ca'][3]
    sim_ring_data_Ann = df_sim['Ring sum list Annexin'][3]
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

    """
    sim_shorten_list_Ca_exposure = np.zeros(shape=(exp_data_shape_t,exp_data_shapeX))
    sim_shorten_list_Ann_exposure = np.zeros(shape=(exp_data_shape_t,exp_data_shapeX))

    i = 0
    norm = 0
    for i in range(exp_data_shape_t):
        
        for x in range(exp_data_shapeX):
            if i < exp_data_shape_t-1:
                for t in range(int(time_check_vec[i]),int(time_check_vec[i+1])):
                    sim_shorten_list_Ann_exposure[i][x] += sim_ring_data_Ca[t][x]
                    sim_shorten_list_Ann_exposure[i][x] += sim_ring_data_Ann[t][x]
                
            if i == exp_data_shape_t:
                t = sim_T_tot
                sim_shorten_list_Ann_exposure[i][x] += sim_ring_data_Ca[t][x]
                sim_shorten_list_Ann_exposure[i][x] += sim_ring_data_Ann[t][x]
                
    """    
            
    
        
    fig, ax = plt.subplots(2,2)
    cmap_type = "RdBu"

    vmin_val_Ca , vmax_val_Ca = 1e-7 , 1e-5
    vmin_val_Ann, vmax_val_Ann = 1e-4 , 1e-5

    for i in range(exp_data_shape_t-1):
        fig.canvas.manager.window.showMaximized()
        t = int(time_check_vec[i])
        t_show = round(t*sim_dt,3)
        T_final = round(time_check_vec[len(time_check_vec)-1]*sim_dt,3)

        ax[0,0].plot(vec,sim_shorten_list_Ca[i]/(max(sim_shorten_list_Ca[i])),label="simulation")
        ax[0,0].plot(vec,real_data_Ca.loc[i]/(max(real_data_Ca.loc[i])),label="Experiment")
        ax[0,0].set_title(f"Calcium rings t={t_show}s of {T_final}")
        ax[0,0].legend()

        ax[0,1].plot(vec,sim_shorten_list_Ann[i]/(max(sim_shorten_list_Ann[i])),label="simulation")
        ax[0,1].plot(vec,real_data_Ann.loc[i]/(max(real_data_Ann.loc[i])),label="Experiment")
        ax[0,1].set_title(f"Annexin rings t={t_show}s of {T_final}")
        ax[0,1].legend()


        pos0 = ax[1,0].matshow(Free_Ca[t],cmap=cmap_type,vmin=vmin_val_Ca,vmax=vmax_val_Ca)
        ax[1,0].set_title("free Ca")

        ToT_ann =Bound_annexin[t] + Free_annexin[t]
        pos0 = ax[1,1].matshow(ToT_ann ,cmap=cmap_type,vmin=vmin_val_Ann,vmax=vmax_val_Ann)
        ax[1,1].set_title("free annexin")
        
        if i == 0:
            fig.colorbar(pos0,ax=ax[1,0],shrink=0.7)
            fig.colorbar(pos0,ax=ax[1,1],shrink=0.7)
            plt.show()
            fig, ax = plt.subplots(2,2)
            fig.colorbar(pos0,ax=ax[1,0],shrink=0.7)
            fig.colorbar(pos0,ax=ax[1,1],shrink=0.7)
        
        plt.draw()
        plt.pause(0.01)
        ax[0,0].clear()
        ax[0,1].clear()
        ax[1,0].clear()
        ax[1,1].clear()

        
            
    


    
    """
    exposure_compare = False
    if exposure_compare == False:
        fig1 = plt.figure()
        for i in range(exp_data_shape_t):

            plt.plot(vec,sim_shorten_list_Ca[i]/(max(sim_shorten_list_Ca[i])),label="simulation Ca")
            plt.plot(vec,real_data_Ca.loc[i]/(max(real_data_Ca.loc[i])),label="Experiment Ca")
            
            plt.title(f"Calcium : {i} of {exp_data_shape_t}")
            plt.legend()

            plt.draw()
            plt.pause(0.1)
            plt.clf()

        
        fig2 = plt.figure()

        for i in range(exp_data_shape_t):

            plt.plot(vec,sim_shorten_list_Ann[i]/(max(sim_shorten_list_Ann[i])),label="simulation Ca")
            plt.plot(vec,real_data_Ann.loc[i]/(max(real_data_Ann.loc[i])),label="Experiment Ca")
            
            plt.title(f"Annexins : {i} of {exp_data_shape_t}")
            plt.legend()

            plt.draw()
            plt.pause(0.1)
            plt.clf()

    if exposure_compare == True:
        fig1 = plt.figure()
        for i in range(exp_data_shape_t):

            plt.plot(vec,sim_shorten_list_Ca_exposure[i]/(max(sim_shorten_list_Ca_exposure[i])),label="simulation Ca")
            plt.plot(vec,real_data_Ca.loc[i]/(max(real_data_Ca.loc[i])),label="Experiment Ca")
            
            plt.title(f"Calcium exposure: {i} of {exp_data_shape_t}")
            plt.legend()

            plt.draw()
            plt.pause(0.1)
            plt.clf()
            exit()

        
        fig2 = plt.figure()

        for i in range(exp_data_shape_t):

            plt.plot(vec,sim_shorten_list_Ann_exposure[i]/(max(sim_shorten_list_Ann_exposure[i])),label="simulation Ca")
            plt.plot(vec,real_data_Ann.loc[i]/(max(real_data_Ann.loc[i])),label="Experiment Ca")
            
            plt.title(f"Annexins exposure : {i} of {exp_data_shape_t}")
            plt.legend()

            plt.draw()
            plt.pause(0.1)
            plt.clf()"""

if __name__ == "__main__":
    main_compare()
    print("\n \n \n ----------- DONE ------------- \n \n ")