import matplotlib.pyplot as plt
import numpy as np

from Constants import constants
import pandas as pd
from Data_extraction_funcs import Ring_sum,Ring_sum_quick

def main_ring_summing(
        fig_save_path
        ,fig_folder_path
        ,df_name
        ,hole_pos = ""
        ,sum_quick = False
                ):
    
    A = constants()
    """
    c_in,c_out,D_Ca_cyto,T_tot,len_size,
    dx,dy,k1,k2,c_in_annexin,bound_annexin_start,
    D_Annexin_cyto,dt,close_time,c_pump,holesize,
    dR,R,x0,y0,wall_val,inside_val,outside_val,
    open_val,Real_sim_time, real_close_time
    """

    #data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"
    #ref_struct_name = "ref_struct__filenum4.txt"
    if type(hole_pos) == str:
        hole_pos = hole_pos.strip()

    if hole_pos != "":
        x0,y0 = hole_pos

    Real_time_steps_data = 234
    #Real_sim_time = 1
    save_data = True
    
    
    df = pd.read_pickle(fig_folder_path + df_name)

    print(df.info())

    # constants
    inside_val = df['inside val'][0]

    Free_Ca = df['Free Calcium'][0]
    Bound_Ca = df['Bound Calcium'][0]
    
    print(np.shape(Free_Ca))

    Free_Annexins = df['Free Annexins'][0]
    Bound_annexins = df['Bound Annexins'][0]

    ref_structure = df['Reference Matrix'][0]
    
    if hole_pos == "":
        x0 = df['hole position X'][0]
        y0 = df['hole position Y'][0]
    
    sim_time_steps = df['time steps']

    time_step_vec_data = np.linspace(0,sim_time_steps,Real_time_steps_data)
    
    if sum_quick == True:
        Ysize,Xsize = np.shape(ref_structure)
        time_step_vec_data = np.linspace(0,sim_time_steps,Real_time_steps_data)
        short_Free_Ca = np.zeros(shape=(Real_time_steps_data,Ysize,Xsize))
        short_Bound_Ca = np.zeros(shape=(Real_time_steps_data,Ysize,Xsize))
        short_Free_Annexins = np.zeros(shape=(Real_time_steps_data,Ysize,Xsize))
        short_Bound_Annexins = np.zeros(shape=(Real_time_steps_data,Ysize,Xsize))

        for frame in range( len(time_step_vec_data)-1):
            t = int(time_step_vec_data[frame])
            short_Free_Ca[frame] = Free_Ca[t]
            short_Bound_Ca[frame] = Bound_Ca[t]
            short_Free_Annexins[frame] = Free_Annexins[t]
            short_Bound_Annexins[frame] = Bound_annexins[t]
            
        Ring_sum_list_Ca, Ring_sum_list_Annexin = Ring_sum(
            ref_grid = ref_structure
            ,sim_grid_free_Ca = short_Free_Ca
            ,sim_grid_bound_Ca = short_Bound_Ca
            ,sim_grid_free_Annexin = short_Free_Annexins
            ,sim_grid_bound_Annexin = short_Bound_Annexins
            ,hole_pos = [x0,y0]
            ,num_of_rings = 10
            ,inside_val = inside_val
        )
    else:
        Ring_sum_list_Ca, Ring_sum_list_Annexin = Ring_sum(
            ref_grid = ref_structure
            ,sim_grid_free_Ca = Free_Ca
            ,sim_grid_bound_Ca = Bound_Ca
            ,sim_grid_free_Annexin = Free_Annexins
            ,sim_grid_bound_Annexin = Bound_annexins
            ,hole_pos=[x0,y0]
            ,num_of_rings = 10
            ,inside_val = inside_val
        )

    
    if save_data == True:
        if 'Ring sum list Ca' in df.columns:
            df = df.drop(columns=['Ring sum list Ca'])
        if 'Ring sum list Annexin' in df.columns:
            df = df.drop(columns=['Ring sum list Annexin'])

    

        df2 = pd.DataFrame({
            'Ring sum list Ca': [Ring_sum_list_Ca],
            'Ring sum list Annexin': [Ring_sum_list_Annexin]
        })

    

        df = df.append(df2,ignore_index=True)
        print("\n \n \n")

    
        print(df.info())

        #fig_folder_path =  fig_save_path + f"simtime={Real_sim_time}\\"

        #fig_name = f"Simulation_data_simtime={Real_sim_time}.pkl"
        df.to_pickle(fig_folder_path + df_name)

    print("\n \n ---------- Done ---------- \n \n")





if __name__ == "__main__":
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val,Real_sim_time, real_close_time = constants()
    
    Real_sim_time = 121

    fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
    fig_folder_path =  fig_save_path + f"simtime={Real_sim_time}\\"
    fig_name = f"Simulation_data_simtime={Real_sim_time}.pkl"

    #fig_folder_path =  fig_save_path + f"Cell structure simtime={Real_sim_time}\\"
    #fig_name = f"Cell structure Simulation_data_simtime={Real_sim_time}.pkl"
    main_ring_summing(
        fig_save_path = fig_save_path
        ,fig_folder_path = fig_folder_path
        ,df_name = fig_name
        ,hole_pos=""#[35,5]
    )