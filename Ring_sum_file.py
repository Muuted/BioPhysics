import matplotlib.pyplot as plt
import numpy as np

from Constants import constants
import pandas as pd
from Data_extraction_funcs import Ring_sum

def main_ring_summing(
        fig_save_path
        ,fig_folder_path
        ,df_name
        ,hole_pos = ""
        ,sum_quick = False
                ):
        
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

    Free_Annexins = df['Free Annexins'][0]
    Bound_annexins = df['Bound Annexins'][0]

    ref_structure = df['Reference Matrix'][0]
    
    if hole_pos == "":
        x0 = df['hole position X'][0]
        y0 = df['hole position Y'][0]
    
    sim_time_steps = df['time steps'][0]    
    
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

        df.to_pickle(fig_folder_path + df_name)

    print("\n \n ---------- Done ---------- \n \n")



if __name__ == "__main__":
    const_list = constants()
    c_in ,c_out, D_Ca_cyto, T_tot, len_size, dx, dy, k1, k2 = const_list[0:9]
    c_in_annexin ,bound_annexin_start ,D_Annexin_cyto = const_list[9:12]
    dt ,close_time, c_pump, holesize ,dR ,R ,x0 ,y0 = const_list[12:20]
    wall_val ,inside_val ,outside_val ,open_val = const_list[20:24]
    Real_sim_time, real_close_time = const_list[24:26]
    ref_struct_name_cell ,fig_save_path = const_list[26:28]
    fig_folder_path ,video_save_path ,fig_name_df, data_path = const_list[28:32]

    main_ring_summing(
        fig_save_path = fig_save_path
        ,fig_folder_path = fig_folder_path
        ,df_name = fig_name_df
        ,hole_pos= ""
        ,sum_quick= True
    )