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

def plotting_ref_ring_structs():
    import matplotlib as mpl
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val,Real_sim_time, real_close_time = constants()
    Real_sim_time = 120
    

    data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"


    fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
    fig_folder_path =  fig_save_path + f"Cell structure simtime={Real_sim_time}\\"
    video_save_path = fig_folder_path + f"video_folder\\"     

    fig_name_df = f"Cell structure Simulation_data_simtime={Real_sim_time}.pkl"

    ref_structure = init_ref_circle(
            boxlen=len_size
            ,Radius=R,dRadius=dR ,offsets=[x0,y0]
            ,inside_val=inside_val
            ,outside_val=outside_val
            ,wall_val=wall_val
            )
    
    open_close_membrane(
        Grid=ref_structure
        ,Radius=R, dRadius=dR
        ,offsets=[x0,y0],holesize=holesize
        ,open_val=open_val
        ,wall_val=wall_val
        ,open_wall=True
                        )

    Ysize, Xsize = np.shape(ref_structure)
    for y in range(Ysize):
        for x in range(Xsize):
            if ref_structure[y][x] == outside_val:
                ref_structure[y][x] = 0

            if ref_structure[y][x] == inside_val:
                ref_structure[y][x] = 1

            if ref_structure[y][x] == open_val:
                ref_structure[y][x] = 2
            
            if ref_structure[y][x] == wall_val:
                ref_structure[y][x] = 3
    
    c = mpl.colors.ListedColormap(['black', 'blue', 'red', 'yellow'])
    n = mpl.colors.Normalize(vmin=0,vmax=3) 
    
    mpl.rcParams['xtick.labelsize'] = 12.0
    mpl.rcParams['ytick.labelsize'] = 12.0
    plt.matshow(ref_structure
                ,cmap=c,norm=n
                )
    plt.colorbar()
    plt.title(
        f"Simulation Ring structure \n "
        +f"The cell wall is shown in yellow \n "
        +f" The opening in the membrane is shown in red"
        ,fontsize=15
        )
    
    plt.show()

def plotting_ref_cell_structs():
    import matplotlib as mpl
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val,Real_sim_time, real_close_time = constants()
    Real_sim_time = 120
    
    ref_fig = 27

    if ref_fig == 4:
        hole_pos1 = [34,4]
    if ref_fig == 27:
        hole_pos1 = [37,4]

    data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"
    ref_struct_name_cell = f"ref_struct_from_Ca_filenum{ref_fig}.txt"

    fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
    fig_folder_path =  fig_save_path + f"Cell structure simtime={Real_sim_time}\\"
    video_save_path = fig_folder_path + f"video_folder\\"     

    fig_name_df = f"Cell structure Simulation_data_simtime={Real_sim_time}.pkl"

    ref_structure,outside_val,inside_val,wall_val,hole_pos = make_ref_structure(
        path=data_path
        ,ref_name=ref_struct_name_cell
        ,hole_pos=hole_pos1
    )
    open_close_membrane2(
        Grid=ref_structure
        ,Xholesize=1
        ,Yholesize=4
        ,open_val=open_val
        ,wall_val=wall_val
        ,open_wall_bool= True
        ,offsets=hole_pos
        )

    Ysize, Xsize = np.shape(ref_structure)
    for y in range(Ysize):
        for x in range(Xsize):
            if ref_structure[y][x] == outside_val:
                ref_structure[y][x] = 0

            if ref_structure[y][x] == inside_val:
                ref_structure[y][x] = 1

            if ref_structure[y][x] == open_val:
                ref_structure[y][x] = 2
            
            if ref_structure[y][x] == wall_val:
                ref_structure[y][x] = 3
    
    c = mpl.colors.ListedColormap(['black', 'blue', 'red', 'yellow'])
    n = mpl.colors.Normalize(vmin=0,vmax=3) 
    
    mpl.rcParams['xtick.labelsize'] = 12.0
    mpl.rcParams['ytick.labelsize'] = 12.0
    plt.matshow(ref_structure
                ,cmap=c,norm=n
                )
    plt.colorbar()
    plt.title(
        f"Simulation structure of the cell \n "
        +f"The cell wall is shown in yellow \n "
        +f" The opening in the membrane is shown in red"
        ,fontsize=15
        )
    
    plt.show()

def plotting_time_evolution():
    
    Real_sim_time = 120
    
    data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"
    ref_struct_name_cell = "ref_struct_from_Ca_filenum4.txt"

    fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
    fig_folder_path =  fig_save_path + f"Cell structure simtime={Real_sim_time}\\"
    video_save_path = fig_folder_path + f"video_folder\\"     

    fig_name_df = f"Cell structure Simulation_data_simtime={Real_sim_time}.pkl"
    pass



if __name__ == "__main__":

    #plotting_ref_cell_structs()
    plotting_ref_ring_structs()