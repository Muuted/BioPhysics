from Circle_funcs import stabil_condi
import matplotlib.pyplot as plt
from Circle_funcs import Annexin_stablilization
import numpy as np
from Data_extraction_funcs import Determining_open_n_close_time

def constants(
          print_vals = True
            ):
    
    Avogadro = 6.02214076e23 # 1/mol
    # time and step size, and diffusion constant
    picture_size = 83e-6 # meters
    Real_sim_time =  110 #seconds
    real_close_time =  10 #s
    
    ref_fig_num_list = [4,27]
    ref_fig_num = 4
    ref_struct_name_cell = f"ref_struct_from_Ca_filenum{ref_fig_num}.txt"
    data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"
    fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
    fig_folder_path =  fig_save_path + f"Cell structure {ref_fig_num} simtime={Real_sim_time}\\"
    video_save_path = fig_folder_path + f"video_folder\\"     
    fig_name_df = f"Cell structure {ref_fig_num} Simulation_data_simtime={Real_sim_time}.pkl"
    Ca_data_experiment = f"Ring_sum_data_Ca_filenum{ref_fig_num}.txt"
    Annexin_data_experiment = f"Ring_sum_data_Annexins_filenum{ref_fig_num}.txt"
    
    
    frame_open, frame_close, data_number_of_frames = Determining_open_n_close_time(
        data_path=data_path
        ,Ca_data_exp=Ca_data_experiment
        ,Annexin_data_exp=Annexin_data_experiment
        ,showplot=False
    )
    

    data_real_time = 120
    data_time_pr_frame = data_real_time/data_number_of_frames
    
    data_hole_open_time = frame_open*data_time_pr_frame
    real_close_time = (frame_close - frame_open)*data_time_pr_frame
    
  
    len_size = 80 # number of grid points
    dx, dy = picture_size/len_size ,picture_size/len_size # meters

    # Concentrations for Ca^2+ inside; outside cell
    c_out =2e-3# 2*10**(-3)#e-3 #M Concentration outside of cell
    c_in = 100e-9#*10*(-9)#e-9 #M Concetration inside cell

    D_Ca_cyto  = 2.7e-11#*10**(-11)#e-11 #meters^2/second
    D_Annexin_cyto  = 2.7e-11 #5.0e-11 #meters^2/second     # 4.4e-13  or 4.4e-12 maybe?

    dt = stabil_condi(dt=0.1,dx=dx,dy=dy
                    ,D_list=[ D_Ca_cyto , D_Annexin_cyto ]
                    )
    
    
    T_tot = int(Real_sim_time/dt)

    pump_multiplyer = 30 # found from sim testing.
    c_pump = dt*(5e-8)*pump_multiplyer # Molar/second
    c_pump_base = (3.028415305204001e-09)/(5.104687e-03)
    c_pump =  c_pump_base*dt
    c_pump = float('{:0.2e}'.format(c_pump))
    

    # Anxexin constants.  
    A_total_density = (2.0e6)/(1.0e-18) #number/meter^3 from article directly
    A_total_conc = (A_total_density/Avogadro)*(0.1/100.0)/1000.0 #should get 0.1% of the 3D volume concentration for all proteins

    # association and dissociation constants
    k1 = 1e4 # 1/(Ms)
    k2 = 1e-1 #1/s 

    
    bound_annexin_start = k1*A_total_conc*c_in/(k1*c_in + k2)
    c_in_annexin = k2*A_total_conc/(k1*c_in+k2)
    #c_in_annexin = A_total_conc - bound_annexin_start 


    #Closing mechanism time
    close_time = int(real_close_time/dt)

    #size of cell, hole in cell and the center's placement.
    holesize = 1
    dR = int(2)
    R = int(len_size*0.4)# - 3*dR) 
    
    x0,y0 = int(len_size/2), int(len_size/2)

    if ref_fig_num == 4:
            x0,y0 = [34,4]
    if ref_fig_num == 27:
        x0,y0 = [37,4]

    # Defining areas of interest, values for regonition.
    wall_val = 100
    open_val = 20
    inside_val = 10
    outside_val = -5
    
    if print_vals == True:
        print(
            f" \n \n"
            + "------------- Constant used in Simulation -------------- \n "
            + f"    Real simulation time = {Real_sim_time} s \n "
            + f"    hole close time = {int(real_close_time)} s \n " 
            + f"    dx=dy = {dx:e} m \n "
            + f"    dt = {dt:e} s \n "
            + f"    [Ca] outside cell = {c_out:e} \n "
            + f"    [Ca] inside cell = {c_in:e} \n "
            + f"    Diffuse constant Ca in cell = {D_Ca_cyto} m^2/s \n "
            + f"    Diffuse constant Annexin in cell = {D_Annexin_cyto} m^2/s \n "
            + f"    c_pump = {c_pump} M/s \n "
            + f"    Total initial conc Annexins = {A_total_conc:e} M \n "
            + f"    c_in_ann = {c_in_annexin:e} M \n "
            + f"    and the bound_ann_start = {bound_annexin_start:e} M \n "
            + f"    k1 = {k1:0.1e} 1/Ms \n "
            + f"    and k2 = {k2:0.1e} 1/s \n "
            + f"    total number of sim steps = {T_tot:e} steps \n "
            #+ f"    sim equil time = {1/(k1*c_in_annexin + k2)} \n "
            #+ f"    size of grid = ({len_size},{len_size}) \n "
            f" ------------------------------------------------------ \n \n "
        )

    args_list = [
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2
        ,c_in_annexin
        ,bound_annexin_start
        ,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,Real_sim_time, real_close_time
        ,ref_struct_name_cell ,fig_save_path
        ,fig_folder_path ,video_save_path 
        ,fig_name_df ,data_path
        ,Ca_data_experiment ,Annexin_data_experiment
        ,frame_open ,frame_close
    ]
    
    return args_list






if __name__ == "__main__":
    A = constants()