import matplotlib.pyplot as plt
import matplotlib
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
from Cell_structure_sim import main_cell_structure_sim

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

def plotting_time_evolution_stabilization():
    const_list = constants()
    c_in ,c_out, D_Ca_cyto, T_tot, len_size, dx, dy, k1, k2 = const_list[0:9]
    c_in_annexin ,bound_annexin_start ,D_Annexin_cyto = const_list[9:12]
    dt ,close_time, c_pump, holesize ,dR ,R ,x0 ,y0 = const_list[12:20]
    wall_val ,inside_val ,outside_val ,open_val = const_list[20:24]
    Real_sim_time, real_close_time = const_list[24:26]
    ref_struct_name_cell ,fig_save_path = const_list[26:28]
    fig_folder_path ,video_save_path ,fig_name_df, data_path = const_list[28:32]
    Ca_data_exp ,Annexin_data_exp = const_list[32:34]
    frame_open ,frame_close = const_list[34:36]

    free_ann_equil = c_in_annexin
    bound_ann_equil = bound_annexin_start

    c_in_annexin += bound_annexin_start
    bound_annexin_start = 0

    

    Sim_data_list = main_cell_structure_sim(
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,open_hole=False
        ,ref_bakteria= ref_struct_name_cell
        ,data_path= data_path
                    )
    
    ref_structure,Free_Ca,Free_annexin,Bound_annexin,Bound_Ca= Sim_data_list

    points_inside = points_inside_cell(
        refgrid=ref_structure
        ,inside_val=inside_val
    )
    sumAfree,sumAbound,sumAtot = sum_annexin(
        A_free=Free_annexin
        ,A_bound=Bound_annexin
    )

    sumAfree =[i/points_inside for i in sumAfree]
    sumAbound =[i/points_inside for i in sumAbound]
    A_f, A_b = Annexin_stablilization(
        k1 =k1,k2 = k2
        ,A_tot = c_in_annexin + bound_annexin_start
        ,c_in= c_in
        ,realtime= Real_sim_time
        ,dt = dt
        )
    
    time_vec1 = np.linspace(
        start=0
        ,stop=Real_sim_time
        ,num=len(A_f)
    )
    time_vec2 = np.linspace(
        start=0
        ,stop=Real_sim_time
        ,num=len(sumAfree)
    )
    
    fig, ax = plt.subplots(1,3)
    fig.canvas.manager.window.showMaximized()

    ax[0].plot(time_vec1,A_f,label="analytical solution",marker='*')#,color="blue")
    ax[0].plot(time_vec2,sumAfree,label="simulation",linestyle='-',color="red")
    ax[0].hlines(y=free_ann_equil,xmin=0,xmax=Real_sim_time,label="Free annexin equilibrium",color="black",linestyle="--")
    ax[0].set_title("Concentration of \n free annexins",fontsize=15)
    ax[0].set_xlabel("Time (s) \n a)",fontsize=15)
    ax[0].tick_params(labelsize="large")
    ax[0].yaxis.offsetText.set_fontsize(15)
    ax[0].legend(fontsize=15)
    

    ax[1].plot(time_vec1,A_b,label="analytical solution",marker='*')#,color="blue")
    ax[1].plot(time_vec2,sumAbound,label="simulation",linestyle='-',color="red")
    ax[1].hlines(y=bound_ann_equil,xmin=0,xmax=Real_sim_time,label="bound annexin equilibrium",color="black",linestyle="--")
    ax[1].set_title("Concentration of \n bound annexins",fontsize=15)
    ax[1].set_xlabel("Time (s) \n b)",fontsize=15)
    ax[1].tick_params(labelsize="large")
    ax[1].yaxis.offsetText.set_fontsize(15)
    ax[1].legend(fontsize=15) 

    A_tot = [A_f[i] + A_b[i] for i in range(len(A_f))]
    sumAtot = [i/points_inside for i in sumAtot]
    offset= 0.0001
    ax[2].plot(time_vec1,A_tot,label="analytical solution",marker='*')
    ax[2].plot(time_vec2,sumAtot,label="simulation",linestyle='-',color="red")
    ax[2].set_title("Total Concentration \n of annexins",fontsize=15)
    ax[2].set_ylim(min(A_tot)*(1-offset),max(A_tot)*(1+offset) )
    ax[2].set_xlabel("Time (s) \n c)",fontsize=15)
    ax[2].tick_params(labelsize="large")
    ax[2].yaxis.offsetText.set_fontsize(15)
    ax[2].legend(fontsize=15) 

    plt.show()


def visual_rings():
    const_list = constants()
    c_in ,c_out, D_Ca_cyto, T_tot, len_size, dx, dy, k1, k2 = const_list[0:9]
    c_in_annexin ,bound_annexin_start ,D_Annexin_cyto = const_list[9:12]
    dt ,close_time, c_pump, holesize ,dR ,R ,x0 ,y0 = const_list[12:20]
    wall_val ,inside_val ,outside_val ,open_val = const_list[20:24]
    Real_sim_time, real_close_time = const_list[24:26]
    ref_struct_name_cell ,fig_save_path = const_list[26:28]
    fig_folder_path ,video_save_path ,fig_name_df, data_path = const_list[28:32]
    Ca_data_exp ,Annexin_data_exp = const_list[32:34]
    frame_open ,frame_close = const_list[34:36]

    
    Sim_data_list = main_cell_structure_sim(
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,open_hole=False
        ,ref_bakteria= ref_struct_name_cell
        ,data_path= data_path
                    )
    
    ref_structure,Free_Ca,Free_annexin,Bound_annexin,Bound_Ca= Sim_data_list

    Ring_sum_list_Ca, Ring_sum_list_Annexins, Visual_grid = Ring_sum(
        ref_grid=ref_structure
        ,sim_grid_free_Ca=Free_Ca
        ,sim_grid_bound_Ca=Bound_Ca
        ,sim_grid_free_Annexin=Free_annexin
        ,sim_grid_bound_Annexin=Bound_annexin
        ,hole_pos=[x0,y0]
        #,num_of_rings=30
        ,inside_val=inside_val
        ,return_visplot=True
    )

    
    #fig.canvas.manager.window.showMaximized()
    plt.matshow(Visual_grid)
    plt.title("")
    plt.show()


def matshow_compare_times():
    const_list = constants(print_vals=False)
    c_in ,c_out, D_Ca_cyto, T_tot, len_size, dx, dy, k1, k2 = const_list[0:9]
    c_in_annexin ,bound_annexin_start ,D_Annexin_cyto = const_list[9:12]
    dt ,close_time, c_pump, holesize ,dR ,R ,x0 ,y0 = const_list[12:20]
    wall_val ,inside_val ,outside_val ,open_val = const_list[20:24]
    Real_sim_time, real_close_time = const_list[24:26]
    ref_struct_name_cell ,fig_save_path = const_list[26:28]
    fig_folder_path ,video_save_path ,fig_name_df, data_path = const_list[28:32]
    Ca_data_exp ,Annexin_data_exp = const_list[32:34]
    frame_open ,frame_close = const_list[34:36]

    matlab_figures_path= "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\figures\\progress comparison\\"
    

    dt_frame = 120/234
    frame_list =[frame_open,25,32,40]
    time_list=[7.984,12.475,15.968,19.96]
    
    matlab_fig_name = [f"matlab_cell_frame={i}.txt" for i in frame_list]

    time_list = [int((i-frame_open)*dt_frame) for i in frame_list]
    t_index_list = [int(i/dt) for i in time_list]
    

    df_sim = pd.read_pickle(fig_folder_path + fig_name_df)
    print(df_sim.info())

    
    Free_Ca = df_sim['Free Calcium'][0]
    Free_annexin = df_sim['Free Annexins'][0]
    Bound_annexin = df_sim['Bound Annexins'][0]
    ref_struct = df_sim['Reference Matrix'][0]

    data_cell_1 = pd.read_csv(matlab_figures_path + matlab_fig_name[0])    
    exp_data_shapeY, exp_data_shapeX = np.shape(data_cell_1)

    #matlab_figures = np.zeros(shape=(len(frame_list),exp_data_shapeY,exp_data_shapeX))

    """for i in range(len(frame_list)):
        print(f" where are at file = {i} now")
        data_cell_1 = pd.read_csv(
            matlab_figures_path + matlab_fig_name[i]
            ,dtype=float
            )    
        exp_data_shapeY, exp_data_shapeX = np.shape(data_cell_1)
        if i == 0:
            matlab_figures = np.zeros(shape=(len(frame_list),exp_data_shapeY,exp_data_shapeX))
        #real_data_cell_1 = np.zeros(shape=(exp_data_shapeY,exp_data_shapeX))
        for y in range(exp_data_shapeY):
            for x in range(exp_data_shapeX):
                matlab_figures[i][y][x] = data_cell_1.loc[y][x]
    """        
        #matlab_figures.append(real_data_cell_1)
    
    
    fig, ax = plt.subplots(1,len(frame_list))
    fig.canvas.manager.window.showMaximized()

    #fig.colorbar()
    #fig.suptitle("hello")

    vmins = 0 #1e-9
    vmaxs = 3e-4
    
    x_label_list =['i)','j)','k)','l)']

    for i in range(len(frame_list)):
        matlab_figures ="matlab_cell_frame={frame_list[i]}.svg"

        k = t_index_list[i]
        t = time_list[i]

        removed_Free_Ca = remove_conc_outside(
            ref_grid=ref_struct
            ,grid=Free_Ca[k]
            ,outside_val=outside_val
        )
        
        pos = ax[i].matshow(
            #Bound_annexin[k]
            removed_Free_Ca
            ,vmin=vmins
            ,vmax=vmaxs
            ,cmap="gray"
                      )
        ax[i].set_title(
            r"simulation at time$ \approx $" + f"{t}s+{int(frame_open*dt_frame)}s"
            ,fontsize=15
            )
        ax[i].set_xlabel(
            x_label_list[i]
            ,fontsize=15
            )
            
        ax[i].set_xticks([])
        ax[i].set_yticks([])
    
    plt.show()



if __name__ == "__main__":
    
    #plotting_ref_cell_structs()
    #plotting_ref_ring_structs()
    #plotting_time_evolution_stabilization()
    #visual_rings()
    matshow_compare_times()