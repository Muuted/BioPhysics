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
    
    ref_fig = 4

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
    
    micro_meter = 20e-6
    pixels_for_20mum = micro_meter/dx
    x_start = 1
    x_line = [x_start, x_start + int(pixels_for_20mum)]
    y_pos = 78
    y_width = 4
    y_line = [y_pos,y_pos]

    plt.plot(
    x_line,y_line
    ,color="white",linewidth=y_width
    )

    plt.text(
        x=(x_line[1]-x_line[0])/3
        ,y =y_pos - y_width/2 + 1
        ,s = r"$20 \mu m $"
        ,color="white"
        ,fontsize = 15
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

    micro_meter = 20e-6
    pixels_for_20mum = micro_meter/dx
    x_start = 1
    x_line = [x_start, x_start + int(pixels_for_20mum)]
    y_pos = 78
    y_width = 4
    y_line = [y_pos,y_pos]
    
    #fig.canvas.manager.window.showMaximized()
    
    plt.matshow(Visual_grid)
    plt.title(
        #f"Rings for radial distribution of either Calcium or Annexins \n "
        f"Visual representation of the radial distribution \n "
        +f" Where the opening is placed af (x,y)={x0,y0}"
        #+f", the white star"
        ,fontsize=15
        )
    plt.ylabel(r"y",fontsize=15)    
    plt.xlabel(r"x",fontsize=15)
    plt.colorbar()

    plt.plot(
        x_line,y_line
        ,color="white",linewidth=y_width)
    
    """plt.plot(
    [34, 34],[4,4]
    ,color="white"
    ,linewidth=y_width
    ,marker='*'
    )"""

    plt.text(
        x=(x_line[1]-x_line[0])/3
        ,y =y_pos - y_width/2 + 1
        ,s = r"$20 \mu m $"
        ,color="white"
        ,fontsize = 15
        )
    
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
    

    time_list = [int((i-frame_open)*dt_frame) for i in frame_list]
    t_index_list = [int(i/dt) for i in time_list]
    

    df_sim = pd.read_pickle(fig_folder_path + fig_name_df)
    print(df_sim.info())

    
    Free_Ca = df_sim['Free Calcium'][0]
    Free_annexin = df_sim['Free Annexins'][0]
    Bound_annexin = df_sim['Bound Annexins'][0]
    ref_struct = df_sim['Reference Matrix'][0]
    dx = df_sim['dx'][0]

    
    fig, ax = plt.subplots(1,len(frame_list))
    fig.canvas.manager.window.showMaximized()

    #fig.colorbar()
    #fig.suptitle("hello")

    vmins = 0 #1e-9
    vmaxs = 3e-4
    
    x_label_list =['i)','j)','k)','l)']

    micro_meter = 20e-6
    pixels_for_20mum = micro_meter/dx
    x_start = 5
    x_line = [x_start, x_start + int(pixels_for_20mum)]
    y_pos = 75
    y_width = 4
    y_line = [y_pos,y_pos]


    for i in range(len(frame_list)):
        matlab_figures ="matlab_cell_frame={frame_list[i]}.svg"

        ax[i].plot(x_line,y_line,color="y",linewidth=y_width)

        ax[i].text(
            x=(x_line[1]-x_line[0])/3
            ,y =y_pos - y_width/2
            ,s = r"$20 \mu m $"
            ,color="y"
            ,fontsize = 12
            )
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


def radial_compare_times():
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
    Real_time_steps_data = 235

    """   The simulation data loaded  """
    df_sim = pd.read_pickle(fig_folder_path + fig_name_df)
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

    data_Ca = pd.read_csv(data_path + Ca_data_exp)
    data_Ann= pd.read_csv(data_path + Annexin_data_exp)
    
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

    vec = np.linspace(
        start=0
        ,stop=exp_data_shapeX
        ,num=exp_data_shapeX
        )    

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
    cmap_type = "gist_ncar"
    #cmap_type = "nipy_spectral"
    #cmap_type = "gnuplot"
    
    cmap_type = matplotlib.cm.get_cmap('gist_ncar')
    #cmap_type.set_under('black')
    
    vmin_val_Ca , vmax_val_Ca = 1e-7 , 6e-4
    vmin_val_Ann, vmax_val_Ann = 2e-6 , 8e-6
    normalize = True
    
    end_i = exp_data_shape_t - 1 - frame_open - 1

    min_sim_ann = min([min(sim_ring_data_Ann[t]) for t in range(Tsize-1)])

    frame_list =[ int(i*(234/120)) for i in range(0,16,2)]
    x_label_list =["a)","a)","b)","b)","c)","c)","d)","d)"]
    white_space_str = "                                  "      
    fig.suptitle(
        f"The radial distribution of the concentration of Calcium ions from"
        +f" the opening in the membrane. \n Shown at different times in the interval"
        +f" of approximately {0}s to {round(frame_list[len(frame_list)-1]*(120/234),1)}s"
        ,fontsize=20
        )
    i = 0
    for h in range(2):
        for k in range(2):
            for _ in range(2):
                fig.canvas.manager.window.showMaximized()

                j = frame_list[i] + frame_open
                t_real = frame_list[i]*(120/234)
                t = frame_list[i]
                t_show = round(t_real,1)

                ax[h,k].plot(vec,sim_ring_data_Ca[t]
                                ,label=f"simulation (Molar) t"+r"$ \approx $"+f"{t_show}"
                                ,linestyle="-"
                                ,marker = "o"
                                )
                ax[h,k].plot(vec,real_data_Ca[j]
                                ,label=f"Experiment t"+r"$ \approx $"+f"{t_show}"
                                ,linestyle="--"
                                ,marker = "*"
                                )
                ax[h,k].set_xlabel(x_label_list[i] + white_space_str + f"Ring"
                                    ,fontsize = 15
                                    ,x=0.70
                                    )
                ax[h,k].set_ylabel(r"[Ca]       "
                                    ,rotation='horizontal'
                                    ,fontsize=15,y=0.45
                                    )
                ax[h,k].set_ylim(np.min(sim_ring_data_Ca),max_ca_sim*1.1)
                ax[h,k].legend(
                    fontsize=15
                )

                i+= 1
    
    
    fig, ax = plt.subplots(2,2)
    fig.suptitle(
        f"The radial distribution of the concentration of total annexin proteins from"
        +f" the opening in the membrane. \n Shown at different times in the interval"
        +f" of approximately {0}s to {round(frame_list[len(frame_list)-1]*(120/234),1)}s"
        ,fontsize=20
        )
    i = 0
    for h in range(2):
        for k in range(2):
            for _ in range(2):
                fig.canvas.manager.window.showMaximized()

                j = frame_list[i] + frame_open
                t_real = frame_list[i]*(120/234)
                t = frame_list[i]
                t_show = round(t_real,1)

                ax[h,k].plot(vec,sim_ring_data_Ann[t]
                                ,label=f"simulation (Molar) t"+r"$ \approx $"+f"{t_show}"
                                ,linestyle="-"
                                ,marker = "o"
                                )
                ax[h,k].plot(vec,real_data_Ann[j]
                                ,label=f"Experiment t"+r"$ \approx $"+f"{t_show}"
                                ,linestyle="--"
                                ,marker = "*"
                                )
                ax[h,k].set_xlabel(x_label_list[i] + white_space_str + f"Ring"
                                    ,fontsize = 15
                                    ,x=0.70
                                    )
                ax[h,k].set_ylabel(r"[Ca]       "
                                    ,rotation='horizontal'
                                    ,fontsize=15,y=0.45
                                    )
                #ax[h,k].set_ylim(np.min(sim_ring_data_Ca),max_ca_sim*1.1)
                ax[h,k].set_ylim(min_sim_ann,max_ann_sim*1.1)
                ax[h,k].legend(
                    fontsize=15
                )

                i+= 1    
    

    plt.show()


if __name__ == "__main__":
    
    #plotting_ref_ring_structs()
    plotting_ref_cell_structs()
    #plotting_time_evolution_stabilization()
    #visual_rings()
    #matshow_compare_times()
    #radial_compare_times()