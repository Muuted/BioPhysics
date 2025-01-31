import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#from Constants import constants

def make_ref_structure(
        path
        ,ref_name
        ,hole_pos
        ):
    

    ref_string = "\\" + "\\"
    for i in range(3):
        if path[len(path)-2:len(path)-1] != ref_string:
            if ref_name[0:2] != path[len(path)-2]:
                ref_name = "\\" + ref_name #just making sure the name is do able
    #ref_struct_name = "ref_struct__filenum4.txt"

    ref_structure = pd.read_csv(path + ref_name,delimiter=',',header=None)

    y0,x0 = ref_structure.shape
    val_list = []

    for y in range(y0):
        for x in range(x0):
            a = ref_structure.loc[y][x]
            if a not in val_list:
                val_list.append(a)

    val_list.sort()

    if len(val_list) > 3:
        print(f"val list to finde inside,outside and wall vals are to big \n its len(val_list)={len(val_list)}")
        exit()

    outside_val = val_list[0]
    inside_val = val_list[1]
    wall_val = val_list[2]

    py_ref_struct = np.zeros(shape=(y0,x0))

    for y in range(y0):
        for x in range(x0):
            py_ref_struct[y][x] = ref_structure.loc[y][x]

   
    
    min_x_wall_pos, max_x_wall_pos = 2*x0,0
    min_y_wall_pos, max_y_wall_pos = 2*y0,0
    for y in range(y0):
        for x in range(x0):
            if py_ref_struct[y][x] == wall_val:

                if min_x_wall_pos > x:
                    min_x_wall_pos = x
                if min_y_wall_pos > y :
                    min_y_wall_pos = y
                
                if max_x_wall_pos < x :
                    max_x_wall_pos = x
                if max_y_wall_pos < y :
                    max_y_wall_pos = y
    
    final_py_ref_struct = py_ref_struct[min_y_wall_pos-3:max_y_wall_pos+3,min_x_wall_pos-3:max_x_wall_pos+3]

    """for y in range(np.shape(final_py_ref_struct)[0]):
        for x in range(np.shape(final_py_ref_struct)[1]):
            if x == 0 or x == np.shape(final_py_ref_struct)[1]-1:
                final_py_ref_struct[y][x] = wall_val
            if y == 0 or y == np.shape(final_py_ref_struct)[0]-1:
                final_py_ref_struct[y][x] = wall_val
    """


    print(f"shape of final structure = {np.shape(final_py_ref_struct)}")
    
    if type(hole_pos) == str:
        hole_pos = hole_pos.strip()
    
    if hole_pos == "" :
        plt.matshow(final_py_ref_struct)
        plt.show(block=False)
        hole_pos_x = int(input("Choose the x value for the hole = "))
        hole_pos_y = int(input("Choose the y value for the hole = "))
        hole_pos = [hole_pos_x,hole_pos_y]

    return final_py_ref_struct,outside_val,inside_val,wall_val,hole_pos


def sum_annexin(A_free,A_bound):
    sumAfree = []
    sumAbound = []
    sumtot = []
    
    shape = np.shape(A_free) # shoulde be identical
    sumfree = 0
    sumbound = 0
    for t in range(shape[0]):

        sumfree = np.sum(A_free[t])
        sumbound = np.sum(A_bound[t])
        
        sumAfree.append(sumfree)
        sumAbound.append(sumbound)
        sumtot.append(sumfree + sumbound)

    return sumAfree, sumAbound, sumtot


def Ring_sum(
        ref_grid
        ,sim_grid_free_Ca,sim_grid_bound_Ca
        ,sim_grid_free_Annexin,sim_grid_bound_Annexin
        ,hole_pos:tuple
        #,num_of_rings:int
        ,inside_val:int
        ,return_visplot = False
            )-> tuple:

    num_of_rings = 30
    
    Tsize, Ysize, Xsize = np.shape(sim_grid_free_Ca)
    x0,y0 = hole_pos 

    Radius_vec = np.linspace(
        start=0
        ,stop=Ysize
        ,num=num_of_rings
        )

    Ring_sum_list_Ca = np.zeros(shape=(Tsize,len(Radius_vec)))
    Ring_sum_list_Annexins = np.zeros(shape=(Tsize,len(Radius_vec)))

    Visual_grid = np.zeros(shape=(Ysize,Xsize))
    
    Vis_val = 10
    dVis_val = 10

    for t in range(Tsize):
        if t%(int(Tsize/10)) == 0:
            #print(f"time={t} of {T_tot}")
            print(f" Ring summing Progress :  {int((t/Tsize)*100)} %")   
        for R in range(0,len(Radius_vec)-1):

            R1 = int(Radius_vec[R+1])
            R2 = int(Radius_vec[R])
            sum_Ca = 0
            sum_Annexin = 0
            norm = 0
            for y in range(Ysize):
                for x in range(Xsize):
                    background_Ca = sim_grid_free_Ca[0][y][x] + sim_grid_bound_Ca[0][y][x]
                    background_Ann = sim_grid_free_Annexin[0][y][x] + sim_grid_bound_Annexin[0][y][x]
                    if 0 <= x < Xsize and 0 <= y < Ysize:
                        if ref_grid[y][x] == inside_val:
                            r = np.sqrt( (x - x0)**2 + (y - y0)**2 )

                            if R2 <= r <= R1 :
                                norm += 1
                                sum_Ca += sim_grid_free_Ca[t][y][x] + sim_grid_bound_Ca[t][y][x] - background_Ca
                                sum_Annexin += sim_grid_free_Annexin[t][y][x] + sim_grid_bound_Annexin[t][y][x] - background_Ann
                                if t == 0:
                                    Visual_grid[y][x] = Vis_val

            if norm == 0:
                norm = 1
            Vis_val += dVis_val
            Ring_sum_list_Ca[t][R] = sum_Ca/norm
            Ring_sum_list_Annexins[t][R] = sum_Annexin/norm

        #if t == 0:
         #   plt.matshow(Visual_grid)
          #  plt.title("Visual Grid, does this look correct? \n if yes, just close figure and sum will continue")
            #plt.show()

        
        
    if return_visplot == True:
        return Ring_sum_list_Ca, Ring_sum_list_Annexins, Visual_grid
    else:
        return Ring_sum_list_Ca, Ring_sum_list_Annexins 


def sum_in_cell(ref_Grid,Matrix_Free,Matrix_Bound,inside_val:int)->list:
    print(f"shape of Matrix = {np.shape(Matrix_Free)}")
    T_tot,Ysize , Xsize  = np.shape(Matrix_Free)
    sidelen = np.shape(Matrix_Free)[1]
    conc_time = []

    #Ysize , Xsize = np.shape(ref_Grid)
    for t in range(T_tot):
        totsum = 0
        for y in range(Ysize):
            for x in range(Xsize):
                if ref_Grid[y][x] == inside_val:
                    totsum += Matrix_Free[t][y][x]
                    totsum += Matrix_Bound[t][y][x]

        conc_time.append(totsum)
        
    return conc_time

def remove_conc_outside(ref_grid,grid,outside_val):

    Y_size,X_size = np.shape(ref_grid)
    removed_conc_grid = np.zeros(shape=(Y_size,X_size))

    for y in range(Y_size):
        for x in range(X_size):
            if ref_grid[y][x] == outside_val:
                removed_conc_grid[y][x] = 0
            else:
                removed_conc_grid[y][x] = grid[y][x]
    return removed_conc_grid    


def Area_under_graph(
        graph_data
        ,ring_dist = 1
        ):

    Tsize,rings = np.shape(graph_data)

    area_time = np.zeros(shape=Tsize)
    
    for t in range(Tsize):
        Area = 0
        for i in range(rings-1):
            Area += (ring_dist/2)*(graph_data[t][i] + graph_data[t][i+1])

        area_time[t] = Area

    return area_time


def scale_the_data(
        sim_ring_data_Ca
        ,sim_ring_data_Ann
        ,real_data_Ca
        ,real_data_Ann
                ):
 
    exp_data_shape_t,exp_data_shapeX = np.shape(real_data_Ca)

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

    
def Determining_open_n_close_time(
        data_path
        ,Ca_data_exp
        ,Annexin_data_exp
        ,showplot = False
        ):
    
  

    """  --------- Experimental data loaded  -------------- """

    data_Ca = pd.read_csv(data_path + Ca_data_exp)
    data_Ann= pd.read_csv(data_path + Annexin_data_exp)
    
    exp_data_shape_t, exp_data_shapeX = np.shape(data_Ca)        

    real_data_Ca = np.zeros(shape=(exp_data_shape_t,exp_data_shapeX))
    real_data_Ann = np.zeros(shape=(exp_data_shape_t,exp_data_shapeX))

    for t in range(exp_data_shape_t):
        for R in range(exp_data_shapeX):
            real_data_Ca[t][R] = data_Ca.loc[t][R]
            real_data_Ann[t][R] = data_Ann.loc[t][R]


    tot_ca = [np.sum(real_data_Ca[t][:]) for t in range(exp_data_shape_t)]

    number_of_frames = len(tot_ca)

    avg_start = 0
    N = 5
    for i in range(N):
        avg_start += abs((tot_ca[i+1] - tot_ca[i]))/N

    max_conc =0
    a_max = 0
    frame_open = -1
    frame_close = 0
    for i in range(len(tot_ca)-1):
        a = tot_ca[i+1] - tot_ca[i]
        if a > 50*avg_start and frame_open == -1:
            a_max = a
            frame_open = i
        
        if max_conc < tot_ca[i]:
            max_conc = tot_ca[i]
            frame_close = i
    
    print(f"open frame = {frame_open} and the close frame = {frame_close}")
    print(f"total frame = {len(tot_ca)}")
    if showplot == True:
        tot_ca = [i/max(tot_ca) for i in tot_ca]
        time_vec = np.linspace(
            start=0
            ,stop=120
            ,num=len(tot_ca)
            )
        plt.figure()
        plt.plot(time_vec,tot_ca,marker='*',color="red",label=r" $ [Ca_{tot}] $")
        plt.vlines(
            x=time_vec[frame_open],ymin=min(tot_ca),ymax=max(tot_ca)
            ,label=f"hole opens " +r"t$ \approx $" + f"{int(time_vec[frame_open])}s"
            )
        plt.vlines(
            x=time_vec[frame_close],ymin=min(tot_ca),ymax=max(tot_ca),linestyles='--'
            ,label="hole closes " + r"t$ \approx $" + f"{int(time_vec[frame_close])}s"
            )
        
        plt.title(
            f"Finding the opening time and closing time of the hole in the cell \n"
            #+r"open time$ \approx $"
            #+ f"{int(time_vec[frame_open])}s  and close time" 
            #+r"$ \approx $"
            #+f"{int(time_vec[frame_close])}s \n"
            +r"closure time  $ t_{closure} \approx $"
            +f"{int(time_vec[frame_close] - time_vec[frame_open])}s"
            ,fontsize=15)
        plt.ylabel(
            r" $ \frac{[Ca_{tot}]}{max([Ca_{tot}])} $ " 
            + "         "
            ,fontsize=20
            ,rotation=0
                   )
        plt.xlabel( f" Time (seconds)",fontsize=15)
        
        plt.legend(fontsize=15)
        plt.show()

    return frame_open, frame_close, number_of_frames
    
if __name__ == "__main__":
    """
    const_list = constants()
    c_in ,c_out, D_Ca_cyto, T_tot, len_size, dx, dy, k1, k2 = const_list[0:9]
    c_in_annexin ,bound_annexin_start ,D_Annexin_cyto = const_list[9:12]
    dt ,close_time, c_pump, holesize ,dR ,R ,x0 ,y0 = const_list[12:20]
    wall_val ,inside_val ,outside_val ,open_val = const_list[20:24]
    Real_sim_time, real_close_time = const_list[24:26]
    ref_struct_name_cell ,fig_save_path = const_list[26:28]
    fig_folder_path ,video_save_path ,fig_name_df, data_path = const_list[28:32]
    Ca_data_exp ,Annexin_data_exp = const_list[32:34]

    Determining_open_n_close_time(
        data_path=data_path
        ,Ca_data_exp=Ca_data_exp
        ,Annexin_data_exp=Annexin_data_exp
    )
    """