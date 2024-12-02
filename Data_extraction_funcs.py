import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def make_ref_structure(
        path,ref_name
        ,hole_pos
        ):
    

    ref_string = "\\" + "\\"
    for i in range(3):
        if path[len(path)-2:len(path)-1] != ref_string:
            if ref_name[0:2] != path[len(path)-2]:
                ref_name = "\\" + ref_name #just making sure the name is do able
    #ref_struct_name = "ref_struct__filenum4.txt"

    ref_structure = pd.read_csv(path + ref_name,delimiter=',',header=None)

    y0,x0= ref_structure.shape
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
        ,num_of_rings:int  
        ,inside_val:int
            )-> tuple:


    Tsize, Ysize, Xsize = np.shape(sim_grid_free_Ca)
    x0,y0 = hole_pos #position of the hole in the membrane

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

                    if 0 <= x < Xsize and 0 <= y < Ysize:
                        if ref_grid[y][x] == inside_val:
                            r = np.sqrt( (x - x0)**2 + (y - y0)**2 )

                            if R2 <= r <= R1 :
                                norm += 1
                                sum_Ca += sim_grid_free_Ca[t][y][x] + sim_grid_bound_Ca[t][y][x]
                                sum_Annexin += sim_grid_free_Annexin[t][y][x] + sim_grid_bound_Annexin[t][y][x]
                                if t == 0:
                                    Visual_grid[y][x] = Vis_val

            if norm == 0:
                norm = 1
            Vis_val += dVis_val
            Ring_sum_list_Ca[t][R] = sum_Ca/norm
            Ring_sum_list_Annexins[t][R] = sum_Annexin/norm

        if t == 0:
            plt.matshow(Visual_grid)
            plt.title("Visual Grid, does this look correct? \n if yes, just close figure and sum will continue")
            plt.show(block=False)
            
        
        

    return Ring_sum_list_Ca, Ring_sum_list_Annexins 




def sum_in_cell(ref_Grid,Matrix_Free,Matrix_Bound,inside_val:int)->list:
    print(f"shape of Matrix = {np.shape(Matrix_Free)}")
    T_tot  = np.shape(Matrix_Free)[0]
    sidelen = np.shape(Matrix_Free)[1]
    conc_time = []

    for t in range(T_tot):
        totsum = 0
        for y in range(sidelen):
            for x in range(sidelen):
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

if __name__ == "__main__":
    pass