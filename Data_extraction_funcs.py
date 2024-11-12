import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def make_ref_structure(path,ref_name):
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
    wall_vall = val_list[2]

    py_ref_struct = np.zeros(shape=(y0,x0))
    for y in range(y0):
        for x in range(x0):
            py_ref_struct[y][x] =ref_structure.loc[y][x]

    #plt.matshow(py_ref_struct)
    #plt.show()
    #xoffset = input("Choose the x value for the hole")
    xoffset = 70
    for y in range(y0):
        a = py_ref_struct[y][int(xoffset)]
        if a == wall_vall:
            yoffset = y
            break
    return py_ref_struct,outside_val,inside_val,wall_vall,xoffset,yoffset


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


def Ring_sum(ref_Grid ,offsets:tuple 
             ,Radius:int, dRadius:int
             ,num_of_rings:int
             )-> tuple:
    x0,y0 = offsets[0],offsets[1]
    Grid = np.copy(ref_Grid)
    grid_shape = np.shape(ref_Grid)
    Visual_Grid = np.zeros(shape=(grid_shape[0],grid_shape[1]))
    xstart = x0 + Radius + 2*dRadius
    xend = x0 - Radius - 2*dRadius
    
    Ring_radius = []
    Ring_sums = []
    
    boxlen = np.shape(Grid)[0]
    R_steps = int(2*Radius/num_of_rings) # gives 10 ringes for summing.

    # first we want to prepare the matrix, so that
    # the out side values arent summed.
    for x in range(boxlen):
        for y in range(boxlen):
            r = np.sqrt( (x-x0)**2 + (y-y0)**2)
            if r > Radius:
                Grid[y][x] = 0

    visual_sum = 0
    for R in range(0,int(2*Radius + x0 ),R_steps):
        sum = 0
        visual_sum += 100
        Z = 1 # normalizing constant
        for x in range(xstart,xend,-1):#(xstart,xend,-1):
            for y in range(boxlen):
                if x < boxlen-1 and y < boxlen-1:
                    r1 = np.sqrt((x-xstart)**2 + (y-y0)**2)
                    r2 = np.sqrt((x-x0)**2 + (y-y0)**2)
                    if R - R_steps < r1 <= R and r2 < Radius + dRadius/2:
                        sum += Grid[y][x] 
                        Z += 1
                        Visual_Grid[y][x] = visual_sum
        
        Ring_sums.append(sum/Z)
        Ring_radius.append(R)

    return Ring_sums, Ring_radius, Grid,Visual_Grid


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