import numpy as np
import matplotlib.pyplot as plt


def sum_annexin(A_free,A_bound):
    sumAfree = []
    sumAbound = []
    sumtot = []
    
    shape = np.shape(A_free) # shoulde be identical

    for t in range(shape[0]):
        sumfree = 0
        sumbound = 0
        totsum = 0
        #for x in range(shape[1]):
         #   for y in range(shape[2]):
        sumfree = np.sum(A_free[t])
        sumbound = np.sum(A_bound[t])
        
        totsum += sumfree + sumbound
        
        sumAfree.append(sumfree)
        sumAbound.append(sumbound)
        sumtot.append(totsum)

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


def sum_in_cell(ref_Grid,Matrix,inside_val:int)->list:
    print(f"shape of Matrix = {np.shape(Matrix)}")
    T_tot  = np.shape(Matrix)[0]
    sidelen = np.shape(Matrix)[1]
    conc_time = []

    for t in range(T_tot):
        totsum = 0
        for y in range(sidelen):
            for x in range(sidelen):
                if ref_Grid[y][x] == inside_val:
                    totsum += Matrix[t][y][x]

        conc_time.append(totsum)
        
    return conc_time