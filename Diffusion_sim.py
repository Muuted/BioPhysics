import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from diff_sim_funcs import dCondt, check_volume


# Concentrations for Ca^2+ inside; outside cell
c_out = 2e-3 #M Concentration outside of cell
c_in = 100e-9 #M Concetration inside cell

# time and step size, and diffusion constant
dx, dy = 1,1
dt = 0.01
D  = 0.1
holesize = 2

# Creation of our Grid
T_tot = 1000
Ysize = 20
Xsize = 2*Ysize + 1
Grid_diff = np.zeros(shape=(2,Ysize,Xsize)) # replication of figure 1.b

# Initilizing our Grid. leaving a barrier in the middle,
# that seperates the outside(left) and inside(right) of the cell
for x in range(1,Xsize-1):
    for y in range(1,Ysize-1):
        if x < Ysize:
            Grid_diff[0][y][x] = c_out
        if x > Ysize:
            Grid_diff[0][y][x] = c_in


tot_conc1 = check_volume(grid=Grid_diff[0])

for t in range(0,T_tot-1):

    t1, t2 = t%2, (t+1)%2 # to change the matricies back and forth

    for x in range(1,Xsize - 1):
        for y in range(1,Ysize - 1):

            if x != Ysize:
                dCdt = dCondt(C=Grid_diff[t1]
                              ,x=x,y=y,dx=dx,dy=dy, h=1
                              )
                Grid_diff[t2][y][x]= Grid_diff[t1][y][x] + D*dt*dCdt
            
            if x == Ysize:
                if(Ysize - holesize)/2 < y < (Ysize + holesize)/2:
                    dCdt = dCondt(C=Grid_diff[t1]
                                  ,x=x,y=y,dx=dx,dy=dy,h=1
                                  )
                    Grid_diff[t2][y][x] = Grid_diff[t1][y][x] + D*dt*dCdt

tot_conc2 = check_volume(grid=Grid_diff[0])

print(f"conc2/conc1 ={tot_conc2/tot_conc1}")
plt.matshow(Grid_diff[0])
plt.show()
