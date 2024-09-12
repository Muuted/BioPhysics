import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from diff_sim_funcs import dCondt, check_volume, boundary_conditions



# Concentrations for Ca^2+ inside; outside cell
c_out = 100#2e-3 #M Concentration outside of cell
c_in = 100#100e-9 #M Concetration inside cell
conc_list = []

# time and step size, and diffusion constant
dx, dy = 1,1
dt = 0.1
D  = 0.1
holesize = 2

# Creation of our Grid
T_tot = 100
Ysize = 20
Xsize = 2*Ysize + 1
Grid_diff = np.zeros(shape=(T_tot,Ysize,Xsize)) # replication of figure 1.b

# Initilizing our Grid. leaving a barrier in the middle,
# that seperates the outside(left) and inside(right) of the cell
for x in range(0,Xsize):
    for y in range(0,Ysize):
        if x < Ysize:
            Grid_diff[0][y][x] = c_out
        if x > Ysize:
            Grid_diff[0][y][x] = c_in


plt.matshow(Grid_diff[0])
plt.show()

conc_list.append(
    #check_volume(grid=Grid_diff[0])
    np.sum(Grid_diff[0])
    )

for t in np.arange(0,T_tot-1):
    t1, t2 = t%2, (t+1)%2 # to change the matricies back and forth
    for x in range(0,Xsize):
        for y in range(0,Ysize):

            pos_list = boundary_conditions(x=x,y=y
                                        ,Xsize=Xsize,Ysize=Ysize
                                        ,holesize=holesize
                                        ,open=False
                                                   )
            if x != Ysize + 1:
                dCdt = dCondt(C=Grid_diff[t]
                              ,pos=pos_list,dx=dx,dy=dy, h=1
                              )
                Grid_diff[t+1][y][x]= Grid_diff[t][y][x] + D*dt*dCdt
            
            if x == Ysize + 1:
                if(Ysize - holesize)/2 < y < (Ysize + holesize)/2:
                    dCdt = dCondt(C=Grid_diff[t]
                                  ,pos=pos_list,dx=dx,dy=dy,h=1
                                  )
                    Grid_diff[t+1][y][x] = Grid_diff[t][y][x] + D*dt*dCdt
    if t%100 == 0:
        conc_list.append(
            np.sum(Grid_diff[t])
                   )



plt.figure()
plt.plot(conc_list)
plt.show()
#print(f"conc2/conc1 ={tot_conc2/tot_conc1}")
plt.matshow(Grid_diff[T_tot-1])
plt.show()
