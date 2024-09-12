import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from diff_sim_funcs import dCdt


# Concentrations for Ca^2+ inside; outside cell
c_out = 5#2e-3 #M Concentration outside of cell
c_in = -1#3e-9 #M Concetration inside cell


# Creation of our Grid
T_tot = 10
Ysize = 5
Xsize = 2*Ysize + 1
Grid_diff = np.zeros(shape=(T_tot,Ysize,Xsize)) # replication of figure 1.b

# Initilizing our Grid. leaving a barrier in the middle,
# that seperates the outside(left) and inside(right) of the cell
for x in range(1,Xsize-1):
    for y in range(1,Ysize-1):
        if x < Ysize:
            Grid_diff[0][y][x] = c_out
        if x > Ysize:
            Grid_diff[0][y][x] = c_in


#print(Grid_diff[0])
N = 100
test_grid = np.zeros(shape=(2,N,N))
t = 0
T= 100
test_grid[0][int(N/2)][int(N/2)] = 100 

D,dx,dy,dt = 1,1,1,1


plt.matshow(test_grid[0])
plt.show()

while t < T:
    if t%10:
        print(f"time ellapsed={t}, of total time={T}")
    
    for x in range(1,N-1):
        for y in range(1,N-1):
            change_val = dt*D*dCdt(C=test_grid[(t+1)%2],D=D,x=x,y=y,dx=dx,dy=dy,dt=dt)
            test_grid[t%2][y][x] = test_grid[(t+1)%2][y][x] + change_val

    t += 1
plt.matshow(test_grid[0])
plt.show()
