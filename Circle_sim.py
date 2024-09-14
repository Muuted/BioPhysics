import matplotlib.pyplot as plt
import numpy as np
from Circle_funcs import init_circle,cicle_boundary,circle_dCondt
from Circle_funcs import replenish

# Concentrations for Ca^2+ inside; outside cell
c_out = 2e-3 #M Concentration outside of cell
c_in = 100e-9 #M Concetration inside cell
conc_list = []

# mechanisms change.
open_condi = True

# time and step size, and diffusion constant
T_tot = 100
len_size = 20
dx, dy = 1,1
dt = 0.1
D  = 0.5
c_pump = 0

# Anxexin constants.
k1,k2 = 1,1

#size of cell, hole in cell and the center's placement.
holesize = 3
R = len_size*0.7/2
dR = 2
center = int(len_size/2)
wall_val = 100

# Creation of our Grids
Free_annexin = np.zeros(shape=(T_tot,len_size,len_size))
Bound_annexin = np.zeros(shape=(T_tot,len_size,len_size))

Grid,ref_Grid = init_circle(boxlen=len_size,time=T_tot
                            ,Radius=R,dRadius=2
                            ,c_in=c_in,c_out=c_out
                            ,refval=100,center=center
                            )

Grid[0][center][center] = 100   

# Initilizing our Grid. leaving a barrier in the middle,
# that seperates the outside(left) and inside(right) of the cell

for t in np.arange(0,T_tot-1):    
    for x in range(0,len_size):
        for y in range(0,len_size):

            if ref_Grid[y][x] != wall_val:
                pos = cicle_boundary(x=x,y=y,boxlen=len_size
                                        ,ref_matrix=ref_Grid
                                        ,refval=wall_val)

                dcdt = circle_dCondt(C=Grid[t],pos=pos,dx=dx,dy=dy)

                Grid[t+1][y][x] = Grid[t][y][x] + dt*D*dcdt

                radii = np.sqrt(    (x-center)**2 + (y-center)**2 )
                if radii < R and Grid[t+1][y][x] > c_pump:
                    Grid[t+1][y][x] += -c_pump # the pumping mechanism
                                               #, for only inside the cell

            
    if t%(int(T_tot*0.1)) == 0:
        conc_list.append(
            np.sum(Grid[t])
                )


plt.matshow(ref_Grid)
plt.matshow(Grid[T_tot-1])
plt.matshow(Grid[0])

plt.figure()
plt.plot(conc_list)
plt.show()
