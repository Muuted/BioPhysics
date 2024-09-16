import matplotlib.pyplot as plt
import numpy as np
from Circle_funcs import init_ref_circle    ,cicle_boundary ,circle_dCondt
from Circle_funcs import init_conc  ,open_close_membrane    ,stabil_condi


# Concentrations for Ca^2+ inside; outside cell
c_out = 2e-3 #M Concentration outside of cell
c_in = 100e-9 #M Concetration inside cell
conc_list = []
conc_time_list = []

# time and step size, and diffusion constant
T_tot = 1000
len_size = 40
dx, dy = 0.1,0.1
D  = 0.1
dt = stabil_condi(dt=0.1,dx=dx,dy=dy,D=D)
c_pump = c_in

# mechanisms change.
open_hole = True
close_time = 300

# Anxexin constants.
k1,k2 = 1,1
c_in_annexin = 1

#size of cell, hole in cell and the center's placement.
holesize = 3
R = len_size*0.5/2
dR = 2
x0,y0 = int(len_size/2), int(len_size/2)
wall_val = 100
inside_val = 10
outside_val = -5
open_val = 20

# Creation of our Grids

ref_structure = init_ref_circle(
    boxlen=len_size
    ,Radius=R,dRadius=2 ,offsets=[x0,y0]
    ,inside_val=inside_val
    ,outside_val=outside_val
    ,wall_val=wall_val
    )

Free_Ca = init_conc(
    ref_grid=ref_structure
    ,time=T_tot
    ,c_in=c_in,c_out=c_out
    ,inside_val=inside_val
    ,outside_val=outside_val
                    )
Free_annexin = init_conc(
    ref_grid=ref_structure
    ,time=T_tot
    ,c_in= c_in_annexin
    ,c_out=0
    ,inside_val=inside_val
    ,outside_val=outside_val
                        )

Bound_annexin = np.zeros(shape=(T_tot,len_size,len_size))
Bound_ca = np.zeros(shape=(T_tot,len_size,len_size))
#Free_Ca[0][y0][x0] = 100   


ref_structure = open_close_membrane(
    Grid=ref_structure
    ,Radius=R, dRadius=dR
    ,offsets=[x0,y0],holesize=holesize
    ,open_val=open_val
    ,wall_val=wall_val
    ,open_wall=True
)


for t in np.arange(0,T_tot-1): 
    if t%(T_tot/10) == 0:
        print(f"time={t} of {T_tot}")   
    t1, t2 = t%2, (t+1)%2
    t1, t2 = t, t+1
    for x in range(0,len_size):
        for y in range(0,len_size):
            if ref_structure[y][x] != wall_val:
                pos = cicle_boundary(x=x,y=y,boxlen=len_size
                                        ,ref_matrix=ref_structure
                                        ,refval=wall_val
                                    )

                dcdt = circle_dCondt(C=Free_Ca[t],pos=pos,dx=dx,dy=dy)

                Free_Ca[t+1][y][x] = Free_Ca[t][y][x] + dt*D*dcdt

                radii = np.sqrt( (x-x0)**2 + (y-y0)**2 )

                if radii < R and Free_Ca[t+1][y][x] > c_pump:
                    Free_Ca[t+1][y][x] += -c_pump # the pumping mechanism
                                               #, for only inside the cell

            if x < 3 or y < 3 or x > len_size-4 or y > len_size-4:
                Free_Ca[t+1][y][x] = c_out 
                #We assume infinite supply, i.e the value
                #is allways the same the walls

            if t >= close_time and open_hole==True:
                ref_structure = open_close_membrane(
                    Grid=ref_structure
                    ,Radius=R, dRadius=dR
                    ,offsets=[x0,y0],holesize=holesize
                    ,open_val=open_val
                    ,wall_val=wall_val
                    ,open_wall=False
                )
                open_hole = False
                print(f"wall closure time t={t}")
            
    #if t%(int(T_tot*0.1)) == 0:
    conc_list.append(np.sum(Free_Ca[t]))
    conc_time_list.append(t)



#plt.figure()
#plt.plot(Free_Ca[90][y0][:],label="t=90,y=y0")
#plt.legend()


plt.matshow(Free_Ca[T_tot-1])
plt.title("free ca end")


plt.figure()
plt.plot(conc_time_list,conc_list)
plt.xlabel("timestep")
plt.title("concentration over time")


plt.show()