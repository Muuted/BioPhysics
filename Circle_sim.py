import matplotlib.pyplot as plt
import numpy as np
import random
from Circle_funcs import init_ref_circle    ,cicle_boundary ,circle_dCondt
from Circle_funcs import init_conc  ,open_close_membrane    ,stabil_condi
from Circle_funcs import Ring_sum, circle_dAfreedt, circle_dAbounddt,sum_annexin

# Concentrations for Ca^2+ inside; outside cell
c_out = 2e-3 #M Concentration outside of cell
c_in = 100e-9 #M Concetration inside cell
conc_list = []
conc_time_list = []

# time and step size, and diffusion constant
T_tot = 100
len_size = 40
dx, dy = 0.1 ,0.1
D  = 0.1
dt = stabil_condi(dt=0.1,dx=dx,dy=dy,D=D)
c_pump = c_in

# mechanisms change.
open_hole = True
close_time = T_tot*0.9

# Anxexin constants.
k1 ,k2 = 0.1 ,0.1
c_in_annexin = c_in
prob_free_ann = 0
prob_bound_ann = 0

#size of cell, hole in cell and the center's placement.
holesize = 3
R = int(len_size*0.5/2)
dR = int(2)
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

                Free_Ca[t+1][y][x] = Free_Ca[t][y][x] + dt*dcdt

                radii = np.sqrt( (x-x0)**2 + (y-y0)**2 )

                if radii < R :
                    if prob_free_ann < random.random():
                        dAfreedt = circle_dAfreedt(
                            A_free= Free_annexin[t]
                            ,A_bound=Bound_annexin[t]
                            ,C=Free_Ca[t]
                            ,pos=pos ,dx=dx ,dy=dy ,k1=k1 ,k2=k2,D=D
                        )
                        Free_annexin[t+1][y][x] = Free_annexin[t][y][x] + dt*dAfreedt
                    
                    if prob_free_ann < random.random():
                        dAbounddt = circle_dAbounddt(
                            A_free= Free_annexin[t]
                            ,A_bound=Bound_annexin[t]
                            ,C=Free_Ca[t]
                            ,pos=pos ,k1=k1 ,k2=k2
                        )
                        Bound_annexin[t+1][y][x] = Bound_annexin[t][y][x] + dt*dAbounddt


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


data_list = Ring_sum(
    ref_Grid=Free_Ca[T_tot-1]
    ,offsets=[x0,y0]
    ,Radius=R   ,dRadius=dR
    ,num_of_rings=3
)
Ring_sums = data_list[0]
Ring_radius = data_list[1]
conc_removed_grid = data_list[2]
visual_rings = data_list[3]

sumfree,sumbound,sumtot = sum_annexin(
    A_free=Free_annexin
    ,A_bound=Bound_annexin
)

timevec = np.linspace(0,T_tot,len(sumfree))

fig, axs = plt.subplots(5)

TT0 = 10
i = 0
axs[i].plot(sumfree[:TT0],'-*',label="free")
axs[i].plot(sumbound[:TT0],label="bound")
axs[i].set_ylim(-max(sumbound[:TT0]),max(sumbound[:TT0])*1.1)
#axs[i].set_ylim(min(sumfree[:TT0]),max(sumfree[:TT0])*1.1)
axs[i].legend()
axs[i].set_title(f"conc of free/bound Annexin"+"\n"
                 +f"k1={k1} , k2={k2} , dt={dt} ,dx={dx}, dy={dy}, D={D}"
                 )

TT0 = 20
i = 1
axs[i].plot(sumfree[:TT0],'-*',label="free")
axs[i].plot(sumbound[:TT0],label="bound")
axs[i].set_ylim(-max(sumbound[:TT0]),max(sumbound[:TT0])*1.1)
#axs[i].set_ylim(min(sumfree[:TT0]),max(sumfree[:TT0])*1.1)
axs[i].legend()

TT0 = 40
i = 2
axs[i].plot(sumfree[:TT0],'-*',label="free")
axs[i].plot(sumbound[:TT0],label="bound")
axs[i].set_ylim(-max(sumbound[:TT0]),max(sumbound[:TT0])*1.1)
#axs[i].set_ylim(min(sumfree[:TT0]),max(sumfree[:TT0])*1.1)
axs[i].legend()

TT0 = 80
i=3
axs[i].plot(sumfree[:TT0],'-*',label="free")
axs[i].plot(sumbound[:TT0],label="bound")
axs[i].set_ylim(-max(sumbound[:TT0]),max(sumbound[:TT0])*1.1)
#axs[i].set_ylim(min(sumfree[:TT0]),max(sumfree[:TT0])*1.1)
axs[i].legend()

TT0 = len(sumfree)
i=4
axs[i].plot(sumfree[:TT0],'-*',label="free")
axs[i].plot(sumbound[:TT0],label="bound")
axs[i].set_ylim(-max(sumbound[:TT0]),max(sumbound[:TT0])*1.1)
#axs[i].set_ylim(min(sumfree[:TT0]),max(sumfree[:TT0])*1.1)
axs[i].legend()
axs[i].set_xlabel("time")

plt.show()
exit()



plt.figure()
plt.plot(Ring_radius,Ring_sums)
plt.title("Ring sums try, normalized by num of points")

plt.matshow(conc_removed_grid)
plt.title("concentration outside inner cell removed")

plt.matshow(visual_rings)
plt.title("Visualizing the rings")

plt.matshow(Free_annexin[0])
plt.title("free annexin start")
plt.matshow(Free_annexin[T_tot-1])
plt.title("free annexin end")

plt.matshow(Bound_annexin[T_tot-1])
plt.title("bound annexin")
plt.show()



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
