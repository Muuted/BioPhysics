from Circle_funcs import stabil_condi

# Concentrations for Ca^2+ inside; outside cell
c_out = 2e-3 #M Concentration outside of cell
c_in = 100e-9 #M Concetration inside cell
D_Ca_cyto  = 2.7e-11 #meters^2/second
#D_Ca_water  = 7.93e-10 #meters^2/second


# time and step size, and diffusion constant
picture_size = 83e-6 # meters
Real_sim_time = 30 #seconds

T_tot = 6000
len_size = 40 # number of grid points
dx, dy = picture_size/len_size ,picture_size/len_size # meters

# Anxexin constants.
k1 ,k2 = 0.1 ,0.1
c_in_annexin = 1e-9
bound_annexin_start = 0#k1*c_in_annexin*c_in/k2
A_b_init = k1*c_in_annexin*c_in/k2

D_Annexin_cyto  = 5.0e-11 #meters^2/second
#D_Annexin_water  = 7.18e-11 #meters^2/seconds



dt = stabil_condi(dt=0.1,dx=dx,dy=dy
                  ,D_list=[D_Ca_cyto,D_Annexin_cyto]
                  )
print("dt=",dt)

#T_tot = int(Real_sim_time/dt)
print("total number of sim steps =",T_tot)


# mechanisms change.
close_time = T_tot*0.9

c_pump = c_in*dt #5e-8*dt # meters/second
print(f"c_pump={c_pump}")
#size of cell, hole in cell and the center's placement.
holesize = 3
dR = int(2)
R = int(len_size*0.5 - 3*dR)

x0,y0 = int(len_size/2), int(len_size/2)
wall_val = 100
inside_val = 10
outside_val = -5
open_val = 20
