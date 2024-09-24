from Circle_funcs import stabil_condi

# Concentrations for Ca^2+ inside; outside cell
c_out = 2e-3 #M Concentration outside of cell
c_in = 100e-9 #M Concetration inside cell

# time and step size, and diffusion constant
picture_size = 83e-6 # meters
T_tot = 100
len_size = 40 # number of grid points
dx, dy = picture_size/len_size ,picture_size/len_size # meters
D_Ca_cyto  = 2.7e-11 #meters^2/second
D_Ca_water  = 2.7e-11 #meters^2/second
dt = stabil_condi(dt=0.1,dx=dx,dy=dy,D=D_Ca_cyto)

c_pump = 5e-7 # meters/second

# mechanisms change.
close_time = T_tot*0.9

# Anxexin constants.
k1 ,k2 = 0.1 ,0.1
c_in_annexin = c_in
prob_free_ann = 0.5
prob_bound_ann = 0.5
D_Annexin_cyto  = 2.7e-11 #meters^2/second
D_Annexin_water  = 2.7e-11 #meters^2/seconds

#size of cell, hole in cell and the center's placement.
holesize = 3
R = int(len_size*0.5/2)
dR = int(2)
x0,y0 = int(len_size/2), int(len_size/2)
wall_val = 100
inside_val = 10
outside_val = -5
open_val = 20