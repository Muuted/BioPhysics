from Circle_funcs import stabil_condi
import matplotlib.pyplot as plt
from Circle_funcs import Annexin_stablilization
import numpy as np
def constants():
    Avogadro = 6.02214076e23 # 1/mol
    # time and step size, and diffusion constant
    picture_size = 83e-6 # meters
    Real_sim_time = 120 #seconds
    real_close_time = 10 #s
    
    len_size = 80 # number of grid points
    dx, dy = picture_size/len_size ,picture_size/len_size # meters

    # Concentrations for Ca^2+ inside; outside cell
    c_out =2e-3# 2*10**(-3)#e-3 #M Concentration outside of cell
    c_in = 100e-9#*10*(-9)#e-9 #M Concetration inside cell

    D_Ca_cyto  = 2.7e-11#*10**(-11)#e-11 #meters^2/second
    D_Annexin_cyto  = 5.0e-11#e-11 #meters^2/second

    dt = stabil_condi(dt=0.1,dx=dx,dy=dy
                    ,D_list=[ D_Ca_cyto , D_Annexin_cyto ]
                    )
    
    
    T_tot = int(Real_sim_time/dt)

    pump_multiplyer = 30 # found from sim testing.
    c_pump = dt*(5e-8)*pump_multiplyer # Molar/second


    # Anxexin constants.  
    A_total_density = ((2e6/10))/(1e-18) #number/meter^3 from article directly
    A_total_conc = (A_total_density/Avogadro)*(0.1/100)#should get 0.1% of the 3D volume concentration for all proteins

    k1 = 30e3 # 1/(Ms)
    k2 = 1e-3 #1/s

    
    bound_annexin_start = k1*A_total_conc*c_in/(k1*c_in + k2)
    c_in_annexin = A_total_conc - bound_annexin_start 
     

    #Closing mechanism time
    close_time = int(real_close_time/dt)

    #size of cell, hole in cell and the center's placement.
    holesize = 1
    dR = int(2)
    R = int(len_size*0.5 - 3*dR) 
    
    x0,y0 = int(len_size/2), int(len_size/2)

    # Defining areas of interest, values for regonition.
    wall_val = 100
    inside_val = 10
    outside_val = -5
    open_val = 20

    print(
        f" \n \n"
        + "------------- Constant used in Simulation -------------------- \n "
        + f"    Real simulation time = {Real_sim_time} s \n "
        + f"    hole close time = {int(close_time*dt)} s \n " 
        + f"    dx=dy = {dx:e} m \n "
        + f"    dt = {dt:e} m \n "
        + f"    [Ca] outside cell = {c_out:e} \n "
        + f"    [Ca] inside cell = {c_in:e} \n "
        + f"    Diffuse constant Ca in cell = {D_Ca_cyto} m^2/s \n "
        + f"    Diffuse constant Annexin in cell = {D_Annexin_cyto} m^2/s \n "
        + f"    c_pump = {c_pump} M/s \n "
        + f"    Total initial conc Annexins = {A_total_conc:e} M \n "
        + f"    c_in_ann = {c_in_annexin:e} M \n "
        + f"    and the bound_ann_start = {bound_annexin_start:e} M \n "
        + f"    k1 = {k1:e} 1/Ms \n "
        + f"    and k2 = {k2:e} 1/s \n "
        + f"    total number of sim steps = {T_tot:e} steps \n "
        #+ f"    sim equil time = {1/(k1*c_in_annexin + k2)} \n "
        + f"    size of grid = ({len_size},{len_size}) \n "
        f" ------------------------------------------------------------- \n \n "
    )

    args_list = [
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2
        ,c_in_annexin
        ,bound_annexin_start
        ,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,Real_sim_time, real_close_time
    ]
    return args_list






if __name__ == "__main__":
    A = constants()