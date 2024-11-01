from Circle_funcs import stabil_condi

def constants():
    Avogadro = 6.02214076*10**(23) # 1/mol
    # time and step size, and diffusion constant
    picture_size = 83e-6#*10**(-6)#e-6 # meters
    Real_sim_time = 30 #seconds

    T_tot = 1600
    len_size = 50 # number of grid points
    dx, dy = picture_size/len_size ,picture_size/len_size # meters

    # Concentrations for Ca^2+ inside; outside cell
    c_out =2e-3# 2*10**(-3)#e-3 #M Concentration outside of cell
    c_in = 100e-9#*10*(-9)#e-9 #M Concetration inside cell
    D_Ca_cyto  = 2.7e-11#*10**(-11)#e-11 #meters^2/second
    #D_Ca_water  = 7.93e-10 #meters^2/second
    D_Annexin_cyto  = 5.0*10**(-11)#e-11 #meters^2/second
    #D_Annexin_water  = 7.18e-11 #meters^2/seconds
    dt = stabil_condi(dt=0.1,dx=dx,dy=dy
                    ,D_list=[ D_Ca_cyto , D_Annexin_cyto ]
                    )
    print("dt=",dt)
    c_pump = dt*30*(5e-8) # meters/second
    print(f"c_pump={c_pump}")

    
    # Anxexin constants.
        # older numbers
    #A_init_density_upper = (4*10**(6))/(10**(-18)) #number/meter^3 from article directly
    #A_init_density_lower = ((2/10)*10**(6))/(10**(-18)) #number/meter^3 from article directly

#    A_init_density_upper = (4*10**(6))/(10**(-18)) #number/meter^3 from article directly
 #   A_init_density_lower = ((2/10)*10**(6))/(10**(-18)) #number/meter^3 from article directly
  #  A_init_density_upper = (4e6)/(1e-18) #number/meter^3 from article directly
    A_init_density_lower = ((2e6/10))/(1e-18) #number/meter^3 from article directly
   # A_init_conc_upper = (A_init_density_upper/Avogadro)*(0.1/100)#should get 0.1% of the 3D volume concentration for all proteins
    A_init_conc_lower= (A_init_density_lower/Avogadro)*(0.1/100)#should get 0.1% of the 3D volume concentration for all proteins
    #print(f"Upper conc of Annexins in start ={A_init_conc_upper}")
    print(f"Lower conc of Annexins in start ={A_init_conc_lower}")


    #k2 = c_pump

    #K_low = k2/c_out
    #K_up = k2/c_in
    k1 = 30e3 # 1/(Ms)
    k2 = 1e-3 #1/s
    #print(f" k2/c_out={K_low} << k1 << k2/c_in={K_up}")
    #k1 = k2*(c_in + c_out)/(2*c_in*c_out)
    print(f"k1 ={k1} \n and k2={k2}")
    c_in_annexin = A_init_conc_lower # 1e-9 this was our initial guess.
    bound_annexin_start = k1*c_in_annexin*c_in/k2
    A_b_init = k1*c_in_annexin*c_in/k2
    #T_tot = int(Real_sim_time/dt)
    print("total number of sim steps =",T_tot)
    

    # mechanisms change.
    close_time = T_tot*0.5

    #size of cell, hole in cell and the center's placement.
    holesize = 3
    dR = int(2)
    R = int(len_size*0.5 - 3*dR) 

    x0,y0 = int(len_size/2), int(len_size/2)
    wall_val = 100
    inside_val = 10
    outside_val = -5
    open_val = 20

    args_list = [
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2
        ,c_in_annexin
        ,bound_annexin_start
        ,A_b_init
        ,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
    ]
    return args_list

if __name__ == "__main__":
    constants()