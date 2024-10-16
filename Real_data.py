import matplotlib.pyplot as plt
import numpy as np
from Constants import constants
data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\GCaMP6s-CAAX (PM anchor)"


if __name__ == "__main__":
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val = constants()

    c = 0 #k2*k1*c_in-k1*c_in
    b= -k1*c_in-k2
    print(f"b={b}")
    a = 1

    d = (b)**2 -4*1*(c)

    print(f"the auxiliary eq  d={d} \n finding the 2 roots")

    root1 = (-b + np.sqrt(d))/(2*a)

    root2 = (-b - np.sqrt(d))/(2*a)

    print(f" root1 ={root1}")
    print(f" root2 ={root2}")

    




