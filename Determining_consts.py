import matplotlib.pyplot as plt
import numpy as np
from Constants import *
from Circle_funcs import *
from Data_extraction_funcs import *
from Constants import constants
from Circle_sim import main_circle_sim
import time as tm
import plotly.graph_objects as go
import pandas as pd


def Annexin_equil_ratio():
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val,Real_sim_time, real_close_time = constants()
    
    num_steps = 100
    A_tot = c_in_annexin + bound_annexin_start
    realtime = 3000

    Ca_linspace = np.linspace(c_in,c_out,num_steps)

    Annexin_ratio = []
    for c_in in Ca_linspace:
        A_f, A_b = Annexin_stablilization(
            k1=k1,k2=k2
            ,A_tot=A_tot
            ,c_in=c_in
            ,realtime=realtime,dt=dt
            )
        
        A_f_end = A_f[len(A_f)-1]
        A_b_end = A_b[len(A_b)-1]
        
        Annexin_ratio.append(A_b_end/A_f_end)
    

    plt.plot(Ca_linspace,Annexin_ratio)
    plt.show()

if __name__ == "__main__":
    
    Annexin_equil_ratio()
