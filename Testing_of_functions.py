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

def test_reference_struct():
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val = constants()
    ref_grid = init_ref_circle(
        boxlen=len_size
        ,Radius=R,dRadius=dR
        ,offsets=[x0,y0]
        ,inside_val=inside_val
        ,outside_val=outside_val
        ,wall_val=wall_val
        )
    Test_truth = True
    N =5
    for i in range(N):
        if i%2 == 0:
            Test_truth = True
        if i%2 != 0:
            Test_truth = False
        open_close_membrane(
            Grid=ref_grid
            ,Radius=R,dRadius=dR
            ,offsets=[x0,y0]
            ,holesize=holesize
            ,open_val=open_val
            ,wall_val=wall_val
            ,open_wall=Test_truth
            )
        
        plt.matshow(ref_grid)
        plt.title("Testing that we can change back and forth"
                  +"\n"+f" of an open structure "+"\n"+"and closed structure "f"i={i} of {N}"
                  )
        plt.show()


def test_Ca_diff_corner_closed_hole():
    print("Testing diffusion of Calcium, with closed cell")
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val = constants()
    c_pump = 0
    #T_tot = 500
    c_in_annexin = 0
    bound_annexin_start = 0
    Sim_data_list = main_circle_sim(
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,open_hole=False
                    )
    ref_structure,Free_Ca,Free_annexin,Bound_annexin,Bound_Ca = Sim_data_list

    Ca_sumfree,Ca_sumbound,Ca_sumtot = sum_annexin(
        A_free=Free_Ca
        ,A_bound=Bound_Ca
    )

    plt.matshow(Free_Ca[0])
    plt.title(r"init $ [ Ca_{tot} ] $ distribution")

    plt.matshow(Free_Ca[T_tot-1])
    plt.title(r"final $ [ Ca_{tot} ] $distribution")

    plt.figure()
    plt.plot(Ca_sumfree)
    plt.title(r"$ [ Ca_{tot} ] $ over time, for closed hole")

    plt.show()
    

def test_Ca_diff_corner_open_hole():
    print("Testing calcium diffusive with open cell")
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val = constants()
    #c_pump = 0
    #T_tot = 1500
    #T_tot = int(30/dt)
    c_in_annexin = 0
    bound_annexin_start = 0

    Sim_data_list = main_circle_sim(
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,open_hole=True
                    )
    ref_structure,Ca,Free_annexin,Bound_annexin,Bound_Ca = Sim_data_list

    Ca_sumfree,Ca_sumbound,Ca_sumtot = sum_annexin(
        A_free=Ca
        ,A_bound=Bound_Ca
    )

    plt.matshow(Ca[0])
    plt.title("init Ca distribution")

    plt.matshow(Ca[T_tot-1])
    plt.title("final Ca distribution")

    plt.figure()
    plt.plot(Ca_sumfree,'r',label="[Ca]")
    plt.vlines(x=close_time,ymin=min(Ca_sumfree),ymax=max(Ca_sumfree),label=f"close hole time={close_time}")
    plt.title(r"$ [ Ca_{tot} ] $ over time, openhole, and $ C_{pump} = $" +f"{c_pump}"
              + "\n"+ f"holeclosed at t={close_time}steps and realtime={int(T_tot*dt)}s")
    
    plt.legend()
    plt.show()


def test_annexin_diff_closed_hole():
    print("Testing Annexin diffusion, closed cell")
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val = constants()
    #T_tot= 1500
    bound_annexin_start = 0
    c_pump = 0
    close_time = 0
    c_in = c_out
    T_tot = 1600
    Sim_data_list = main_circle_sim(
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,open_hole=False
                    )
    ref_structure,Free_Ca,Free_annexin,Bound_annexin,Bound_Ca = Sim_data_list

    A_sumfree,A_sumbound,A_sumtot = sum_annexin(
        A_free=Free_annexin
        ,A_bound=Bound_annexin
    )

    print(f"maximum free annexin ={max(A_sumfree)}")
    
    plt.figure()
    plt.plot(A_sumtot,'r',label="[Annexin]")
    plt.vlines(
        x=close_time,ymin=min(A_sumtot),ymax=max(A_sumtot)
        ,label=f"close hole time={close_time}"
               )
    plt.xlabel("time steps")
    plt.ylabel(r" $ A_f + A_b $ ")
    plt.title(r" $ A_{tot} $ and  $ c_{pump} $ =0")
    plt.legend()

    plt.figure()
    plt.plot(A_sumbound,label="simulated")
    plt.ylabel(r"$ A_b (t)  $")
    plt.xlabel(r"time steps $ c_{pump} $ =0")
    plt.title("Bound annexin")
    plt.legend()

    plt.figure()
    plt.plot(A_sumfree,label="simulated")
    plt.ylabel(r"$ A_f (t) $")
    plt.xlabel(r"time steps $ c_{pump} $ =0")
    plt.title("Free annexin")
    plt.legend()


    plt.show()

    
def test_annexin_diff_open_hole():
    print("Testing Annexin diffusion, open cell")
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val = constants()
    #T_tot= 1500
    bound_annexin_start = 0
    c_pump = 0
    
    Sim_data_list = main_circle_sim(
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,open_hole=True
                    )
    ref_structure,Free_Ca,Free_annexin,Bound_annexin,Bound_Ca = Sim_data_list

    A_sumfree,A_sumbound,A_sumtot = sum_annexin(
        A_free=Free_annexin
        ,A_bound=Bound_annexin
    )

    plt.matshow(Free_annexin[0])
    plt.title("init Annexin distribution")

    plt.matshow(Free_annexin[T_tot-1])
    plt.title("final Annexin distribution")

    plt.figure()
    plt.plot(A_sumtot,'r',label="[Ca]")
    plt.vlines(x=close_time,ymin=min(A_sumtot),ymax=max(A_sumtot)
               ,label="close hole time")
    plt.title("A_tot concentration over time")
    plt.legend()
    plt.show()


def Finding_the_pump_value():
    print(r"Finding the value for the $ c_{pump} $ term")
    from Circle_sim import main_circle_sim
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val = constants()
    Real_time = 80 #seconds.
    T_tot = int(Real_time/dt) # number of time steps.
    #c_pump = c_pump*30
    close_time = T_tot*0.1
    
    time1 = tm.time()
    Sim_data_list = main_circle_sim(
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
                    )
    time2 = tm.time()

    print(f"Sim took {(time2-time1)/60}min")
    ref_structure = Sim_data_list[0]
    Free_Ca = Sim_data_list[1]

    conc_Ca_time = sum_in_cell(
        ref_Grid=ref_structure
        ,Matrix=Free_Ca
        ,inside_val=inside_val
                               )

    real_time_vec = np.linspace(0,Real_time,T_tot)

    plt.figure()
    plt.plot(real_time_vec,conc_Ca_time/max(conc_Ca_time))
    plt.title(f"[Ca] in cell, " + r"$ c_{pump} $" +f"={c_pump}")
    plt.xlabel("time [s]")
    plt.ylabel("Concentration/max(conc)")
    plt.show()
    

def test_analytical_vs_sim_dAf_dAb():
    print(r"Testing the analytical vs simulated evolution of $ A_f (t) $ and $ A_b (t) $")
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val = constants()
    T_tot= 500
    #bound_annexin_start = 0
    #c_pump = 0
    close_time = 0


    Sim_data_list = main_circle_sim(
        c_in,c_out,D_Ca_cyto,T_tot,len_size
        ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto
        ,dt,close_time,c_pump,holesize,dR,R,x0,y0
        ,wall_val,inside_val,outside_val,open_val
        ,open_hole=False
                    )
    
    ref_structure,Free_Ca,Free_annexin,Bound_annexin,Bound_Ca = Sim_data_list

    A_sumfree,A_sumbound,A_sumtot = sum_annexin(
        A_free=Free_annexin
        ,A_bound=Bound_annexin
    )

    print(f"maximum free annexin ={max(A_sumfree)}")
    A_f_stability,A_b_stability = Annexin_stablilization(
        k1=k1,k2=k2
        ,A_fo=c_in_annexin
        ,realtime=int(T_tot*dt)
        ,c_in=c_in
        ,dt=dt
    )
    point_count = points_inside_cell(refgrid=ref_structure,inside_val=inside_val)
    #initial_free_Annexin = np.sum(Free_annexin[0])

    #point_count = point_count*dx*dy
    time_vec = np.linspace(start=0,stop=T_tot,num=len(A_f_stability))
    #A_f_stability = [i/max(A_f_stability) for i in A_f_stability]
    A_sumfree = [i/point_count for i in A_sumfree]
    #A_b_stability = [i/max(A_b_stability) for i in A_b_stability]
    A_sumbound = [i/point_count for i in A_sumbound]
    
    
    plt.figure()
    plt.plot(A_sumtot,'r',label="[Annexin]")
    plt.vlines(
        x=close_time,ymin=min(A_sumtot),ymax=max(A_sumtot)
        ,label="close hole time"
               )
    plt.xlabel("time steps")
    plt.ylabel(r" $ A_f + A_b $ ")
    plt.title(r" $ A_{tot} $ and  $ c_{pump} $ =0")
    plt.legend()

    
    print(f"area of cell ={point_count}")

    plt.figure()
    plt.plot(time_vec,A_f_stability,label="equation")
    plt.plot(A_sumfree,label="simulated")
    plt.xlabel("time steps")
    plt.ylabel(r" $ \frac{A_f(t)}{max(A_f)} $ ")
    plt.title(r"Free annexin  $ c_{pump} $ =0")
    plt.legend()


    print("max abstab=",max(A_b_stability))
    
    plt.figure()
    plt.plot(time_vec,A_b_stability,label="equation")
    plt.plot(A_sumbound,label="simulated")
    plt.ylabel(r"$ \frac{ A_b (t) }{max( A_f )} $")
    plt.xlabel(r"time steps $ c_{pump} $ =0")
    plt.title("Bound annexin")
    plt.legend()
    

    point_sum_Af = []

    for t in range(T_tot):
        point_sum_Af.append(
            Free_annexin[t][int(len_size/2)][int(len_size/2)]#*dx*dy
        )

    a = abs(point_sum_Af[len(point_sum_Af)-1] - A_f_stability[len(A_f_stability)-1])
    print(f"A_sumtot : max-min sim ={max(A_sumtot)-min(A_sumtot)}")
    print(f"diff at end of pointsumAf and afstabi = {a}")

    plt.figure()
    plt.plot(point_sum_Af,'*',label="point sum sim")
    plt.plot(time_vec,A_f_stability,label="analytical sol")
    plt.title(f"point sum Free annexin x,y={int(len_size/2)}")
    plt.legend()

    diff_Af,diff_Ab = [], []
    for i in range(len(A_f_stability)):
        diff_Af.append(
            A_f_stability[i] - A_sumfree[i]
        )
        diff_Ab.append(
            A_b_stability[i] - A_sumbound[i]
        )

    plt.figure()
    plt.plot(diff_Af)
    plt.title(r"$ A_{f,equation} - A_{f,sim} $")

    plt.figure()
    plt.plot(diff_Ab)
    plt.title(r"$ A_{b,equation} - A_{b,sim} $")
    
    plt.show()


if __name__ =="__main__":
    #test_reference_struct()
    #test_Ca_diff_corner_closed_hole()
    #test_Ca_diff_corner_open_hole()
    #test_annexin_diff_closed_hole()
    #test_annexin_diff_open_hole()
    #Finding_the_pump_value()
    test_analytical_vs_sim_dAf_dAb()
    