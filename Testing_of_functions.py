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
    ref_grid = init_ref_circle(
        boxlen=len_size
        ,Radius=R,dRadius=dR
        ,offsets=[x0,y0]
        ,inside_val=inside_val
        ,outside_val=outside_val
        ,wall_val=wall_val
        )
    
    Ca = init_conc(
        ref_grid=ref_grid
        ,time=T_tot
        ,c_in=c_in,c_out=c_out
        ,inside_val=inside_val
        ,outside_val=outside_val
                        )
    
    #Ca = np.zeros(shape=(T_tot,len_size,len_size))
    #Ca[0][int(len_size/2)][int(len_size/2)] = 100 #middle
    #Ca[0][5][5] = 1 #top left corner
    #Ca[0][int(len_size-5)][5] = 1 # top right corner
    #Ca[0][5][int(len_size-5)] = 1 # bottom left corner
    #Ca[0][int(len_size-5)][int(len_size-5)] = 1 # bottom right corner   
    
    open_hole = True
    Ca_concentration = []
    for t in np.arange(0,T_tot-1): 
        Ca_concentration.append(np.sum(Ca[t]))
        if t%(T_tot/10) == 0:
            print(f"time={t} of {T_tot}")   
        t1, t2 = t%2, (t+1)%2
        t1, t2 = t, t+1
        for x in range(0,len_size):
            for y in range(0,len_size):
                #Bound_Ca[t+1][y][x] += Bound_Ca[t][y][x]
                if ref_grid[y][x] == wall_val or ref_grid[y][x] == outside_val:
                    Ca[t+1][y][x] = Ca[t][y][x]

                if ref_grid[y][x] == inside_val or ref_grid[y][x] == open_val:
                    radii = np.sqrt( (x-x0)**2 + (y-y0)**2 )
                    pos = cicle_boundary(x=x,y=y,boxlen=len_size
                                            ,ref_matrix=ref_grid
                                            ,refval=wall_val
                                        )

                    circle_dCondt(
                        C=Ca,pos=pos
                        ,const=[t,dt,dx,dy,R,dR,radii]
                        ,D=D_Ca_cyto
                    )


    plt.matshow(Ca[0])
    plt.title("init Ca distribution")

    plt.matshow(Ca[T_tot-1])
    plt.title("final Ca distribution")

    plt.figure()
    plt.plot(Ca_concentration)
    plt.title("Ca total concentration over time")

    plt.show()
    

def test_Ca_diff_corner_open_hole():
    ref_grid = init_ref_circle(
        boxlen=len_size
        ,Radius=R,dRadius=dR
        ,offsets=[x0,y0]
        ,inside_val=inside_val
        ,outside_val=outside_val
        ,wall_val=wall_val
        )
    
    Ca = init_conc(
        ref_grid=ref_grid
        ,time=T_tot
        ,c_in=c_in,c_out=c_out
        ,inside_val=inside_val
        ,outside_val=outside_val
                        )
    
    #Ca = np.zeros(shape=(T_tot,len_size,len_size))
    #Ca[0][int(len_size/2)][int(len_size/2)] = 100 #middle
    #Ca[0][5][5] = 1 #top left corner
    #Ca[0][int(len_size-5)][5] = 1 # top right corner
    #Ca[0][5][int(len_size-5)] = 1 # bottom left corner
    #Ca[0][int(len_size-5)][int(len_size-5)] = 1 # bottom right corner   
    
    open_hole = True
    open_close_membrane(
                        Grid=ref_grid
                        ,Radius=R, dRadius=dR
                        ,offsets=[x0,y0],holesize=holesize
                        ,open_val=open_val
                        ,wall_val=wall_val
                        ,open_wall=open_hole
                    )
    
    Ca_concentration = []
    close_time = int(T_tot*0.5)
    for t in np.arange(0,T_tot-1): 
        Ca_concentration.append(np.sum(Ca[t]))
        if t%(T_tot/10) == 0:
            print(f"time={t} of {T_tot}")   
        t1, t2 = t%2, (t+1)%2
        t1, t2 = t, t+1
        for x in range(0,len_size):
            for y in range(0,len_size):
                #Bound_Ca[t+1][y][x] += Bound_Ca[t][y][x]
                if ref_grid[y][x] == wall_val or ref_grid[y][x] == outside_val:
                    Ca[t+1][y][x] = Ca[t][y][x]

                if ref_grid[y][x] == inside_val or ref_grid[y][x] == open_val:
                    radii = np.sqrt( (x-x0)**2 + (y-y0)**2 )
                    pos = cicle_boundary(x=x,y=y,boxlen=len_size
                                            ,ref_matrix=ref_grid
                                            ,refval=wall_val
                                        )

                    circle_dCondt(
                        C=Ca,pos=pos
                        ,const=[t,dt,dx,dy,R,dR,radii]
                        ,D=D_Ca_cyto
                    )
        if t >= close_time and open_hole==True:
                    #ref_structure = 
                    open_close_membrane(
                        Grid=ref_grid
                        ,Radius=R, dRadius=dR
                        ,offsets=[x0,y0],holesize=holesize
                        ,open_val=open_val
                        ,wall_val=wall_val
                        ,open_wall=False
                    )
                    open_hole = False
                    print(f"wall closure time t={t}")

    plt.matshow(Ca[0])
    plt.title("init Ca distribution")

    plt.matshow(Ca[T_tot-1])
    plt.title("final Ca distribution")

    plt.figure()
    plt.plot(Ca_concentration,'r',label="[Ca]")
    plt.vlines(x=close_time,ymin=min(Ca_concentration),ymax=max(Ca_concentration)
               ,label="close hole time")
    plt.title("Ca total concentration over time")
    plt.legend()
    plt.show()


def test_annexin_diff_closed_hole():
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val = constants()
    T_tot= 500
    bound_annexin_start = 0
    c_pump = 0

    
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


    """
    plt.matshow(Free_annexin[0])
    plt.title("init Free annexin distribution")

    plt.matshow(Free_annexin[T_tot-1])
    plt.title("final Free annexin distribution")
    """
    

    
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

    plt.figure()
    plt.plot(A_sumbound,label="simulated")
    plt.ylabel(r"$ \frac{ A_b (t) }{max( A_f )} $")
    plt.xlabel(r"time steps $ c_{pump} $ =0")
    plt.title("Bound annexin")
    plt.legend()
    plt.show()

    
def test_annexin_diff_open_hole():
    ref_grid = init_ref_circle(
        boxlen=len_size
        ,Radius=R,dRadius=dR
        ,offsets=[x0,y0]
        ,inside_val=inside_val
        ,outside_val=outside_val
        ,wall_val=wall_val
        )
    
    A_f = init_conc(
        ref_grid=ref_grid
        ,time=T_tot
        ,c_in=c_in,c_out=c_out
        ,inside_val=inside_val
        ,outside_val=outside_val
                        )
    Ca = np.zeros(shape=(T_tot,len_size,len_size))
    Ca_b = np.zeros(shape=(T_tot,len_size,len_size))
    A_f = np.zeros(shape=(T_tot,len_size,len_size))
    A_b = np.zeros(shape=(T_tot,len_size,len_size))
    A_f[0][int(len_size/2)][int(len_size/2)] = 100 #middle
    A_f[0][5][5] = 1 #top left corner
    A_f[0][int(len_size-5)][5] = 1 # top right corner
    A_f[0][5][int(len_size-5)] = 1 # bottom left corner
    A_f[0][int(len_size-5)][int(len_size-5)] = 1 # bottom right corner   
    
    open_hole = True
    open_close_membrane(
                        Grid=ref_grid
                        ,Radius=R, dRadius=dR
                        ,offsets=[x0,y0],holesize=holesize
                        ,open_val=open_val
                        ,wall_val=wall_val
                        ,open_wall=open_hole
                    )
    
    A_f_conc = []
    A_b_conc = []
    A_tot_conc = []
    close_time = int(T_tot*0.5)
    for t in np.arange(0,T_tot-1): 
        f,b = np.sum(A_f[t])    ,np.sum(A_b[t])
        A_f_conc.append(f)
        A_b_conc.append(b)
        A_tot_conc.append(f+b)

        if t%(T_tot/10) == 0:
            print(f"time={t} of {T_tot}")   
        t1, t2 = t%2, (t+1)%2
        t1, t2 = t, t+1
        for x in range(0,len_size):
            for y in range(0,len_size):
                #Bound_Ca[t+1][y][x] += Bound_Ca[t][y][x]
                if ref_grid[y][x] == wall_val or ref_grid[y][x] == outside_val:
                    Ca[t+1][y][x] = Ca[t][y][x]
                    A_f[t+1][y][x] = A_f[t][y][x]
                    A_b[t+1][y][x] = A_b[t][y][x]


                if ref_grid[y][x] == inside_val or ref_grid[y][x] == open_val:
                    radii = np.sqrt( (x-x0)**2 + (y-y0)**2 )
                    pos = cicle_boundary(x=x,y=y,boxlen=len_size
                                            ,ref_matrix=ref_grid
                                            ,refval=wall_val
                                        )

                    circle_dCondt(
                        C=Ca,pos=pos
                        ,const=[t,dt,dx,dy,R,dR,radii]
                        ,D=D_Ca_cyto
                    )

                    circle_dAdt(
                        A_free= A_f
                        ,A_bound=A_b
                        ,C=Ca
                        ,C_bound=Ca_b
                        ,pos=pos 
                        ,const=[t,dt,dx,dy,k1,k2,R,dR,radii]
                        ,D=D_Annexin_cyto
                    )
        if t >= close_time and open_hole==True:
            open_close_membrane(
                Grid=ref_grid
                ,Radius=R, dRadius=dR
                ,offsets=[x0,y0],holesize=holesize
                ,open_val=open_val
                ,wall_val=wall_val
                ,open_wall=False
            )
            open_hole = False
            print(f"wall closure time t={t}")

    plt.matshow(A_f[0])
    plt.title("init Ca distribution")

    plt.matshow(A_f[T_tot-1])
    plt.title("final Ca distribution")

    plt.figure()
    plt.plot(A_tot_conc,'r',label="[Ca]")
    plt.vlines(x=close_time,ymin=min(A_tot_conc),ymax=max(A_tot_conc)
               ,label="close hole time")
    plt.title("A_tot concentration over time")
    plt.legend()
    plt.show()


def Finding_the_pump_value():
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
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val = constants()
    T_tot= 500
    bound_annexin_start = 0
    c_pump = 0

    
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
    point_count = point_count*dx*dy
    time_vec = np.linspace(start=0,stop=T_tot,num=len(A_f_stability))
    #A_f_stability = [i/max(A_f_stability) for i in A_f_stability]
    A_sumfree = [i/point_count for i in A_sumfree]
    #A_b_stability = [i/max(A_b_stability) for i in A_b_stability]
    A_sumbound = [i/point_count for i in A_sumbound]
    
    fig1 = go.Figure()
    fig1.add_trace(
        go.Scatter(
            x=time_vec
            ,y=A_sumfree#[i/max(A_sumfree) for i in A_sumfree]
            ,mode="lines + markers"
            ,name="simulated"
        )
    )
    fig1.update_xaxes(title_text="timesteps")
    fig1.update_yaxes(title_text="[A_free]")
    fig1.update_yaxes(
        exponentformat="SI"
    )
    fig1.show()

    fig2 = go.Figure()
    fig2.add_trace(
        go.Scatter(
            x=time_vec
            ,y=A_f_stability#[i/max(A_f_stability) for i in A_f_stability]
            ,mode="lines + markers"
            ,name="Analytical sol"
        )
    )
    fig2.update_xaxes(title_text="timesteps")
    fig2.update_yaxes(title_text="[A_free]")
    fig2.show()

    fig3 = go.Figure()
    fig3.add_trace(
        go.Scatter(
            x=time_vec
            ,y=A_sumtot
            ,name="total sum"
        )
    )
    fig3.show()

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
    #plt.plot(A_sumfree,label="simulated")
    plt.xlabel("time steps")
    plt.ylabel(r" $ \frac{A_f(t)}{max(A_f)} $ ")
    plt.title(r"Free annexin  $ c_{pump} $ =0")
    plt.legend()

    plt.figure()
    #plt.plot(time_vec,A_f_stability,label="equation")
    plt.plot(A_sumfree,label="simulated")
    plt.xlabel("time steps")
    plt.ylabel(r" $ \frac{A_f(t)}{max(A_f)} $ ")
    plt.title(r"Free annexin  $ c_{pump} $ =0")
    plt.legend()

    print("max abstab=",max(A_b_stability))
    
    plt.figure()
    plt.plot(time_vec,A_b_stability,label="equation")
    #plt.plot(A_sumbound,label="simulated")
    plt.ylabel(r"$ \frac{ A_b (t) }{max( A_f )} $")
    plt.xlabel(r"time steps $ c_{pump} $ =0")
    plt.title("Bound annexin")
    plt.legend()
    
    plt.figure()
    #plt.plot(time_vec,A_b_stability,label="equation")
    plt.plot(A_sumbound,label="simulated")
    plt.ylabel(r"$ \frac{ A_b (t) }{max( A_f )} $")
    plt.xlabel(r"time steps $ c_{pump} $ =0")
    plt.title("Bound annexin")
    plt.legend()

    point_sum_Af = []

    for t in range(T_tot):
        point_sum_Af.append(
            Free_annexin[t][int(len_size/2)][int(len_size/2)]*dx*dy
        )
    a = abs(point_sum_Af[len(point_sum_Af)-1] - A_f_stability[len(A_f_stability)-1])
    print(f"A_sumtot : max-min sim ={max(A_sumtot)-min(A_sumtot)}")
    print(f"diff at end of pointsumAf and afstabi = {a}")

    plt.figure()
    plt.plot(point_sum_Af,'*',label="point sum sim")
    plt.plot(time_vec,A_f_stability,label="analytical sol")
    plt.title(f"point sum Free annexin x,y={int(len_size/2)}")
    plt.legend()

    plt.show()


if __name__ =="__main__":
    #test_reference_struct()
    #test_Ca_diff_corner_closed_hole()
    #test_Ca_diff_corner_open_hole()
    #test_annexin_diff_closed_hole() # also tests the analytical sol for dA_f and dA_b
    #test_annexin_diff_open_hole()
    #Finding_the_pump_value()
    test_analytical_vs_sim_dAf_dAb()
    pass