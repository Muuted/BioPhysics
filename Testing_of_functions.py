import matplotlib.pyplot as plt
import numpy as np
from Constants import *
from Circle_funcs import *
from Data_extraction_funcs import *
from Constants import constants
import time as tm

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
    #Ca[0][5][5] = 1 #top left corner
    #Ca[0][int(len_size-5)][5] = 1 # top right corner
    #Ca[0][5][int(len_size-5)] = 1 # bottom left corner
    #Ca[0][int(len_size-5)][int(len_size-5)] = 1 # bottom right corner   
    
    open_hole = False
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
    c_pump = c_pump*30
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
    plt.title("conc of Ca over time, inside the cell")
    plt.xlabel("time [s]")
    plt.ylabel("Concentration/max(conc)")
    plt.show()
    

if __name__ =="__main__":
    #test_reference_struct()
    #test_Ca_diff_corner_closed_hole()
    #test_Ca_diff_corner_open_hole()
    #test_annexin_diff_closed_hole()
    #test_annexin_diff_open_hole()
    Finding_the_pump_value()