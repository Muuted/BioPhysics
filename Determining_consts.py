import matplotlib.pyplot as plt
import numpy as np
from Constants import *
from Circle_funcs import *
from Data_extraction_funcs import *
from Constants import constants
import time as tm
import os
import pandas as pd
from Ring_sum_file import main_ring_summing
from Make_movie import Make_video2
from Compare_data import main_compare
from Circle_sim import main_circle_sim

def Finding_the_pump_value():
    print(r"Finding the value for the $ c_{pump} $ term")
    c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val,Real_sim_time, real_close_time = constants()
    Real_time = 80 #seconds.
    ref_time = int(55/dt)
    T_tot = int(Real_time/dt) # number of time steps.
    c_in_annexin = 0
    bound_annexin_start = 0
    time1 = tm.time()
    real_time_vec = np.linspace(0,Real_time,T_tot)
    d_c_pump = 0.5 # change of the pump value
    to_small = False
    to_big = False
    N = 100
    for i in range(N):
        print(
            f"\n ----------------------- \n"
            +f" we are doing the i={i} of {N-1}  \n "
            +f" and c_pump={c_pump}"
              )

        Sim_data_list = main_circle_sim(
            c_in,c_out,D_Ca_cyto,T_tot,len_size
            ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto
            ,dt,close_time,c_pump,holesize,dR,R,x0,y0
            ,wall_val,inside_val,outside_val,open_val
                        )
        time2 = tm.time()

        print(f" \n Sim took {(time2-time1)%60} min \n ")

        ref_structure = Sim_data_list[0]
        Free_Ca = Sim_data_list[1]

        conc_Ca_time = sum_in_cell(
            ref_Grid=ref_structure
            ,Matrix_Free=Free_Ca
            ,Matrix_Bound= np.zeros(shape=(T_tot,len_size,len_size))
            ,inside_val=inside_val
                                    )
        

        if conc_Ca_time[ref_time]/max(conc_Ca_time) < 0.1:
            c_pump *=  (1 - d_c_pump)
            to_small = True
            print(f"\n ----------------------- \n")
            print(f"to small   c_pump={c_pump}")


        if conc_Ca_time[ref_time]/max(conc_Ca_time) > 0.2:
            c_pump *= (1 + d_c_pump)
            to_big = True
            print(f"\n \n ----------------------- \n \n")
            print(f"to big   c_pump={c_pump}")


        if to_big and to_small == True:
            d_c_pump *= 0.5
            to_big = False
            to_small = False
            print(f"we changed the d_c_pump")               
        
        if 0.1 < conc_Ca_time[ref_time]/max(conc_Ca_time)  < 0.2:
            print(f"\n ----------------- ------ \n")
            print(f"We found the c_pump value it is c_pump = {c_pump}")
            print(f"\n ----------------- ------ \n")

            df = pd.DataFrame({
                'graph': [conc_Ca_time],
                'Found C_pump': c_pump
            })

            fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
            fig_folder_path =  fig_save_path #+ f"simtime={Real_sim_time}\\"
            df_name = "Found C_pump value.pkl"
            df.to_pickle(fig_save_path + df_name)
            break
            
        if i == N-1:
            print(f"\n \n ----------------- ------ \n \n")
            print("Never found the correct value")
            print(f"\n \n ----------------- ------ \n \n")
        
        #print(f"value at {len(conc_Ca_time)-1} s, conc_Ca_time[time_compare] = {conc_Ca_time[len(conc_Ca_time)-1]}")

    plt.figure()
    plt.plot(real_time_vec,conc_Ca_time/max(conc_Ca_time),label= r" Total [ $ Ca^{2+} $ ] ")
    plt.plot(55,0.2,label="target",marker='o',color='red')
    plt.title(
        r"Total [ $ Ca^{2+} $ ] in cell "
        +f" and " 
        +r"$ c_{pump} $ $ \approx $ " 
        + f"{c_pump:e}"
        ,fontsize=15
        )
    plt.xlabel("time [s]",fontsize=15)
    plt.ylabel(
        r" $ \frac{ [Ca]_{tot} }{ max( [Ca]_{tot} ) } $ "
        +f"              "
        ,rotation='horizontal'
        ,fontsize=20#,y=0.45
    )
    plt.legend()
    plt.show()


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



def Find_pump_via_area():
    print(r"Finding the value for the $ c_{pump} $ term")

    const_list = constants()
    c_in ,c_out, D_Ca_cyto, T_tot, len_size, dx, dy, k1, k2 = const_list[0:9]
    c_in_annexin ,bound_annexin_start ,D_Annexin_cyto = const_list[9:12]
    dt ,close_time, c_pump, holesize ,dR ,R ,x0 ,y0 = const_list[12:20]
    wall_val ,inside_val ,outside_val ,open_val = const_list[20:24]
    Real_sim_time, real_close_time = const_list[24:26]
    ref_struct_name_cell ,fig_save_path = const_list[26:28]
    fig_folder_path ,video_save_path ,fig_name_df, data_path = const_list[28:32]
    Ca_data_exp ,Annexin_data_exp = const_list[32:34]

    Real_time = 80 #seconds.
    ref_time = int(55/dt)
    T_tot = int(Real_time/dt) # number of time steps.
    c_in_annexin = 0
    bound_annexin_start = 0
    time1 = tm.time()
    real_time_vec = np.linspace(0,Real_time,T_tot)
    d_c_pump = 0.5 # change of the pump value
    to_small = False
    to_big = False
    N = 100
    for i in range(N):
        print(
            f"\n ----------------------- \n"
            +f" we are doing the i={i} of {N-1}  \n "
            +f" and c_pump={c_pump}"
              )

        """Sim_data_list = main_circle_sim(
            c_in,c_out,D_Ca_cyto,T_tot,len_size
            ,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,D_Annexin_cyto
            ,dt,close_time,c_pump,holesize,dR,R,x0,y0
            ,wall_val,inside_val,outside_val,open_val
                        )
        ref_structure,Free_Ca,Free_annexin,Bound_annexin,Bound_Ca= Sim_data_list
        

        df = pd.DataFrame({
            'Free Calcium': [Free_Ca],
            'Bound Calcium': [Bound_Ca],
            'Free Annexins': [Free_annexin],
            'Bound Annexins': [Bound_annexin],
            'Reference Matrix': [ref_structure],
            'hole position X': x0 ,
            'hole position Y': y0,
            'Sim time (s)': Real_sim_time,
            'time steps': T_tot,
            'k1': k1,
            'k2' :k2,
            'dt': dt,
            'free annexin start': c_in_annexin,
            'bound annexin start': bound_annexin_start,
            'D_Ca': D_Ca_cyto,
            'D_annexin': D_Annexin_cyto,
            'hole closure time': close_time,
            'c_pump': c_pump,
            'hole size': holesize,
            'dx':dx,
            'dy':dy,
            'side length':len_size,
            'inside val': inside_val,
            'outside val': outside_val,
            'wall val': wall_val
                        })
        
        if not os.path.exists(fig_folder_path):
            os.makedirs(fig_folder_path)

        
        df.to_pickle(fig_folder_path + fig_name_df)
        
        main_ring_summing(
            fig_save_path=fig_save_path
            ,fig_folder_path=fig_folder_path
            ,df_name=fig_name_df
            ,hole_pos= ""
            ,sum_quick=True
        )
        """
        """ --------- Loading Simulated data ----------"""
        df_sim = pd.read_pickle(fig_folder_path + fig_name_df)
        print(df_sim.info())
        for i in range(5):
            a = df_sim['Ring sum list Ca'][i]
            if isinstance(a , np.ndarray):
                sim_ring_data_Ca = df_sim['Ring sum list Ca'][i]
                sim_ring_data_Ann = df_sim['Ring sum list Annexin'][i]
                print(f"i = {i}")
                break

        ref_structure = df_sim['Reference Matrix'][0]
        Free_Ca = df_sim['Free Calcium'][0]
        Bound_Ca = df_sim['Bound Calcium'][0]

        conc_over_time_Free_Ca = sum_in_cell(
            ref_Grid=ref_structure
            ,Matrix_Free=Free_Ca
            ,Matrix_Bound=Bound_Ca
            ,inside_val=inside_val
                                 )
        
        """  --------- Experimental data loaded  -------------- """

        data_Ca = pd.read_csv(data_path + Ca_data_exp)
        data_Ann= pd.read_csv(data_path + Annexin_data_exp)
        
        exp_data_shape_t, exp_data_shapeX = np.shape(data_Ca)        

        real_data_Ca = np.zeros(shape=(exp_data_shape_t,exp_data_shapeX))
        real_data_Ann = np.zeros(shape=(exp_data_shape_t,exp_data_shapeX))

        for t in range(exp_data_shape_t):
            for R in range(exp_data_shapeX):
                real_data_Ca[t][R] = data_Ca.loc[t][R]
                real_data_Ann[t][R] = data_Ann.loc[t][R]




        ca_sim = []
        ca_exp = []
        for t in range(234):
            sim_ca_conc=0
            exp_ca_conc =0
            for ring in range(10):
                sim_ca_conc += sim_ring_data_Ca[t][ring]
                exp_ca_conc += real_data_Ca[t][ring]
            ca_sim.append(sim_ca_conc)
            ca_exp.append(exp_ca_conc)
            
        vec = np.linspace(
            start=0
            ,stop = len(ca_sim)
            ,num= len(conc_over_time_Free_Ca)
        )
        conc_over_time_Free_Ca = [conc_over_time_Free_Ca[t] - conc_over_time_Free_Ca[0] for t in range(len(conc_over_time_Free_Ca))]
        plt.figure()
        plt.plot(ca_sim/max(ca_sim),'.-',label="sim")
        plt.plot(ca_exp[16:]/max(ca_exp),'.-',label="exp")
        plt.plot(vec,conc_over_time_Free_Ca/max(conc_over_time_Free_Ca),label="conc time")
        plt.legend()
        plt.show()

        """ Scaling the Experimental data   """

        scale_the_data(
            sim_ring_data_Ca=sim_ring_data_Ca
            ,sim_ring_data_Ann=sim_ring_data_Ann
            ,real_data_Ca=real_data_Ca
            ,real_data_Ann=real_data_Ann
        )

        area_time_Ca_sim = Area_under_graph(
            graph_data=sim_ring_data_Ca
            ,ring_dist=1
        )

        area_time_Ca_exp = Area_under_graph(
            graph_data=real_data_Ca
            ,ring_dist=1
        )

        print(np.shape(sim_ring_data_Ca))
        k = 0
        for t in range(np.shape(sim_ring_data_Ca)[0]):
            if sim_ring_data_Ca[t][0] > 0.00035:
                k = t
                print(t)
        
            
        plt.figure()
        plt.plot(area_time_Ca_sim/max(area_time_Ca_sim),label="sim")
        plt.plot(area_time_Ca_exp[16:]/max(area_time_Ca_exp),label="exp")
        plt.legend()
        plt.show()



    pass
if __name__ == "__main__":
    
    #Annexin_equil_ratio()
    #Finding_the_pump_value()
    Find_pump_via_area()
