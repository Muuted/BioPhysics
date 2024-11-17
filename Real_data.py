import matplotlib.pyplot as plt
import numpy as np

from Constants import constants
import pandas as pd


def main():
    data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"
    ref_struct_name = "ref_struct__filenum4.txt"



    Real_sim_time = 100
    fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
    fig_folder_path =  fig_save_path + f"simtime={Real_sim_time}\\"
    fig_name = f"Simulation_data_simtime={Real_sim_time}.pkl"

    df = pd.read_pickle(fig_folder_path+fig_name)


    print(df.info())

if __name__ == "__main__":
    main()