
import matplotlib.pyplot as plt
import numpy as np

from Constants import constants
import pandas as pd
from Data_extraction_funcs import *

def main_compare():
    Real_time_steps_data = 235
    Real_sim_time = 120

    fig_save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\Python_simulation_data\\"
    fig_folder_path =  fig_save_path + f"simtime={Real_sim_time}\\"
    fig_name = f"Simulation_data_simtime={Real_sim_time}.pkl"

    df_sim = pd.read_pickle(fig_folder_path+fig_name)

    


if __name__ == "__main__":
    main_compare()