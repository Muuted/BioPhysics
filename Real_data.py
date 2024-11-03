import matplotlib.pyplot as plt
import numpy as np
import scipy.io as spio
import h5py
from Constants import constants


def extract_data_mean_intensity(data_path,datafilename,filenum:str):
    
    file = data_path + datafilename + "\\" + filenum + ".fig" 
    #print(spio.whosmat(file2))
    Timeseries2 = spio.loadmat(file)
    Mean_intensity_T2 = Timeseries2["__function_workspace__"]

    plt.figure()
    plt.plot(Mean_intensity_T2[0])
    plt.show()

def extract_data_Radial_mena_intensity(data_path,datafilename):

    file = data_path + datafilename
    #print(spio.whosmat(file))
    data = h5py.File(file,mode="r")
    keywords = ['#refs#', 'Frame', 'ICAGFP'
                , 'IMaskCAGFP', 'IMaskMeanCAGFP'
                , 'MeanRadialIntensitiesCAGFP'
                , 'Radius', 'Time'
                ]
    
    #with h5py.File(file,mode="r") as f:
     #   print(f.keys())

    dd = data[keywords[0]]
    print(dd)
    
    
    

data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"

Calcium_data= "Ca_time_images_filenum4.dat"
Annexin_data = "Annexin_time_images_filenum4.dat"
