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




import pandas as pd

file_path = data_path + Calcium_data
df = pd.read_csv(file_path)


#print(df.loc[0][0]) # this is how you get pixel data
#print(df.loc[1][1]) # this is how you get pixel data
Ysize,Xsize = df.shape
sum = 0
i = 0
for y in range(Ysize):
    for x in range(Xsize):
        
        sum += df.loc[0][0] #df.loc[0][1]
        i +=1
avg_val = sum/i
print(f"average val =")
print(avg_val)
exit()

df_yshape, df_xshape = df.shape

def sum_small_circle(x0,y0,xpos,ypos,matrix,circ_radius):

    sum = 0
    lattice_points = 0
    for y in range(ypos-circ_radius,ypos+circ_radius):
        for x in range(xpos-circ_radius,xpos+circ_radius):

            r = np.sqrt((x-x0)**2 + (y-y0)**2)

            if r < circ_radius:
                sum += matrix.loc[y][x]
                lattice_points += 1

    return sum,lattice_points

def average_val(matrix):

    Ysize,Xsize = matrix.shape
    sum = 0
    i = 0
    for y in range(Ysize):
        for x in range(Xsize):
            sum += matrix.loc[y][x]
            i +=1
    
    return sum/i

Ysize,Xsize = df.shape
sum = 0
i = 0
for y in range(Ysize):
    for x in range(Xsize):
        
        sum += df.loc[0][0] #df.loc[0][1]
        i +=1
avg_val = sum/i

#avg_val.type()
#print(f"average val ={avg_val}")

exit()

for y in range(df_yshape):
    for x in range(df_xshape):

            df.loc[y][x]



