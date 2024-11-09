import matplotlib.pyplot as plt
import numpy as np
import scipy.io as spio
import h5py
from Constants import constants
import pandas as pd

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
    
    
    
def avg_ring_sum(grid,x0,y0,radius):
    sum = 0
    i = 0
    for y in range(y0-radius , y0 + radius):
        for x in range(x0 - radius , x0 + radius):
            r = np.sqrt( (x-x0)**2 + (y-y0)**2)
            if r < radius:
                sum += grid.loc[y][x]
                i += 1
    return sum/i

data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"

Calcium_data= "Ca_time_images_filenum4.txt"
Annexin_data = "Annexin_time_images_filenum4.dat"


import pandas as pd

ca_data = pd.read_csv(data_path+Calcium_data)

plt.matshow(ca_data)
print(ca_data.shape)

Y0,X0 = ca_data.shape

ref_struct = np.zeros(shape=(Y0,X0))

print("np sum start")
avg_val =np.sum( np.sum(ca_data))/(X0*Y0)
print("np sum done")
"""
sum = 0
i = 0

for y in range(Y0-1):
    for x in range(X0-1):
        sum += ca_data.loc[y][x]
        i += 1

        if i%(a)==0:
            print(f"i={int( i*100/(Y0*X0))} %  - avg sum loop")

avg_val = sum/i
"""
a= int((Y0*X0)*0.01)
i = 0
R = 30
for y in range(R+1,Y0-R-1):
    for x in range(R+1,X0-R-1):
        i += 1
        ring_avg_val = avg_ring_sum(
            grid=ca_data
            ,x0=x
            ,y0=y
            ,radius=R
        )
        if ring_avg_val > avg_val:
            ref_struct[y][x] = 20

        if i%(a)==0:
            print(f"i={ i*100/(Y0*X0)} %  - ref matrix loop")

plt.matshow(ref_struct)
plt.show()
#print(f"the sum/i=")
#print(sum/i)



def compress_img(path,img_name):

    data = pd.read_csv(path + img_name)


if __name__=="__main__":
    name = "Adam"
    print(name)

    