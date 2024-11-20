import matplotlib.pyplot as plt
import numpy as np
import scipy.io as spio
import h5py
from Constants import constants
import pandas as pd
from Data_extraction_funcs import *
#from Circle_sim import *
from Circle_funcs import *
     

def compress_img(path, img_name,ysize_goal):
    ref_string = "\\" + "\\"
    for i in range(3):
        if path[len(path)-2:len(path)-1] != ref_string:
            if img_name[0:2] != path[len(path)-2]:
                img_name = "\\" + img_name #just making sure the name is do able
    
    full_img = pd.read_csv(path+img_name,delimiter=',',header=None)
    #full_img = np.genfromtxt(path+img_name,delimiter=',')
    full_img.to_numpy
    Ysize,Xsize = full_img.shape

    compression = int(Ysize/ysize_goal)
    xsize_goal = int(Xsize/compression)

    img_compressed = np.zeros(shape=(ysize_goal,xsize_goal))
    comp_y,comp_x = 0,0
    for y in range(compression,Ysize,compression):
        if y >= Ysize:
                y = Ysize - compression

        for x in range(compression,Xsize,compression):
            if x >= Xsize:
                x = Xsize - compression

            sum = 0
            i = 0
            if y%(int(Ysize)*0.1) == 0:
                 print(f"y={y} of {Ysize}")
            for y_sum in range(y - compression, y + compression):
                    if y_sum >= Ysize-1:
                        break
                    for x_sum in range(x - compression, x + compression):
                        if x_sum >= Xsize-1:
                            break
                        sum += full_img.loc[y_sum][x_sum]
                        i += 1
            img_compressed[comp_y][comp_x] = sum/i
            comp_x += 1
            if comp_x >= np.shape(img_compressed)[1]-1:
                 comp_x = np.shape(img_compressed)[1]-1
        comp_y += 1
        if comp_y >= np.shape(img_compressed)[0]-1:
            comp_y = np.shape(img_compressed)[0]-1

    plt.matshow(full_img)
    plt.title("full image")

    plt.matshow(img_compressed)
    plt.title("compressed image")

    plt.show()

"""
data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"
Calcium_data= "Ca_time_images_filenum4.txt"
"""
"""compress_img(
    path=data_path
    ,img_name=Calcium_data
    ,ysize_goal=80
)"""
"""
ref_struct_name = "ref_struct__filenum4.txt"
c_in,c_out,D_Ca_cyto,T_tot,len_size,dx,dy,k1,k2,c_in_annexin,bound_annexin_start,A_b_init,D_Annexin_cyto,dt,close_time,c_pump,holesize,dR,R,x0,y0,wall_val,inside_val,outside_val,open_val = constants()
py_ref_struct,outside_val,inside_val,wall_vall,xoffset,yoffset = make_ref_structure(
    path=data_path , ref_name=ref_struct_name
)


holesize = 3
open= True
print(f"open_val={open_val} \n outside_val={outside_val} \n wall_val={wall_val} \n inside_val={inside_val}")

ref_structure = init_ref_circle(
            boxlen=len_size
            ,Radius=R,dRadius=dR ,offsets=[x0,y0]
            ,inside_val=inside_val
            ,outside_val=outside_val
            ,wall_val=wall_val
                                        )

for i in range(6):
    open_close_membrane2(
            Grid=ref_structure
            ,offsets=[x0,y0]
            ,holesize=holesize
            ,open_val=open_val
            ,wall_val=wall_val
            ,open_wall_bool=open
                            )
    if i%2 == 0:
         open = False
    if i%2 != 0:
         open = True

    plt.matshow(ref_structure)
    plt.show()"""

N = 2
matrix = np.zeros(shape=(N,N))

df = pd.DataFrame({
     'A': [matrix],
    'B': 10,
    'C': 5
    })

df2 = pd.DataFrame({
     'A':[matrix],
     'D':5
})

df = df.append(df2,ignore_index=True)
#print(df.info())
print(df.info())

matrix2 = np.zeros(shape=(N,N))
for i in range(2):
     for j in range(2):
          matrix2[i][j] = 10



if 'A' in df.columns:
    df = df.drop(columns=['A'])
    print("hello")
else:
     print("not there")

print(df.info())


df2 = pd.DataFrame({
     'A':[matrix2],
     'D':5
})

df = df.append(df2,ignore_index=True)

print(" new ")
print(df.info())

