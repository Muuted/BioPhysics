import matplotlib.pyplot as plt
import numpy as np
import scipy.io as spio
import h5py
from Constants import constants
import pandas as pd
from Data_extraction_funcs import *
#from Circle_sim import *
from Circle_funcs import *
     

A = [1,2,3,4]
print(A[0:2])
print(A[2:6])


exit()
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

import plotly.graph_objects as go

import pandas as pd

url = "https://raw.githubusercontent.com/plotly/datasets/master/gapminderDataFiveYear.csv"
dataset = pd.read_csv(url)

years = ["1952", "1962", "1967", "1972", "1977", "1982", "1987", "1992", "1997", "2002",
         "2007"]

# make list of continents
continents = []
for continent in dataset["continent"]:
    if continent not in continents:
        continents.append(continent)
# make figure
fig_dict = {
    "data": [],
    "layout": {},
    "frames": []
}

# fill in most of layout
fig_dict["layout"]["xaxis"] = {"range": [30, 85], "title": "Life Expectancy"}
fig_dict["layout"]["yaxis"] = {"title": "GDP per Capita", "type": "log"}
fig_dict["layout"]["hovermode"] = "closest"
fig_dict["layout"]["updatemenus"] = [
    {
        "buttons": [
            {
                "args": [None, {"frame": {"duration": 500, "redraw": False},
                                "fromcurrent": True, "transition": {"duration": 300,
                                                                    "easing": "quadratic-in-out"}}],
                "label": "Play",
                "method": "animate"
            },
            {
                "args": [[None], {"frame": {"duration": 0, "redraw": False},
                                  "mode": "immediate",
                                  "transition": {"duration": 0}}],
                "label": "Pause",
                "method": "animate"
            }
        ],
        "direction": "left",
        "pad": {"r": 10, "t": 87},
        "showactive": False,
        "type": "buttons",
        "x": 0.1,
        "xanchor": "right",
        "y": 0,
        "yanchor": "top"
    }
]

sliders_dict = {
    "active": 0,
    "yanchor": "top",
    "xanchor": "left",
    "currentvalue": {
        "font": {"size": 20},
        "prefix": "Year:",
        "visible": True,
        "xanchor": "right"
    },
    "transition": {"duration": 300, "easing": "cubic-in-out"},
    "pad": {"b": 10, "t": 50},
    "len": 0.9,
    "x": 0.1,
    "y": 0,
    "steps": []
}

# make data
year = 1952

print(f"continents shape = {np.shape(continents)}")
print(f"years shape = {np.shape(years)}")

for continent in continents:
    dataset_by_year = dataset[dataset["year"] == year]
    dataset_by_year_and_cont = dataset_by_year[
        dataset_by_year["continent"] == continent]

    data_dict = {
        "x": list(dataset_by_year_and_cont["lifeExp"]),
        "y": list(dataset_by_year_and_cont["gdpPercap"]),
        "mode": "markers",
        "text": list(dataset_by_year_and_cont["country"]),
        "marker": {
            "sizemode": "area",
            "sizeref": 200000,
            "size": list(dataset_by_year_and_cont["pop"])
        },
        "name": continent
    }
    fig_dict["data"].append(data_dict)

# make frames
for year in years:
    frame = {"data": [], "name": str(year)}
    for continent in continents:
        dataset_by_year = dataset[dataset["year"] == int(year)]
        dataset_by_year_and_cont = dataset_by_year[
            dataset_by_year["continent"] == continent]

        data_dict = {
            "x": list(dataset_by_year_and_cont["lifeExp"]),
            "y": list(dataset_by_year_and_cont["gdpPercap"]),
            "mode": "markers",
            "text": list(dataset_by_year_and_cont["country"]),
            "marker": {
                "sizemode": "area",
                "sizeref": 200000,
                "size": list(dataset_by_year_and_cont["pop"])
            },
            "name": continent
        }
        frame["data"].append(data_dict)

    fig_dict["frames"].append(frame)
    slider_step = {"args": [
        [year],
        {"frame": {"duration": 300, "redraw": False},
         "mode": "immediate",
         "transition": {"duration": 300}}
    ],
        "label": year,
        "method": "animate"}
    sliders_dict["steps"].append(slider_step)


fig_dict["layout"]["sliders"] = [sliders_dict]

fig = go.Figure(fig_dict)

fig.show()
exit()
data_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\ISA Biophys\\data eksperimenter\\20191203-Calcium-sensors-ANXA-RFP for Python\\"

Ca_data = "Ring_sum_data_Ca_filenum4.txt"
Annexin_data = "Ring_sum_data_Annexins_filenum4.txt"

real_data_Ca = pd.read_csv(data_path + Ca_data)
real_data_Ann= pd.read_csv(data_path + Annexin_data)
print(real_data_Ca.loc[0][0])


print(f"the shape of ca data = {np.shape(real_data_Ca)}")




vec = np.linspace(0,len(real_data_Ca.loc[0]),len(real_data_Ca.loc[0]))
plt.figure()
plt.plot(vec,real_data_Ca.loc[0])
plt.show()