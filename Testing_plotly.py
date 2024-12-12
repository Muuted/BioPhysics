import plotly.graph_objects as go
import matplotlib.pyplot as plt
import numpy as np
import cv2
import argparse
import glob
from pathlib import Path
import shutil
import os

def Make_video2(
        output_path
        ,input_path
        ,video_name
        ,fps
        ):
       
    
    # Create a list of all the input image files
    FILES = []
    num = 0
    for frame in os.listdir(input_path):
        num += 1#FILES.append(frame)
    
    FILES = [f"{i}.png" for i in range(num-1)]
    
    # Get the filename from the output path
    filename = video_name
    #print(f'Creating video "{filename}" from images "{FILES}"')

    # Load the first image to get the frame size
    frame = cv2.imread(input_path + FILES[0],cv2.IMREAD_UNCHANGED)
    

    print(f" first file name = {FILES[0]}")
    
    height, width, layers = np.shape(frame)
    
    # Create a VideoWriter object to write the video file
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')

    video = cv2.VideoWriter(filename=filename, fourcc=fourcc
                            , fps= fps
                            , frameSize=(width, height)
                            )

    # Loop through the input images and add them to the video
    for image_path in FILES:
        #print(f'Adding image "{image_path}" to video "{output_path}"... ')
        video.write(cv2.imread( input_path + image_path))
        plt.pause(0.05)

    # Release the VideoWriter and move the output file to the specified location
    cv2.destroyAllWindows()
    video.release() 

    shutil.move(filename, output_path)


def Animate_multi_figure(X_frames,Y_frames_sim,Y_frames_data
        ,figure_title: str
        ,ymax,ymin,xmax,xmin
        ,time_vec
        ,dt
        ,data_hole_open
        ,save_path = ""
        ,fig_name = ""
        ):
    time_shape = np.shape(Y_frames_sim)[0]
    

    fig_dict = {
        "data": []
        ,"layout" : {}
        ,"frames": []
    }

    fig_dict["layout"]["xaxis"] = {"range": [xmin, xmax], "title": "Ring"}
    fig_dict["layout"]["yaxis"] = {"range":[ymin,ymax],"title": figure_title}
    fig_dict["layout"]["hovermode"] = "closest"
    fig_dict["layout"]["updatemenus"] = [
        {
            "buttons": [
                {
                    "args": [None, {"frame": {"duration": 500, "redraw": True},
                                    "fromcurrent": True, "transition": {"duration": 300,
                                                                        "easing": "quadratic-in-out"}}],
                    "label": "Play",
                    "method": "animate"
                },
                {
                    "args": [[None], {"frame": {"duration": 0, "redraw": True},
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
        "active": 235,
        "yanchor": "top",
        "xanchor": "left",
        "currentvalue": {
            "font": {"size": 20},
            "prefix": "frame:",
            "visible": True,
            "xanchor": "right"
        },
        "transition": {"duration": 10, "easing": "linear"},
        "pad": {"b": 10, "t": 50},
        "len": 0.9,
        "x": 0.1,
        "y": 0,
        "steps": []
    }
    tot_time = round(time_vec[len(time_vec)-1]*dt,1)


    # Make the data 

    data_dict_sim = {
            "x": X_frames,
            "y": Y_frames_sim[0],
            "mode": "lines",
            "text": ' Sim',
            "marker": {
                "sizemode": "area",
                "sizeref": 200000
                #,"size": list(dataset_by_year_and_cont["pop"])
            },
            "name": ' sim'
        }
    fig_dict["data"].append(data_dict_sim)

    data_dict_data = {
            "x": X_frames,
            "y": Y_frames_data[data_hole_open],
            "mode": "lines",
            "text": ' Sim',
            "marker": {
                "sizemode": "area",
                "sizeref": 200000
                #,"size": list(dataset_by_year_and_cont["pop"])
            },
            "name": ' data'
        }
    fig_dict["data"].append(data_dict_data)


    for frames in range(time_shape - data_hole_open):
        time_now = round(time_vec[frames]*dt,1)
        
        frame = {
            "data": []
            ,"name": str (f" time = {time_now}s of {tot_time}s")

        }
        
        data_dict_sim = {
            "x": X_frames,
            "y": Y_frames_sim[frames],
            "mode": "lines",
            "text": ' Sim',
            "marker": {
                "sizemode": "area",
                "sizeref": 200000
                #,"size": list(dataset_by_year_and_cont["pop"])
            },
            "name": ' sim'
        }
        frame["data"].append(data_dict_sim)
        
        data_dict_data = {
            "x": X_frames,
            "y": Y_frames_data[frames + data_hole_open],
            "mode": "lines",
            "text": ' Sim',
            "marker": {
                "sizemode": "area",
                "sizeref": 200000
                #,"size": list(dataset_by_year_and_cont["pop"])
            },
            "name": ' sim'
        }
        frame["data"].append(data_dict_data)

        fig_dict["frames"].append(frame)
      
        slider_step = {"args": [
            [frames],
            {"frame": {"duration": 0, "redraw": True},
            "mode": "immediate",
            "transition": {"duration": 0}}
            ],
            "label": frames,
            "method": "animate"
        }
        
        sliders_dict["steps"].append(slider_step)
    
    
    #fig_dict["layout"]["sliders"] = [sliders_dict]



    fig = go.Figure(fig_dict)



    fig.show()

    """
    time_shape = np.shape(Y_frames_sim)[0]
    Animate_figures = go.Figure()

    Animate_figures.add_trace(
        go.Scatter(
            x = X_frames
            ,y = Y_frames_sim[0]
            ,mode = 'lines'
            ,name = ' Sim'            
        )
        ,frames =[
            go.Frame(
                data=[go.Scatter(x=X_frames,y=Y_frames_sim[frame])]
                ,layout=go.Layout(
                    title_text=figure_title 
                    +f"   time = {round(time_vec[frame]*dt,1)}s of {round(time_vec[len(time_vec)-1]*dt,1)}s"
                    )
                ) 
                for frame in range(time_shape)
                ]
    )

    Animate_figures.add_trace(
        go.Scatter(
            x = X_frames
            ,y = Y_frames_data[0]
            ,mode = 'lines'
            ,name =  ' Data'
        )
        ,frames =[
            go.Frame(
                data=[go.Scatter(x=X_frames,y=Y_frames_data[frame])]
                ,layout=go.Layout(
                    title_text=figure_title 
                    +f"   time = {round(time_vec[frame]*dt,1)}s of {round(time_vec[len(time_vec)-1]*dt,1)}s"
                    )
                ) 
                for frame in range(time_shape)
                ]
    )   
    
    Animate_figures.update_layout(
            xaxis = dict(
                range = [xmin,xmax]
                ,autorange=False
            )
            ,yaxis = dict(
                range=[ymin,ymax]
                ,autorange=False
            )
            ,title = dict(
                text = figure_title +f"   time = {round(time_vec[0]*dt,1)}s of {round(time_vec[len(time_vec)-1]*dt,1)}s"
            )
            ,updatemenus=[dict(
                type="buttons"
                ,buttons=[dict(
                    label= "Play"
                    ,method = "animate"
                    ,args=[None]
                )]
            )]
        ,showlegend = True
    )
    Animate_figures.show()
    Animate_figures.write_html(save_path+fig_name)
    """

def Animate_figures(
        X_frames,Y_frames_sim,Y_frames_data
        ,figure_title: str
        ,ymax,ymin,xmax,xmin
        ,time_vec
        ,dt
        ,save_path = ""
        ,fig_name = ""
        ):
    time_shape = np.shape(Y_frames_sim)[0]

    animate_fig = go.Figure(
        data=[
            go.Scatter(
                x = X_frames
                ,y = Y_frames_sim[0]
                )
            ]
        ,layout = go.Layout(
            xaxis = dict(
                range = [xmin,xmax]
                ,autorange=False
            )
            ,yaxis = dict(
                range=[ymin,ymax]
                ,autorange=False
            )
            ,title = dict(
                text = figure_title +f"   time = {round(time_vec[0]*dt,1)}s of {round(time_vec[len(time_vec)-1]*dt,1)}s"
            )
            ,updatemenus=[dict(
                type="buttons"
                ,buttons=[dict(
                    label= "Play"
                    ,method = "animate"
                    ,args=[None]
                )]
            )]
        )

        ,frames =[
            go.Frame(
                data=[go.Scatter(x=X_frames,y=Y_frames_sim[frame])]
                ,layout=go.Layout(
                    title_text=figure_title 
                    +f"   time = {round(time_vec[frame]*dt,1)}s of {round(time_vec[len(time_vec)-1]*dt,1)}s"
                    )
                ) 
                for frame in range(time_shape)
                ]
        )
    

    animate_fig.show()
    animate_fig.write_html(save_path+fig_name)