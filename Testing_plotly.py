import plotly.graph_objects as go
import numpy as np
"""

    fig1 = go.Figure()
    fig1.add_trace(
        go.Scatter(
            x=time_vec
            ,y=A_sumfree#[i/max(A_sumfree) for i in A_sumfree]
            ,mode="lines + markers"
            ,name="simulated"
        )
    )
    fig1.update_xaxes(title_text="timesteps")
    fig1.update_yaxes(title_text="[A_free]")
    fig1.update_yaxes(
        exponentformat="SI"
    )
    fig1.show()

    fig2 = go.Figure()
    fig2.add_trace(
        go.Scatter(
            x=time_vec
            ,y=A_f_stability#[i/max(A_f_stability) for i in A_f_stability]
            ,mode="lines + markers"
            ,name="Analytical sol"
        )
    )
    fig2.update_xaxes(title_text="timesteps")
    fig2.update_yaxes(title_text="[A_free]")
    fig2.show()

    fig3 = go.Figure()
    fig3.add_trace(
        go.Scatter(
            x=time_vec
            ,y=A_sumtot
            ,name="total sum"
        )
    )
    fig3.show()
    """



def Animate_figures(
        X_frames,Y_frames
        ,figure_title: str
        ,ymax,ymin,xmax,xmin
        ,time_vec
        ,dt
        ):
    time_shape = np.shape(Y_frames)[0]

    animate_fig = go.Figure(
        data=[
            go.Scatter(
                x = X_frames
                ,y = Y_frames[0]
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
                text = figure_title
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
                data=[go.Scatter(x=X_frames,y=Y_frames[frame])]
                ,layout=go.Layout(
                    title_text=figure_title 
                    +f"   time = {int(time_vec[frame]*dt)}s of  {int(time_vec[len(time_vec)-1]*dt)}s"
                    )
                ) 
                for frame in range(time_shape)
                ]
        )
    
    animate_fig.show()
