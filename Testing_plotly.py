import plotly.graph_objects as go
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