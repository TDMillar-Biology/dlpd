'''
Docstring for visualization.renderers
hold matplotlib and plotly renderers working on primitive line objects from
svmu.models.line:PrimitiveLine
'''

def render_matplotlib(primitives, xlabel, ylabel, ax=None):
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(6,6))
    else:
        fig = ax.figure

    for p in primitives:
        ax.plot(p.x, p.y, color=p.color)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return fig, ax

def render_plotly(primitives, title):
    import plotly.graph_objects as go
    traces = []

    for p in primitives:
        traces.append(go.Scatter(
            x=p.x,
            y=p.y,
            mode='lines',
            line=dict(color=p.color),
            hovertext=p.hover_text,
            hoverinfo='text' if p.hover_text else 'skip'
        ))

    return go.Figure(data=traces)