import numpy as np
import pandas as pd
# from plotnine import *
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
pd.options.plotting.backend = 'plotly'
# pd.options.plotting.backend = 'matplotlib'
# pio.templates.default = 'none'
pio.templates.default = 'plotly_white'

import colorcet as cc
px.defaults.color_continuous_scale = cc.bjy

# px.defaults.color_continuous_scale = px.colors.diverging.Picnic
# px.defaults.color_continuous_scale = px.colors.sequential.Plasma_r
# px.defaults.color_continuous_scale = px.colors.sequential.Plasma

px.defaults.color_discrete_sequence = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
sky_blue = "#56B4E9"
orange = "#E69F00"
green = "#009E73"

# px.defaults.width = 1300
px.defaults.height = 400
height = 400

heatmap = {'kind' : 'imshow', 'aspect' : 'auto', 'color_continuous_midpoint': 0}

def lineup_plots(plots, specs = None):
    
    fig = make_subplots(cols=len(plots), rows=1, specs = specs)

    for i, plot in enumerate(plots):
        for trace in range(len(plot["data"])):
            fig.append_trace(plot["data"][trace], col=i+1, row=1)

    fig.show()

def plot_patches(patches):
    
    df = patches.reset_index(drop=True)

    mask = np.zeros_like(df.corr(), dtype=bool)
    mask[np.triu_indices_from(mask)] = True

    df.T.plot(height = height, **heatmap).show()
    lineup_plots([df.mean().plot.barh(height = height),
                  df.sum(axis=1).value_counts(ascending=False, normalize=True).plot.bar(height = height)
    ])
    df.corr().plot(height = height, **heatmap).show() #.mask(mask).
    # df.cov().plot(height = height, **heatmap).show()  #.mask(mask).
    
    return patches