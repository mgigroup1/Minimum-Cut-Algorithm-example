#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 13:24:42 2023

@author: wenchingmei
"""
# Replace with the path to your Excel file



import pandas as pd
import plotly.graph_objects as go
import plotly.offline as pyo

# Read the nodes and values
# Replace with the path to your Excel file
node_names = pd.read_csv('/Users/wenchingmei/Desktop/bio_code/DMMM-main/sankey plot/sankey plot/node_names.csv').values.flatten().tolist()
node_values = pd.read_csv('/Users/wenchingmei/Desktop/bio_code/DMMM-main/sankey plot/sankey plot/node_values.csv', header=None).values.flatten().tolist()

# Define empty lists for source, target, and values of the links
source = []
target = []
values = []

link_colors = [colors[s % len(colors)] for s in source]
# Fill the lists using the node values
for idx, value in enumerate(node_values):
    source.append(idx)
    target.append(len(node_names) + 1)  # this is a fake node at the end
    values.append(value)

# Define a Sankey plot
fig = go.Figure(go.Sankey(
    node=dict(
        pad=30, 
        thickness=10, 
        line=dict(color="black", width=1), 
        label=node_names + ["Fake Node"]  # add a fake node label at the end
    ),
    link=dict(
        source=source,
        target=target,
        value=values,
        color=link_colors 
    )
))

# Remove the fake node by setting its color to white (or the background color)
fig.update_traces(node=dict(color=['rgba(255,255,255,0)' if i == len(node_names) else 'blue' for i in range(len(node_names) + 1)]))


# Define a list of softer colors using RGBA
colors = [
    'rgba(0, 0, 255, 0.2)',     # blue
    'rgba(0, 128, 0, 0.2)',    # green
    'rgba(255, 0, 0, 0.2)',    # red
    'rgba(128, 0, 128, 0.2)',  # purple
    'rgba(255, 255, 0, 0.2)',  # yellow
    'rgba(255, 165, 0, 0.2)',  # orange
    'rgba(255, 192, 203, 0.2)',# pink
    'rgba(0, 255, 255, 0.2)',  # cyan
    'rgba(255, 0, 255, 0.2)',  # magenta
    'rgba(139, 69, 19, 0.2)'   # brown
]

y_positions = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1]


# Assign colors based on the source node
link_colors = [colors[s % len(colors)] for s in source]

fig = go.Figure(go.Sankey(
    node=dict(
        pad=30,  
        thickness=10,  
        line=dict(color="black", width=1),
        label=node_names,
        y=y_positions
    ),
    textfont=dict(size=14),  # Adjusted font size and set the family to Arial Bold
    link=dict(
        source=source,
        target=target,
        value=values,
        color=link_colors  # Use the custom colors
    )
))

fig.show()

pyo.plot(fig)


##
##
##
# flux dependent graph
# Read the data
file_path = '/Users/wenchingmei/Desktop/min graph based FBA_paper draft/si/DATA SET S1.xlsx'  # Replace with the path to your Excel file
edges = pd.read_excel(file_path, sheet_name= ' E. coli _flux dependent graph')
 # Replace with the path to your Excel file
node_names = pd.read_csv('/Users/wenchingmei/Desktop/bio_code/DMMM-main/sankey plot/sankey plot/node_names.csv').values.flatten().tolist()


# Assuming source_Gluc and sink_Gluc are 1-indexed, we subtract 1 to make them 0-indexed
source = (edges['source_Gluc'] - 1).tolist() 
target = (edges['sink_Gluc'] - 1).tolist()
values = edges['weight_Gluc'].tolist()


# Define a list of softer colors using RGBA
colors = [
    'rgba(0, 0, 255, 0.2)',     # blue
    'rgba(0, 128, 0, 0.2)',    # green
    'rgba(255, 0, 0, 0.2)',    # red
    'rgba(128, 0, 128, 0.2)',  # purple
    'rgba(255, 255, 0, 0.2)',  # yellow
    'rgba(255, 165, 0, 0.2)',  # orange
    'rgba(255, 192, 203, 0.2)',# pink
    'rgba(0, 255, 255, 0.2)',  # cyan
    'rgba(255, 0, 255, 0.2)',  # magenta
    'rgba(139, 69, 19, 0.2)'   # brown
]

y_positions = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1]


# Assign colors based on the source node
link_colors = [colors[s % len(colors)] for s in source]

fig = go.Figure(go.Sankey(
    node=dict(
        pad=30,  
        thickness=10,  
        line=dict(color="black", width=1),
        label=node_names,
        y=y_positions
    ),
    textfont=dict(size=14),  # Adjusted font size and set the family to Arial Bold
    link=dict(
        source=source,
        target=target,
        value=values,
        color=link_colors  # Use the custom colors
    )
))

fig.show()
pyo.plot(fig)

##
##
##
# flux dependent graph with min-cut

import pandas as pd
import plotly.graph_objects as go

# Read the data
file_path = '/Users/wenchingmei/Desktop/min graph based FBA_paper draft/si/DATA SET S1.xlsx'  # Replace with the path to your Excel file
edges = pd.read_excel(file_path, sheet_name= ' E. coli_mc pathfinding')
edges = edges.dropna(subset=['source_Gluc', 'sink_Gluc', 'weight_Gluc'])
 # Replace with the path to your Excel file
node_names = pd.read_csv('/Users/wenchingmei/Desktop/bio_code/DMMM-main/sankey plot/sankey plot/node_names.csv').values.flatten().tolist()

# Assuming source_Gluc and sink_Gluc are 1-indexed, we subtract 1 to make them 0-indexed
source = [(int(s) - 1) for s in edges['source_Gluc'].tolist()]
target = [(int(s) - 1) for s in edges['sink_Gluc'].tolist()]
values = edges['weight_Gluc'].tolist()





import plotly.offline as pyo
# Define a list of softer colors using RGBA
colors = [
    'rgba(0, 0, 255, 0.2)',     # blue
    'rgba(0, 128, 0, 0.2)',    # green
    'rgba(255, 0, 0, 0.2)',    # red
    'rgba(128, 0, 128, 0.2)',  # purple
    'rgba(255, 255, 0, 0.2)',  # yellow
    'rgba(255, 165, 0, 0.2)',  # orange
    'rgba(255, 192, 203, 0.2)',# pink
    'rgba(0, 255, 255, 0.2)',  # cyan
    'rgba(255, 0, 255, 0.2)',  # magenta
    'rgba(139, 69, 19, 0.2)'   # brown
]

y_positions = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1]


# Assign colors based on the source node
link_colors = [colors[s % len(colors)] for s in source]

fig = go.Figure(go.Sankey(
    node=dict(
        pad=30,  
        thickness=10,  
        line=dict(color="black", width=1),
        label=node_names,
        y=y_positions
    ),
    textfont=dict(size=14),  # Adjusted font size and set the family to Arial Bold
    link=dict(
        source=source,
        target=target,
        value=values,
        color=link_colors  # Use the custom colors
    )
))

fig.show()
pyo.plot(fig)
