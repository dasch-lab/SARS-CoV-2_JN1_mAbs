import os
import math

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import PercentFormatter

# Get script path
file_path  = os.path.split(__file__)[0]

# Get the root path
root_path = os.path.abspath(os.path.join(file_path, '..'))

# Make sure the result path exists
result_path = os.path.join(root_path, 'results')
if not os.path.exists(result_path):
  os.path.mkdir(result_path)

# Load the distance frequency file
frequency_path = os.path.join(root_path, 'data', 'MD_distances.tsv')
df_frequency = pd.read_csv(frequency_path, sep='\t')

# Get the number of plots
plot_list = df_frequency['distance'].unique().tolist()

# Prepare the plot
num_cols = 2
num_rows = math.ceil(len(plot_list)/num_cols)
fig, axis = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(4*num_cols, 4*num_rows))
axis = axis.ravel()

# Define the colors
color_map = {
  'WT': {
    'color': 'skyblue',
    'edgecolor': 'white'
  },
  'JN1': {
    'color': 'red',
    'edgecolor': 'white'
  }
}

# Plot the data
for idx, distance in enumerate(plot_list):
    
  # Get the title from the distance name
  title = 'Position {}'.format(distance[1:4])
  
  # Plot the distribution of the 'Value' column
  axis[idx].set_xlabel('Distance (Å)')
  axis[idx].set_ylabel('Frequency')
  axis[idx].set_title(title)
  axis[idx].margins(x=0)
  axis[idx].yaxis.set_major_formatter(PercentFormatter(100))
  
  # Extract unique types
  df_distance = df_frequency[ df_frequency['distance'] == distance]
  
  # Create a bar plot for each type
  types = df_distance['type'].unique()
  
  for typ in types:
    subset = df_distance[df_distance['type'] == typ]
    
    # Make sure last rounded number is plotted
    x_values = subset['bin_start'].tolist()
    x_values.append(math.ceil(max(x_values)))
    y_values = subset['frequency'].tolist()
    y_values.append(0)
    bin_width = subset['bin_end']-subset['bin_start']
    bin_width = bin_width.tolist()
    bin_width.append(x_values[-1]-x_values[-2])

    axis[idx].bar(
      x_values,
      y_values,
      width=bin_width,
      color=color_map[typ]['color'],
      edgecolor=color_map[typ]['edgecolor']
    )
    
# Plot result
fig.tight_layout()
frequency_path = os.path.join(result_path, 'MD_distance.png')
plt.savefig(frequency_path, dpi = 600)
plt.close(fig)