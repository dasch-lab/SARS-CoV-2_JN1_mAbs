import os
import sys
import json
import subprocess

import pandas as pd
import numpy as np

# Get script path
file_path = os.path.split(__file__)[0]

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Create the output path
output_path = os.path.join(file_path, 'results')
if not os.path.exists(output_path):
  os.mkdir(output_path)
  
def calculate_frequency_matrix(strings):
  # Determine the length of the strings
  length = len(strings[0])
  
  # Initialize a dictionary to hold counts of each character at each position
  char_counts = {i: {} for i in range(length)}
  
  # Populate the dictionary with character counts
  for string in strings:
      for position, char in enumerate(string):
          if char in char_counts[position]:
              char_counts[position][char] += 1
          else:
              char_counts[position][char] = 1
  
  # Determine all unique characters across all positions
  all_chars = sorted(set(char for counts in char_counts.values() for char in counts))
  
  # Create a frequency matrix
  frequency_matrix = pd.DataFrame(index=all_chars, columns=range(length)).fillna(0)
  
  for position, counts in char_counts.items():
      for char, count in counts.items():
          frequency_matrix.at[char, position] = count
  
  # Normalize the frequency matrix to convert counts to probabilities
  frequency_matrix = frequency_matrix / len(strings)
  
  # Calculate the information content
  # information_content = -np.sum(frequency_matrix * np.log2(frequency_matrix + 1e-9), axis=0)
  
  # Calculate the information content for each character at each position
  ic_matrix = -frequency_matrix * np.log2(frequency_matrix + 1e-9)
  
  return frequency_matrix, ic_matrix

# Open the JSON file
gisaid_path = os.path.join(file_path, 'data', 'hcov19_lineage_stacked.json')
with open(gisaid_path, 'r') as handle:
  data = json.load(handle)

# Get only the variant frequency data
if 'progressionChartData' not in data:
  raise Exception('Invalid GISAID data')

# We want data from all countries
data = data['progressionChartData']
for i in range(len(data)):
  if data[i]['country'] == 'All':
    data = data[i]['countryData']
    break

# Parse GISAID data
GISAID_list = []
for i in range(len(data)):
  variant = data[i]['vn']
  for j in range(len(data[i]['md'])):
    date = data[i]['md'][j]['ym']
    freq = data[i]['md'][j]['p']
    GISAID_list.append({
      'variant': variant,
      'date': date,
      'freq': freq
    })

# Store frequency
GISAID_df = pd.DataFrame(GISAID_list)

# Fix some issues with names
GISAID_df['variant'] = GISAID_df['variant'].str.replace('BA.2.86.1', 'BA.2.86')

# Store DataFrame as TSV file
GISAID_df.to_csv(os.path.join(file_path, 'results', 'GISAD_frequency.tsv'), sep='\t', index=False)

# Convert to a wide table for supplementary
# Pivot the DataFrame to wide format
GISAID_df_wide = GISAID_df.pivot(index='variant', columns='date', values='freq')
GISAID_df_wide = GISAID_df_wide.reindex(sorted(GISAID_df_wide.columns), axis=1)
GISAID_df_wide.to_csv(os.path.join(file_path, 'results', 'GISAD_frequency_wide.tsv'), sep='\t')

# Load the master table
df_master = pd.read_csv(os.path.join(file_path, 'data', 'master_table.tsv'), sep='\t')

# Select a subset of columns
df_master = df_master[ [
  'IGHV;IGHJ Rearrangement',
  'Neutralization Wuhan IC100 (ng/mL)',
  'Neutralization XBB.1.5  IC100 (ng/mL)',
  'Neutralization EG.5.1.1 IC100 (ng/mL)',
  'Neutralization BA.2.86 IC100 (ng/mL)',
  'Neutralization JN.1 IC100 (ng/mL)'
]]

# Rename columns
df_master = df_master.rename(columns={
  'IGHV;IGHJ Rearrangement': 'germline',
  'Neutralization Wuhan IC100 (ng/mL)': 'WT',
  'Neutralization XBB.1.5  IC100 (ng/mL)': 'XBB.1.5',
  'Neutralization EG.5.1.1 IC100 (ng/mL)': 'EG.5.1.1',
  'Neutralization BA.2.86 IC100 (ng/mL)': 'BA.2.86',
  'Neutralization JN.1 IC100 (ng/mL)': 'JN.1',
})

# Remove empty germline rows
df_master = df_master.dropna(subset=['germline'])

# Convert to long format
df_master = pd.melt(df_master, id_vars=['germline'], var_name='variant', value_name='value')

# Fix some issues with names
df_master['germline'] = df_master['germline'].str.replace('IGHV5-a', 'IGHV5-10-1')

# Store results
df_master.to_csv(os.path.join(file_path, 'results', 'ab_functionality.tsv'), sep='\t', index=False)

# Calculate the web logo and logo matrix
df_sequence = pd.read_csv(os.path.join(file_path, 'data', 'sh_data.tsv'), sep='\t')

# Calculate the frequency matrix and information content
frequency_matrix, information_content = calculate_frequency_matrix(df_sequence['CDRH1'].tolist())
frequency_matrix.to_csv(os.path.join(file_path, 'results', 'CDRH1_sequence_logo.tsv'), sep='\t')

# Create the structure for each antibody in the dataset
print('Building structures')
subprocess.run(['python3', os.path.join(script_dir, 'build.py')], check=True)
# build_structures(os.path.join(file_path, 'data', 'sh_data.tsv'), os.path.join(output_path, 'structure'))