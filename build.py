'''
This script builds antibodies structures from the dataframe
'''

import os
import sys
os.environ["KMP_DUPLICATE_LIB_OK"] = "True"

from ImmuneBuilder import ABodyBuilder2
import pandas as pd

predictor = ABodyBuilder2()

# Get script path
file_path  = os.path.split(__file__)[0]

# The build script
def build_structures(input_path, output_path):
  
  # Get sequence file
  if not os.path.exists(input_path):
    raise Exception('Unable to find sequence data')

  # Create structure path if not exists
  if not os.path.exists(output_path):
    os.mkdir(output_path)

  # Load the list of sequences
  df = pd.read_csv(input_path, sep='\t')

  # Iterate the sequences
  print('Total number of antibodies: {}'.format(len(df)))
  for index, row in df.iterrows():
    name = row['sequence_id']
    heavy = row['heavy']
    light = row['light']
    print('Processing sequence {}'.format(name))
    ab_path = os.path.join(output_path, '{}.pdb'.format(name))
    if os.path.exists(ab_path):
      continue

    try:
      # Run predictor
      antibody = predictor.predict({
        'H': heavy,
        'L': light
      })

      # Store result
      antibody.save(ab_path)
    except Exception as e:
      print(str(e))

if __name__ == "__main__":
  
  # Prepare input and output paths
  sequence_path = os.path.join(file_path, 'data', 'sh_data.tsv')
  output_path = os.path.join(file_path, 'results', 'structure')
  
  # Build the structures
  build_structures(sequence_path, output_path)