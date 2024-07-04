# Reset all
reinitialize everything

# Show contacts function
python
DEBUG=1

def color_h(selection='all'):
        s = str(selection)
        print(s)
        cmd.set_color('color_ile',[0.996,0.062,0.062])
        cmd.set_color('color_phe',[0.996,0.109,0.109])
        cmd.set_color('color_val',[0.992,0.156,0.156])
        cmd.set_color('color_leu',[0.992,0.207,0.207])
        cmd.set_color('color_trp',[0.992,0.254,0.254])
        cmd.set_color('color_met',[0.988,0.301,0.301])
        cmd.set_color('color_ala',[0.988,0.348,0.348])
        cmd.set_color('color_gly',[0.984,0.394,0.394])
        cmd.set_color('color_cys',[0.984,0.445,0.445])
        cmd.set_color('color_tyr',[0.984,0.492,0.492])
        cmd.set_color('color_pro',[0.980,0.539,0.539])
        cmd.set_color('color_thr',[0.980,0.586,0.586])
        cmd.set_color('color_ser',[0.980,0.637,0.637])
        cmd.set_color('color_his',[0.977,0.684,0.684])
        cmd.set_color('color_glu',[0.977,0.730,0.730])
        cmd.set_color('color_asn',[0.973,0.777,0.777])
        cmd.set_color('color_gln',[0.973,0.824,0.824])
        cmd.set_color('color_asp',[0.973,0.875,0.875])
        cmd.set_color('color_lys',[0.899,0.922,0.922])
        cmd.set_color('color_arg',[0.899,0.969,0.969])
        cmd.color("color_ile","("+s+" and resn ile)")
        cmd.color("color_phe","("+s+" and resn phe)")
        cmd.color("color_val","("+s+" and resn val)")
        cmd.color("color_leu","("+s+" and resn leu)")
        cmd.color("color_trp","("+s+" and resn trp)")
        cmd.color("color_met","("+s+" and resn met)")
        cmd.color("color_ala","("+s+" and resn ala)")
        cmd.color("color_gly","("+s+" and resn gly)")
        cmd.color("color_cys","("+s+" and resn cys)")
        cmd.color("color_tyr","("+s+" and resn tyr)")
        cmd.color("color_pro","("+s+" and resn pro)")
        cmd.color("color_thr","("+s+" and resn thr)")
        cmd.color("color_ser","("+s+" and resn ser)")
        cmd.color("color_his","("+s+" and resn his)")
        cmd.color("color_glu","("+s+" and resn glu)")
        cmd.color("color_asn","("+s+" and resn asn)")
        cmd.color("color_gln","("+s+" and resn gln)")
        cmd.color("color_asp","("+s+" and resn asp)")
        cmd.color("color_lys","("+s+" and resn lys)")
        cmd.color("color_arg","("+s+" and resn arg)")
cmd.extend('color_h',color_h)

python end

# Load the reference structure
fetch 7BZ5, WT_RBD, async=0

# Load antibody structures
python

import os

# Retrieve the variable from the command-line arguments
file_path = "./"

# Load JN1
cmd.load(os.path.join(file_path, 'data', '7BZ5_last.pdb'), 'JN1_RBD')

# Align the reference structure
cmd.align('JN1_RBD', 'WT_RBD')

# Remove the reference structure
cmd.remove('WT_RBD')

# Separate RBD from wild-type Ab
cmd.create('WT_ab', 'JN1_RBD and not chain A')
cmd.remove('JN1_RBD and not chain A')

one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C', 'UNK': 'X'}

# Set RBD color
cmd.color('gray80', 'JN1_RBD')
cmd.color('atomic', 'JN1_RBD and not elem C')

# Get all residues for RBD chain A
residue_map = {}
cmd.iterate('JN1_RBD and n. CA', 'residue_map[resv] = resn')

# Set labels for specific residues
cmd.set('float_labels', True)
cmd.set("label_connector", True)   # adds a line
cmd.set("label_connector_mode", 1) # how the line looks label_position
cmd.set("label_position", [-1.5, -1.5, -1.5])  # offset the label
cmd.set("ray_label_connector_flat", True) # cylinders as line? No thanks.
cmd.set("label_outline_color", "white")
cmd.set("label_size", 20)
cmd.set("label_font_id", 13)

RBD_list = [455, 417, 421, 456]
for residue in RBD_list:
  print(residue, residue_map[residue], one_letter[ residue_map[residue] ])
  residue_name = one_letter[ residue_map[residue] ]

  # Show licorice
  cmd.select('RBD_residue', f"JN1_RBD and chain A and resi {residue}")
  cmd.show("sticks", 'RBD_residue')
  cmd.delete('RBD_residue')

# Show original antibody residue Y33
cmd.select('Y33', 'WT_ab and chain H and resi 33')
cmd.show("sticks", 'Y33')
cmd.color("lightpink", "WT_ab")
cmd.color("atomic", "WT_ab and not elem C")

# Hide portion of the original antibody
cmd.select("portion_to_keep", "WT_ab and chain H and resi 27-40")
cmd.remove("WT_ab and not portion_to_keep")
cmd.delete("portion_to_keep")

# Selection of structures
structure_list = ['PLATE5_C08','PLATE3_D09','PLATE3_F08','PLATE6_D09','PLATE6_F12','PLATE2_D05','PLATE5_C10','PLATE3_D05','PLATE3_G06','PLATE4_B02','PLATE2_F01','PLATE4_C06','PLATE4_F06','PLATE4_G03','PLATE5_B04','PLATE5_B06','PLATE5_C05','PLATE5_C06','PLATE5_C09','PLATE5_E02','PLATE5_G09','PLATE6_A12','PLATE6_B05']

# List all files in the directory
pdb_directory = os.path.join(file_path, 'results', 'structure')
if os.path.exists(pdb_directory):
  for pdb_file in os.listdir(pdb_directory):
    if pdb_file.endswith(".pdb"):
      pdb_name = os.path.splitext(pdb_file)[0]
      
      cmd.load(os.path.join(pdb_directory, pdb_file), pdb_name)
      print(pdb_name)

      # Align the loaded structure to the reference structure (AB)
      cmd.align(pdb_name, "WT_ab")

      # Append the loaded object name to the list
      cmd.group("MAD_ab", pdb_name)

      # Select residue 37 of chain L in each structure
      selection = f"{pdb_name}_residue38"
      cmd.select(selection, f"{pdb_name} and chain H and resi 38")
      
      # Show the selected residue as licorice
      cmd.show("sticks", selection)

      # Remove selection
      cmd.delete(selection)

      # This is done for visualization purpose
      selection_to_keep = f"{pdb_name} and chain H and resi 27-40"  # Adjust the residue range as needed
      cmd.select("portion_to_keep", selection_to_keep)
      cmd.remove(f"{pdb_name} and not portion_to_keep")
      cmd.delete("portion_to_keep")

      # Do not show molecules not in list
      if pdb_name not in structure_list:
        cmd.disable(pdb_name)

# Set surface for RBD
cmd.show('surface', 'JN1_RBD')
cmd.set('transparency', 0.3, 'JN1_RBD')

python end

# Set background
bg_color white

# Show cartoons as tube
cartoon tube
set cartoon_tube_radius, 0.3
set cartoon_side_chain_helper, on

# Set representation for MAD_lab antibodies
set transparency_mode, 1
# set stick_transparency, 0.9, MAD_ab
# set cartoon_transparency, 0.7, MAD_ab
color yellow, MAD_ab
color yellow, WT_ab

# Hide hydrogens from the antibody
hide (hydro AND WT_ab) 
hide (hydro AND MAD_ab)

# Dash style
set dash_as_cylinders, 0
set dash_color, black

# Set some styles
set spec_power, 1200
set spec_reflect, 0.20
set ray_opaque_background, 1
set ray_shadow, 0
set ray_texture, 1

set surface_quality, 2
set light_count, 0
set ray_trace_fog,0
set ray_shadows,0
unset depth_cue
bg_color white
set antialias,2
set hash_max, 300
set ray_trace_mode,  1

color_h JN1_RBD

set_view (\
    -0.086982585,   -0.946105182,   -0.311902553,\
    -0.991281986,    0.051119424,    0.121393412,\
    -0.098904766,    0.319746822,   -0.942323029,\
     0.000238499,   -0.000779249,  -38.530986786,\
   -56.881885529,  -36.415988922,    5.768488407,\
   -38.703334808,  115.721618652,   20.000000000 )



