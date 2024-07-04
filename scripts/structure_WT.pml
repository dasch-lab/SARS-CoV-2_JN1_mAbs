
# Reset all
reinitialize everything

# Fetch structure
fetch 7BZ5, WT_RBD, async=0

# Cleanup
remove resn hoh
remove resn nag

# Separate antibody from RBD
create WT_ab, WT_RBD and not chain A
remove WT_RBD and not chain A

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

# Load antibody structures
python

import os

# Retrieve the variable from the command-line arguments
file_path = r"./"

#file_path = r"C:\Users\eandreano\Desktop\ERC PostDoc EA\My Articles\SARS-CoV-2_Antibody response JN.1"

one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C', 'UNK': 'X'}

# Set RBD color
cmd.color('gray80', 'WT_RBD')
cmd.color('atomic', 'WT_RBD and not elem C')

# Get all residues for RBD chain A
residue_map = {}
cmd.iterate('WT_RBD and n. CA', 'residue_map[resv] = resn')

# Set labels for specific residues
cmd.set('float_labels', True)
cmd.set("label_connector_color", 0xff8c69) # connect or bounding box
cmd.set("label_connector", True)   # adds a line
cmd.set("label_connector_mode", 1) # how the line looks label_position
cmd.set("label_position", [-1.5, -1.5, -1.5])  # offset the label
cmd.set("ray_label_connector_flat", True) # cylinders as line? No thanks.
cmd.set("label_outline_color", "white")
cmd.set("label_size", 20)
cmd.set("label_font_id", 13)

RBD_list = [455, 417, 421, 456]
for residue in RBD_list:
  residue_name = one_letter[ residue_map[residue] ]

  # Add label
  #cmd.label(f"WT_RBD and chain A and resi {residue} and name CB", f"'{residue_name}{residue}'")  # Replace 'Y33-H' with the label you want to display
  
  # Show licorice
  cmd.select('RBD_residue', f"WT_RBD and chain A and resi {residue}")
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

# Create a copy of the WT_RBD to show whole residues
copy_selection = ' or '.join([ f"resi {residue}" for residue in RBD_list ])
cmd.create('WT_wire', f"WT_RBD and ({copy_selection})")
cmd.show_as('sticks', 'WT_wire')

# Set surface for RBD
cmd.show('surface', 'WT_RBD')
cmd.set('transparency', 0.1, 'WT_RBD')

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
hide (hydro AND MAD_ab) 
hide (hydro AND WT_ab) 

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

# set surface_quality, 2
# set light_count, 5
# set ambient_occlusion_mode, 1
# set ambient_occlusion_scale, 50
# set ambient, 0.40

color_h WT_RBD

set_view (\
    -0.086982585,   -0.946105182,   -0.311902553,\
    -0.991281986,    0.051119424,    0.121393412,\
    -0.098904766,    0.319746822,   -0.942323029,\
     0.000238499,   -0.000779249,  -38.530986786,\
   -56.881885529,  -36.415988922,    5.768488407,\
   -38.703334808,  115.721618652,   20.000000000 )

# Bottom camera
#set_view (\
#    -0.071750842,   -0.749409556,    0.658181131,\
#    -0.961691737,    0.227019876,    0.153644398,\
#    -0.264564604,   -0.621954262,   -0.737000525,\
#     0.000238499,   -0.000779249,  -38.530986786,\
#   -56.881885529,  -36.415988922,    5.768488407,\
#   -38.703334808,  115.721618652,   20.000000000 )

# Ray trace 
python
#cmd.png(os.path.join(file_path, 'RBD_WT_7BZ5_bottom.png'), width=2000, height=2000, ray=1, quiet=1)
python end

# Set side view
#set_view (\
#    -0.210849613,   -0.969662547,   -0.123530388,\
#    -0.931610227,    0.161071643,    0.325807780,\
#    -0.296028763,    0.183780730,   -0.937328160,\
#     0.000238499,   -0.000779249,  -38.530986786,\
#   -56.881885529,  -36.415988922,    5.768488407,\
#   -38.703334808,  115.721618652,   20.000000000 )

# Ray trace 
python
#cmd.png(os.path.join(file_path, 'RBD_WT_7BZ5_side.png'), width=2000, height=2000, ray=1, quiet=1)
python end

# Set the camera view
# set_view (\
#     -0.186612174,   -0.982069969,   -0.026138239,\
#     -0.901905477,    0.160710946,    0.400917202,\
#     -0.389531583,    0.098388724,   -0.915739119,\
#      0.000238499,   -0.000779249,  -39.411396027,\
#    -56.881885529,  -36.415988922,    5.768488407,\
#    -37.822925568,  116.602020264,   20.000000000 )
# set_view (\
#     -0.293418199,   -0.845420659,    0.446259111,\
#      0.322523743,    0.351898253,    0.878716886,\
#     -0.899933457,    0.401761591,    0.169416189,\
#      0.000238499,   -0.000779249,  -39.411396027,\
#    -56.881885529,  -36.415988922,    5.768488407,\
#    -37.822925568,  116.602020264,   20.000000000 )
# set_view (\
#     -0.253469676,   -0.799344897,    0.544774830,\
#      0.289091438,    0.474839061,    0.831234634,\
#     -0.923133612,    0.368182957,    0.110727467,\
#      0.000238499,   -0.000779249,  -39.411396027,\
#    -56.881885529,  -36.415988922,    5.768488407,\
#    -37.822925568,  116.602020264,   20.000000000 )
# set_view (\
#     -0.168432310,   -0.952280641,   -0.254503906,\
#     -0.898018181,    0.041785635,    0.437967926,\
#     -0.406435728,    0.302319080,   -0.862211347,\
#      0.000238499,   -0.000779249,  -39.411396027,\
#    -56.881885529,  -36.415988922,    5.768488407,\
#    -37.822925568,  116.602020264,   20.000000000 )

