# Supporting data for the paper `JN.1 global spread unleashed by evasion of IGHV3-53/3-66 B cell germlines`

# Supporting data for COVID-19 superhybrid paper

This repository collects the scripts necessary for analysing the JN1 dataset and generating the images for the manuscript. 

## Content ##
```
. 
├── cleanup.py  # Script for generating cleaned-up data
├── scripts     # Scripts generating images
└── data
    ├── JN1_AF2_relaxed.pdb    # JN1 RBD predicted with AlphaFold2
    ├── 7BZ5_last.pdb    # Last frame from MD simulation of 7BZ5 with JN1
    ├── sh_data.tsv    # IGHV3-53/3-66 antibody sequences
    ├── hcov19_lineage_stacked.json    # GISAID lineage statistics
    ├── MD_distances.tsv    # Distances measured between residues in MD simulation 
    └── master_table.tsv               # Antibody master table
```

## Set-up

    conda create -p ./env pip ipython pandas

# Usage #
To generate the data:

    python3 ./cleanup.py

To generate the variant plot:

    Rscript ./scripts/variant.R --frequency ./results/GISAD_frequency.tsv -o ./results/variant.pdf --functionality ./results/ab_functionality.tsv

to generate the MD simulation distance plot:

    python3 ./scripts/distance.py

Generate the MD simulation images, please run the following pymol scripts: 

```console
pymol ./scripts/structure_WT.pml # Structure of Wuhan RBD with b38 Ab
pymol ./scripts/structure_WT_consensus.pml # Structure of Wuhan RBD with IGHV3-53/3-66 germlines
pymol ./scripts/structure_JN1.pml # Structure of JN1 RBD with b38 Ab
pymol ./scripts/structure_JN1_consensus.pml # Structure of JN1 RBD with IGHV3-53/3-66 germlines
```