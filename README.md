# Analyses of the phylogenetic distribution of Lap proteins

This project investigates the phylogenetic distribution of the LapA-like adhesins (which include MapA) and the LapD and LapG regulatory proteins. It uses annotations from the NCBI conserved domain database (CDD) to find LapD and LapG homologs and then uses characteristics common to characterized LapA homologs to identify LapA-like proteins in organisms that encode LapD and LapG. For more information see the corresponding study:

> MapA, a second large RTX adhesin conserved across the Pseudomonads, contributes to biofilm formation by Pseudomonas fluorescens
Alan J. Collins, Alexander B. Pastora, T. Jarrod Smith, George A. O'Toole
Journal of Bacteriology Jul 2020, JB.00277-20; DOI: 10.1128/JB.00277-20

Included in this repository is all code and raw data that were used to identify the LapA-like proteins in bacterial genomes present in the RefSeq database. In addition, a .csv file containing all ideintified LapA-like proteins is included in the root directory above [and here.](https://github.com/GeiselBiofilm/Collins-MapA/blob/master/All_lapg_targets.csv)

## Repository structure

The ['Identify_LapA-like_proteins'](Identify_LapA-like_proteins) directory contains R code to identify organisms that encode LapD, LapG, and LapA homologs, as well as example input files that were used in the published analysis.

The ['LapA-like_Alignments'](LapA-like_Alignments) directory contains python 3 code to process outputs generated by the analysis in the ['Identify_LapA-like_proteins'](Identify_LapA-like_proteins) directory. This analysis involves the alignment of all identified LapA homologs with LapA and MapA and the production of Supplemental files 2 and 3 from the manuscript.

The ['Plotting'](Plotting) directory contains all code used to generate figures 10 and 11 in the manuscript.

The ['Supplemental_File_Alignments']('Supplemental_File_Alignments') directory contains supplemental files referred to in the main text of the manuscript.

In the relevant directories, the data that were used for our analyses are provided in the Data folder. These include the list of LapD and LapG encoding organisms, the list of bacterial genomes present in the GenBank database when our analysis was conducted, and a Newick format encoding of the tree presented in Figure 11 of the manuscript. All other input files can be generated using the scripts provided here.
