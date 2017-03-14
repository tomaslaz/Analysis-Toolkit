# KLMC Analysis Tool Kit
# Author: Tomas Lazauskas, 2015-2017
# www: www.lazauskas.net

A set of scripts to pre/post-process data for/from FHI-aims/GULP/KLMC simulations.

# Requirements 
Python v.2.7.x (x => 9)

Matplotlib v.1.x (x => 5.0)

numpy v.1.x (x => 10.1)

ADUL v.2.0 (For DM_Surface_Energy)

# Scripts

### DA_FHIaims_geometry_stability
(Being developed) Analyses a geometry with FHIaims - pushes a system in the direction of a chosen eigenvector and performs FHIaims single point evaluation.

### DA_Find_Defects
(Being developed) Analyses the change in a system by comparing two xyz files: initial and final configurations.

### DM_Comparison 
(Being developed) Compares structures in terms of their energy ranking between different levels of theory (IPs (GULP) and DFT (FHIaims))

![Comparison example](exampleImages/DM_Comparison.png)

### DM_Convert_Files 
Converts files from one format to another. Works with the most popular formats, such us XYZ, CAR, and GIN.

### DM_Coordination_Bonding 
Analyses systems (xyz format) in terms of avg. bond distance and coordination.

### DM_DOS
Plots DOS (and integrated DOS) graphs [Based on the David Mora Fonz's (UCL) implementation].

![DOS example](exampleImages/DM_DOS.png)

### DM_FHIaims_analysis 
Analyses FHI-aims simulations in terms of runtime, systems' energies etc.

### DM_FHIaims_Spin_Analysis 
Prepares and executes FHIaims spin polarized simulations. 

Reads in structures in xyz format, prepares geometry.in and control.in files by varying the default_initial_moment keyword. The script can also run the prepared simulations one by one. 

### DM_Surface_Energy
Calculates cluster's surface energy. 

Reads in a structure from an xyz file, estimates the area of a cluster using a predefined radius and calculates cluster's surface energy using a provided bulk energy of one atom.

### DM_RDF
Plots radial distribution function.

![RDF example](exampleImages/DM_RDF.png)
### GA_Energy_Evolution 
Plots the energy evolution graph of the n lowest energy structures during a KLMC GA simulation.

![Evolution example](exampleImages/GA_Energy_Evolution.png)

### GA_Energy_Histogram 
Plots energy distribution histogram from the last iteration of a KLMC GA simulation.

![Histogram example](exampleImages/GA_Energy_Histogram.png)

### GA_Family_Tree 
Finds parents, grandparents, etc for a specific GA iteration.

### GA_GM_Iter
Finds the GA iteration number on which the GM has been found.