# KLMC Analysis Tool Kit
# Author: Tomas Lazauskas, 2015-2016
# www: www.lazauskas.net

A set of scripts to pre/post-process data for/from KLMC simulations.

# Requirements 
Python v.2.7.x (x => 9)

Matplotlib v.1.x (x => 5.0)

numpy v.1.x (x => 10.1)

# Scripts
DM_Comparison - compare unique structures in terms of their energy ranking between different levels of theory (IPs (GULP) and DFT (FHIaims)) (Being developed)

DM_Convert_Files - converts files from one format to another. Works with the most popular formats, such us XYZ, CAR, and GIN.

DM_Coordination_Bonding - analyzes systems (xyz format) in terms of avg. bond distance and coordination.

DM_DOS - plots DOS (and integrated DOS) graphs [Based on the David Mora Fonz's (UCL) implementation].

![ScreenShot](exampleImages/DM_DOS.png)

DM_FHI_aims_analysis - analyzes FHI-aims simulations in terms of runtime, systems' energies etc.

DM_RDF - plots radial distribution function.

GA_Energy_Evolution - plots the energy evolution graph of the n lowest energy structures during a KLMC GA simulation.

GA_Energy_Histogram - plots energy distribution histogram from the last iteration of a KLMC GA simulation.

GA_Family_Tree - finds parents, grandparents, etc for a specific GA iteration.

GA_GM_Iter - finds the GA iteration number on which the GM has been found.