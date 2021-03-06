#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
FHI-aims module.

@author Tomas Lazauskas
"""

import copy
import os
import sys
import shutil
import socket
import subprocess

# import IO
# import System
import IO
import System

import numpy as np

_const_def_value = -9999999999.9
_const_control_in = "control.in"
_const_geometry_in = "geometry.in"

_const_nice_day_line = "Have a nice day."

_const_scf_iter = "Begin self-consistency iteration"

_const_curr_spin = "Current spin moment of the entire structure"
_const_spin_N = "| N = N_up - N_down"
_const_spin_S = "| S                 :"
_const_spin_J = "| J                 :"

_const_eigenvalues = "Writing Kohn-Sham eigenvalues."
_const_evs_up = "Spin-up eigenvalues"
_const_evs_down = "Spin-down eigenvalues"

_const_vbm = "Highest occupied state (VBM)"
_const_cbm = "Lowest unoccupied state (CBM)"
_const_occ_num = "Occupation number"
_const_spin_chan = "Spin channel"

_const_homo_lumo = "Overall HOMO-LUMO gap:"

_const_ZPE = "Cumulative ZPE "
_const_without_six = "without first six eigenmodes"
_const_all_freq = "List of all frequencies found:"

_const_total_energy = "| Total energy of the DFT / Hartree-Fock s.c.f. calculation      :"
_const_no_atoms = "Number of atoms"

_const_total_energy_corrected = "Total energy corrected "

def _readAimsStructure(geometryFile, outputFile, relaxed=True, eigenvalues=False):
  """
  Reads in FHI-aims structure
    
  """
  
  system = None
  
  # getting the number of atoms
  success, error, atomsCnt = IO.countLines(geometryFile)
    
  if not success:
    return success, error, system
  
  else:
    system = System.System(atomsCnt)
    
    # if the system is not relaxed, we need to read the atom positions from the geometry.in file
    if not relaxed:
      success, error = _read_aims_geometry(system)
        
    if not success:
      return success, error, system
    
    success, error = _readAimsOutput(outputFile, system, relaxed=relaxed, eigenvalues=eigenvalues)
  
  return success, error, system

def _read_aims_geometry(system, input_file=_const_geometry_in):
  """
  Reads the system structure from the geometry.in type file.
  
  """
  
  error = ""
  success = True
  
  try:
    fin = open(input_file, "r")
    
  except:
    success = False
    error = __name__ + ": Cannot open: " + input_file
  
  atomsLineCnt = 0
  
  for line in fin:
  
    fields = line.strip().split()
    
    if len(fields) > 0:
    
      sym = fields[4].strip()
      
      # adding specie?
      if sym not in system.specieList:
        system.addSpecie(sym)
        
      specInd = system.specieIndex(sym)
      system.specieCount[specInd] += 1
      system.specie[atomsLineCnt] = specInd
    
      # adding pos
      for j in range(3):
        system.pos[atomsLineCnt*3 + j] = float(fields[j+1])
            
      # adding charge
      system.charge[atomsLineCnt] = 0.0
      
      atomsLineCnt += 1
  
  fin.close()
  
  return success, error

def _readAimsOutput(inputFile, system, relaxed=True, eigenvalues=False):
  """
  Reads in FHI-aims output as a system.
  """
  
  error = ""
  success = True
  have_nice_day = False
  
  readAtoms = False
  readCompleted = False
  readEigenvalues = False
  
  read_evs_up = False
  read_evs_down = False
  
  initial_energy_read = False
  initial_energy = 0.0
  
  eigen_values_array = []
  eigen_values_up_array = []
  eigen_values_down_array = []
  
  atomsLineCnt = 0
  
  noOfAtoms = 0
  energy = 0.0
  version = ''
  noOfcores = 0
  runTime = 0.0
  
  success = IO.checkFile(inputFile)
  
  if not success:
    error = __name__ + ": Cannot locate: " + inputFile

    return success, error

  try:
    fin = open(inputFile, "r")
    
  except:
    success = False
    error = __name__ + ": Cannot open: " + inputFile
    
    return success, error
  
  # initialising spin values
  spin_N = _const_def_value
  spin_S = _const_def_value
  spin_J = _const_def_value
  
  vbm = _const_def_value
  vbm_occ_num = _const_def_value
  vbm_spin_chan = _const_def_value
  
  cbm = _const_def_value
  cbm_occ_num = _const_def_value
  cbm_spin_chan = _const_def_value
  
  homo_lumo_gap = _const_def_value
  
  for line in fin:
    
    fields = line.strip().split()
      
    # reading the final atoms positions
    if (readAtoms and (system is not None)):
      
      if atomsLineCnt >= noOfAtoms:
        readCompleted = True
        readAtoms = False
        
      # ignoring the first line
      if atomsLineCnt > 0:
        #print fields
        
        sym = fields[4].strip()
        
        if sym not in system.specieList:
          system.addSpecie(sym)
        
        specInd = system.specieIndex(sym)
        system.specieCount[specInd] += 1
        system.specie[atomsLineCnt-1] = specInd
        
        for j in range(3):
          system.pos[(atomsLineCnt-1)*3 + j] = float(fields[j+1])
        
        system.charge[atomsLineCnt-1] = 0.0
        
      atomsLineCnt += 1
      
    if ((len(fields) > 1) and (fields[0] == "Version")):
      version = fields[1]
        
    if ((len(fields) > 5) and (' '.join(fields[1:4]) == "Number of atoms")):
      noOfAtoms = int(fields[5])
      
      if (system.NAtoms != noOfAtoms):
        success = False
        error = __name__ + ": the number of atoms does not match the original number of atoms"
        
    if ((len(fields) > 5) and (' '.join(fields[1:4]) == "Total energy uncorrected")):
      energy = float(fields[5])
      
      if not initial_energy_read:
        initial_energy = copy.copy(energy)
        initial_energy_read = True
    
    # reset values
    if _const_scf_iter in line:
      vbm = _const_def_value
      vbm_occ_num = _const_def_value
      vbm_spin_chan = _const_def_value
      
      cbm = _const_def_value
      cbm_occ_num = _const_def_value
      cbm_spin_chan = _const_def_value
    
    if _const_spin_N in line:
      spin_N = float(fields[7])
    
    if _const_spin_S in line:
      spin_S = float(fields[3])
    
    if _const_spin_J in line:
      spin_J = float(fields[3])
    
    if _const_vbm in line:
      vbm = float(fields[5])
      
      # stop reading eigenvalues
      readEigenvalues = False
    
    if _const_cbm in line:
      cbm = float(fields[5])
    
    if _const_homo_lumo in line:
      homo_lumo_gap = float(fields[3])
    
    # Highest occupied state (VBM)
    if vbm != _const_def_value and cbm == _const_def_value:
      if _const_occ_num in line:
        vbm_occ_num = float(fields[3])
      
      if _const_spin_chan in line:
        vbm_spin_chan = float(fields[3])
    
    # Lowest unoccupied state (CBM)
    if vbm != _const_def_value and cbm != _const_def_value:
      if _const_occ_num in line:
        cbm_occ_num = float(fields[3])
      
      if _const_spin_chan in line:
        cbm_spin_chan = float(fields[3])
        
    # reading the eigenvalues
    if (eigenvalues and readEigenvalues and (system is not None)):
      
      if read_evs_up and (_const_evs_down in line):
        read_evs_up = False
      
      if read_evs_down and (_const_curr_spin in line):
        read_evs_down = False
        readEigenvalues = False
      
      if len(fields) == 4:
        if read_evs_up:
          eigen_values_up_array.append(fields[3])
          
        elif read_evs_down:
          eigen_values_down_array.append(fields[3])
        
        else:
          eigen_values_array.append(fields[3])
      
      # reading in spin ups?
      if _const_evs_up in line:
        read_evs_up = True
        
      # reading in spin downs?
      if _const_evs_up in line:
        read_evs_down = True
      
    # Checking whether geometry relaxation was performed
    if relaxed:
      if ((len(fields) > 2) and (' '.join(fields[0:3]) == "Final atomic structure:")):
        readAtoms = True
    
    if ((len(fields) == 4) and (fields[0] == "Using")):
      noOfcores = int(fields[1])
    
    if ((len(fields) > 5) and (' '.join(fields[1:4]) == "Total time :")):
      runTime = float(fields[6])
    
    # Checks whether the have a nice day is in the output file. It indicates that the simulation was successful.
    if _const_nice_day_line in line:
      have_nice_day = True
    
    # Start reading eigenvalues
    if _const_eigenvalues in line:
      readEigenvalues = True
      
      eigen_values_array = []
      
      eigen_values_up_array = []
      eigen_values_down_array = []
    
  fin.close()
  
  # saving the eigenvalues
  if eigenvalues:
    if (len(eigen_values_array) > 0):
      system.eigenvalues = np.array(eigen_values_array, np.float128)
    
    if (len(eigen_values_up_array) > 0) : 
      system.evs_up = np.array(eigen_values_up_array, np.float128)
    
    if (len(eigen_values_down_array) > 0):
      system.evs_down = np.array(eigen_values_down_array, np.float128)
    
  if not have_nice_day:
    system.totalEnergy = 99999999.99
    success = False
    error = __name__ + ": it seems that we are not having a good day"
  
  if not success:
    return success, error
  
  if relaxed and not readCompleted:
    success = False
    error = __name__ + ": data has not been read from: " + inputFile + ". It seems FHIaims encountered an error."
    
  elif system is not None:
    
    system.homo_lumo_gap = homo_lumo_gap

    system.vbm = vbm    
    system.vbm_occ_num =vbm_occ_num
    
    system.cbm = cbm
    system.cbm_occ_num = cbm_occ_num
    
    if spin_N != _const_def_value:
      system.spin_N = spin_N
      system.spin_S = spin_S
      system.spin_J = spin_J
      
      system.vbm_spin_chan = vbm_spin_chan
  
      system.cbm_spin_chan = cbm_spin_chan
    
    system.totalEnergy_initial = initial_energy
    system.totalEnergy = energy
    system.energyDefinition = "FHI-aims_" + version
    
    system.noOfcores = noOfcores
    system.runTime = runTime
    
  else:
    success = False
    error = __name__ + ": data has not been read from: " + inputFile
  
  return success, error

def _readAimsFrequenciesFile(fileName):
  """
  Reads in FHI-aims frequencies file 
  
  """
  
  cZPE = 0.0
  wZPE = 0.0
  
  result = None

  if not checkFile(fileName):    
    return result
  
  try:
    fin = open(fileName, "r")
    
  except:
    return result
    
  sectionStarted = False
  sectionEnded = False
  lineCnt = 0
  
  eigenValues = []
  
  for line in fin:
    line = line.strip()
        
    if _const_ZPE in line:
      cZPE = float(line.split()[4])
    
    if _const_without_six in line:
      wZPE = float(line.split()[6])
      
    if _const_all_freq in line:
      sectionStarted = True
    
    if sectionStarted and not sectionEnded:
      
      if len(line) == 0:
        sectionEnded = True
    
    if sectionStarted and not sectionEnded:
    
      if lineCnt > 1:
        array = line.split()
        
        eigenValues.append(array[1])
      
      lineCnt += 1
    
  fin.close()
  
  result = np.array(eigenValues[6:], np.float64)
  
  return result, cZPE, wZPE

def _readAimsOutputEnergyAtoms(inputFile):
  """
  Reads in FHI-aims output as a system.
  
  """
  
  energy = None

  try:
    fin = open(inputFile, "r")
    
  except:    
    return energy
  
  for line in fin:
    
    fields = line.strip().split()
  
    if _const_total_energy in line:
      energy = float(fields[11])
      
    if ((len(fields) > 5) and (' '.join(fields[1:4]) == _const_no_atoms)):
      noOfAtoms = int(fields[5])
      
  fin.close()
  
  return energy, noOfAtoms

def getRelaxStepsEnergies(outputFile):
  """
  Counts the number of relaxation steps
  
  """
  
  relax_steps = []
  error = None
  success = IO.checkFile(outputFile)
  
  if not success:
    error = __name__ + ": Cannot locate: " + outputFile
    return success, error, relax_steps

  try:
    fout = open(outputFile, "r")
    
  except:
    success = False
    error = __name__ + ": Cannot open: " + outputFile
    return success, error, relax_steps
  
  for line in fout:
    # checks if the run was successful
    if _const_nice_day_line in line:
      success = True
    
    # reads in energy
    if _const_total_energy_corrected in line:
      energy = np.float128((line.split(":")[1]).split("eV")[0].strip())
      
      relax_steps.append(energy)

  fout.close()
  
  return success, error, relax_steps
