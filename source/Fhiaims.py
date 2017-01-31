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

import IO
import System

_const_def_value = -9999999999.9
_const_control_in = "control.in"
_const_geometry_in = "geometry.in"

_const_nice_day_line = "Have a nice day."

_const_scf_iter = "Begin self-consistency iteration"

_const_curr_spin = "Current spin moment of the entire structure"
_const_spin_N = "| N = N_up - N_down"
_const_spin_S = "| S                 :"
_const_spin_J = "| J                 :"

_const_vbm = "Highest occupied state (VBM)"
_const_cbm = "Lowest unoccupied state (CBM)"
_const_occ_num = "Occupation number"
_const_spin_chan = "Spin channel"

_const_homo_lumo = "Overall HOMO-LUMO gap:"

def _readAimsStructure(geometryFile, outputFile, relaxed=True):
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
    
    success, error = _readAimsOutput(outputFile, system, relaxed=relaxed)
  
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
      system.specie[atomsLineCnt-1] = specInd
      
      # adding pos
      for j in range(3):
        system.pos[(atomsLineCnt-1)*3 + j] = float(fields[j+1])
      
      # adding charge
      system.charge[atomsLineCnt-1] = 0.0
      
      atomsLineCnt += 1
  
  fin.close()
  
  return success, error

def _readAimsOutput(inputFile, system, relaxed=True):
  """
  Reads in FHI-aims output as a system.
  """
  
  error = ""
  success = True
  have_nice_day = False
  
  readAtoms = False
  readCompleted = False
  
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
    
  fin.close()
          
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
    
    system.totalEnergy = energy
    system.energyDefinition = "FHI-aims_" + version
    
    system.noOfcores = noOfcores
    system.runTime = runTime
    
  else:
    success = False
    error = __name__ + ": data has not been read from: " + inputFile
  
  return success, error
