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


def _readAimsStructure(geometryFile, outputFile):
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
    
    success, error = _readAimsOutput(outputFile, system)
      
  return success, error, system

def _readAimsOutput(inputFile, system):
  """
  Reads in FHI-aims output as a system.
  """
  
  error = ""
  success = True
  
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
    
    if ((len(fields) > 2) and (' '.join(fields[0:3]) == "Final atomic structure:")):
      readAtoms = True
    
    if ((len(fields) == 4) and (fields[0] == "Using")):
      noOfcores = int(fields[1])
    
    if ((len(fields) > 5) and (' '.join(fields[1:4]) == "Total time :")):
      runTime = float(fields[6])
    
  fin.close()
  
  if not success:
    return success, error
  
  if not readCompleted:
    success = False
    error = __name__ + ": data has not been read from: " + inputFile + ". It seems FHIaims encountered an error."
    
  elif system is not None:
    system.totalEnergy = energy
    system.energyDefinition = "FHI-aims_" + version
    
    system.noOfcores = noOfcores
    system.runTime = runTime
    
  else:
    success = False
    error = __name__ + ": data has not been read from: " + inputFile
  
  return success, error

  