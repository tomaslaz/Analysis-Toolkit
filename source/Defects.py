"""
Defects module.

@author Tomas Lazauskas
"""

import copy
import numpy as np

from .c_libs import defects as defects_c

def find_defects(input_system, final_system, vac_radius):
  """
  Compares two systems and tries to identify defects. This is based on the algorithm developed by Chris Scott
  at Loughborough University for LAKMC software package.
  
  """
  
  success = True
  error = ""
  
  includeVacs = 1
  includeInts = 1
  includeAnts = 1
  
  # make temporary list to store defects
  defectCluster = np.zeros(final_system.NAtoms, np.int32)
  
  NDefectsByType = np.zeros(4, np.int32)
  vacancies = np.empty(input_system.NAtoms, np.int32)
  antisites = np.empty(input_system.NAtoms, np.int32)
  onAntisites = np.empty(input_system.NAtoms, np.int32)
  interstitials = np.empty(final_system.NAtoms, np.int32)
  
  forcedSpeciesList = input_system.specieList
  
  # ignore these atoms as defects, by the fault everything is taken into account
  inputStateExclSpecs = np.empty(0, np.int32)
  refExclSpecs = np.empty(0, np.int32)
  
  cellDims3 = np.empty(3, np.float64)
  
  # here is the place where we can adjust the cell dimensions
  cellDims3[0] = final_system.cellDims[0]
  cellDims3[1] = final_system.cellDims[1]
  cellDims3[2] = final_system.cellDims[2]
  
  PBC = np.ones(3, np.int32)
  
  # checking whether the defects are in the same volume
  defectNeighbourRadius = 10.0
  
  input_system.minMaxPos(PBC)
  
  defects_cnt = defects_c.findDefects(includeVacs, includeInts, includeAnts,
                                      defectCluster, NDefectsByType, 
                                      vacancies, interstitials, antisites, onAntisites,
                                      forcedSpeciesList, inputStateExclSpecs, refExclSpecs,
                                      final_system.NAtoms, final_system.specieList, final_system.specie, final_system.pos,
                                      input_system.NAtoms, input_system.specieList, input_system.specie, input_system.pos,
                                      cellDims3, PBC, vac_radius, defectNeighbourRadius, input_system.minPos, input_system.maxPos, 
                                      3, 1)
  
   # resize arrays
  NDef = NDefectsByType[0]
  NVac = NDefectsByType[1]
  NInt = NDefectsByType[2]
  NAnt = NDefectsByType[3]
  vacancies.resize(NVac)
  interstitials.resize(NInt)
  antisites.resize(NAnt)
  onAntisites.resize(NAnt)
  defectCluster.resize(NDef)
  
  
    
  return success, error

