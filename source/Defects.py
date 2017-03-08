"""
Defects module.

@author Tomas Lazauskas
"""

import copy
import numpy as np

from c_libs import defects as defects_c

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
  
  # saving the positions
  # build list of vacancies
  vacSpecie = np.empty(NVac, np.int32)
  vacPos = np.empty(3 * NVac, np.float64)
    
  for i in xrange(len(vacancies)):
    index = vacancies[i]
    
    vacSpecie[i] = input_system.specie[index]
    vacPos[3*i] = input_system.pos[3*index]
    vacPos[3*i+1] = input_system.pos[3*index+1]
    vacPos[3*i+2] = input_system.pos[3*index+2]
    
  # build list of interstitials
  intSpecie = np.empty(NInt, np.int32)
  intPos = np.empty(3 * NInt, np.float64)
  
  for i in xrange(len(interstitials)):
    index = interstitials[i]
    
    intSpecie[i] = final_system.specie[index]
    intPos[3*i] = final_system.pos[3*index]
    intPos[3*i+1] = final_system.pos[3*index+1]
    intPos[3*i+2] = final_system.pos[3*index+2]
      
  # build list of antisites
  antSpecie = np.empty(NAnt, np.int32)
  antPos = np.empty(3 * NAnt, np.float64)
  onAntPos = np.empty(3*NAnt, np.float64)
  onAntSpecie = np.empty(3 * NAnt, np.int32)
  
  for i in xrange(len(antisites)):
    index = antisites[i]
    index2 = onAntisites[i]
     
    antSpecie[i] = input_system.specie[index]
    onAntSpecie[i] = final_system.specie[index2]
    
    antPos[3*i] = input_system.pos[3*index]
    antPos[3*i+1] = input_system.pos[3*index+1]
    antPos[3*i+2] = input_system.pos[3*index+2]
    
    onAntPos[3*i] = final_system.pos[3*index2]
    onAntPos[3*i+1] = final_system.pos[3*index2+1]
    onAntPos[3*i+2] = final_system.pos[3*index2+2]
  
  final_system.NDef = NDef
  final_system.NVac = NVac
  final_system.NInt = NInt
  final_system.NAnt = NAnt
  final_system.vacSpecie = copy.deepcopy(vacSpecie)
  final_system.vacPos = copy.deepcopy(vacPos)
  final_system.intSpecie = copy.deepcopy(intSpecie)
  final_system.intPos = copy.deepcopy(intPos)
  final_system.antSpecie = copy.deepcopy(antSpecie)
  final_system.antPos = copy.deepcopy(antPos)
  final_system.onAntPos = copy.deepcopy(onAntPos)
  final_system.onAntSpecie = copy.deepcopy(onAntSpecie)
  final_system.vacancies = copy.deepcopy(vacancies)
  final_system.interstitials = copy.deepcopy(interstitials)
  final_system.antisites = copy.deepcopy(antisites)
  final_system.onAntisites = copy.deepcopy(onAntisites)
  final_system.defectCluster = copy.deepcopy(defectCluster)
  
  final_system.printDefectsPositions()
  
  return success, error

