"""
GULP module.

@author Tomas Lazauskas
"""

import os
import numpy as np

import IO

_constIniParams = "Cell parameters (Angstroms/Degrees):"
_constIniParamsCellVol = "Initial cell volume"
_constOutOptiAchieved = "**** Optimisation achieved ****"
_constOutFinalCartCoords = "Final cartesian coordinates of atoms :"
_constOutFinalFracCoords = "Final fractional coordinates of atoms :"
_constOutFinalEnergy = "Final energy ="
_constOutFinalParams = "Final cell parameters and derivatives"
_constOutFinalDerivs = "Final internal derivatives :"

def readGulpOutput(system, fileName):
  """
  Reads the gulp output and updates atoms' positions and systems energy
  
  """
  
  success = False
  error = ""
    
  if not os.path.isfile(fileName):
    error = "File [%s] doesn't exist." % (fileName)
    return success, error, system
  
  try:
    f = open(fileName)
  except:
    error = "Cannot read file [%s]" % (fileName)
    return success, error, system
  
  # optimisation achieved?
  success = IO.stringInFile(_constOutOptiAchieved, f)
  
  if not success:
    error = "Optimisation has not been achieved"
  
  optiAchievedSection = False
  finalEnergy_eV = None
  
  iniParamsSectionSt = False
  iniParamsCount = 0
  
  finalCoordsSectionSt = False
  finalParamsSectionSt = False
  finalParamsSectionBreakCount = 0
  finalParamsCount = 0
  
  i = 0
  
  # read in final atom positions
  if success:
    for line in f:
      line = line.strip()
      
      # Looking for the final energy
      if optiAchievedSection:
        if _constOutFinalEnergy in line:
          energyArray = line.split()
          finalEnergy_eV = float(energyArray[3])
          
          # updating the final energy value
          system.totalEnergy = finalEnergy_eV
          system.energyDefinition = "Gulp"
          
          optiAchievedSection = False
      
       # Found the end of the section
      if _constOutFinalDerivs in line:
        finalParamsSectionSt = False
      
      # Found the end of the section
      if _constIniParamsCellVol in line:
        iniParamsSectionSt = False
      
      if finalCoordsSectionSt:
        array = line.split()
        
        if (len(array) == 7):
          atomSpecie = array[1]
          atomType =  array[2]

          # Do atom speces and types match?
          if ((not system.specieList[system.specie[i]].lower() == atomSpecie.lower()) or
              (not system.gulpAtomType[i] == atomType)):
            
            error = "atom types do not match %d : %s != %s" % (i, system.specieList[system.specie[i]], atomSpecie)
            return False, error
          
          # Reading the atom positions
          posx = float(array[3])
          posy = float(array[4])
          posz = float(array[5])
          
          # Updating the atom positions
          system.pos[3*i]   = posx
          system.pos[3*i+1] = posy
          system.pos[3*i+2] = posz
                     
          i += 1
          
          if i == system.NAtoms:
            finalCoordsSectionSt = False
                
      if finalParamsSectionSt:
        
        if "-----" in line:
          finalParamsSectionBreakCount += 1
        
        if finalParamsSectionBreakCount < 2:
          array = line.split()
          
          if len(array) == 6:
            # new cell dimensions
            if finalParamsCount == 0:
              system.cellDims_final[0] = np.float64(array[1])
              
            elif finalParamsCount == 1:
              system.cellDims_final[1] = np.float64(array[1])
            
            elif finalParamsCount == 2:
              system.cellDims_final[2] = np.float64(array[1])
            
            # new cell angles:
            elif finalParamsCount == 3:
              system.cellAngles_final[0] = np.float64(array[1])
            
            elif finalParamsCount == 4:
              system.cellAngles_final[1] = np.float64(array[1])
            
            elif finalParamsCount == 5:
              system.cellAngles_final[2] = np.float64(array[1])

            finalParamsCount += 1
        
        else:
          finalParamsSectionSt = False
      
      if iniParamsSectionSt:
        array = line.split()

        if len(array) == 6:
          # new cell dimensions
          if iniParamsCount == 0:
            system.cellDims_ini[0] = np.float64(array[2])
            system.cellAngles_ini[0] = np.float64(array[5])
            
          elif iniParamsCount == 1:
            system.cellDims_ini[1] = np.float64(array[2])
            system.cellAngles_ini[1] = np.float64(array[5])
            
          elif iniParamsCount == 2:
            system.cellDims_ini[2] = np.float64(array[2])
            system.cellAngles_ini[2] = np.float64(array[5])
          
          iniParamsCount += 1

      # Found the beginning of the section
      if ((_constOutFinalCartCoords in line) or (_constOutFinalFracCoords in line)):
        finalCoordsSectionSt = True
      
      # Found the beginning of the section
      if _constOutFinalParams in line:
        finalParamsSectionSt = True
 
       # Found the beginning of the section
      if _constIniParams in line:
        iniParamsSectionSt = True

      if _constOutOptiAchieved in line:
        optiAchievedSection = True
  
  success = True
  return success, error