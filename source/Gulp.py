"""
GULP module.

@author Tomas Lazauskas
"""

import os

import IO

_constOutOptiAchieved = "**** Optimisation achieved ****"
_constOutFinalCartCoords = "Final cartesian coordinates of atoms :"
_constOutFinalFracCoords = "Final fractional coordinates of atoms :"
_constOutFinalEnergy = "Final energy ="

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
  
  finalCoordsSectionSt = False
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
                
      # Found the beginning of the section
      if ((_constOutFinalCartCoords in line) or (_constOutFinalFracCoords in line)):
        finalCoordsSectionSt = True
      
      if _constOutOptiAchieved in line:
        optiAchievedSection = True
  
  success = True
  return success, error