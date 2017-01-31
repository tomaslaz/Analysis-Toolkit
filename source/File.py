"""
File module.

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

def writeATS(system, outputFile, radius):
  """
  Writes system as an ATS file.
  
  """
  
  # TODO: Need a way set the surface radius for a each atom. Probably add it as 
  #       as a new column in the atoms.in file.
  
  atoms_cnt = 0
  error = ""
  success = True
  
  if system is None:
    success = False
    error = __name__ + ": no data to write"
    
    return success, error
  
  if (len(system.specieList) < 1):
    success = False
    error = __name__ + ": no data to write"
    
    return success, error
  
  try:
    fout = open(outputFile, "w")
    
  except:
    success = False
    error = __name__ + ": Cannot open: " + outputFile
    
    return success, error
  
  group_name = "LE"
  
  for i in range(system.NAtoms):
    atoms_cnt += 1
    
    fout.write("%16f %16f %16f %5f %8s %d \n" % 
               (system.pos[3*i], system.pos[3*i+1], system.pos[3*i+2], 
                radius, group_name, atoms_cnt))
    
  fout.close()
    
  return success, error
  
def writeCAR(system, outputFile):
  """
  Writes system as a CAR file.
  
  """
  
  error = ""
  success = True
  
  if system is None:
    success = False
    error = __name__ + ": no data to write"
    
    return success, error
  
  try:
    fout = open(outputFile, "w")
  except:
    success = False
    error = __name__ + ": Cannot open: " + filePath
    
    return success, error
  
  fout.write("%s\n" % "!BIOSYM archive 3")
 
  if (system.PBC[0] <> 0 or system.PBC[0] <> 0 or system.PBC[0] <> 0):
    success = False
    error = __name__ + ": PBC are not implemented for CAR"
  
    return success, error
  
  fout.write("%s\n" % "PBC=OFF")
    
  fout.write("\n")
  fout.write("%s\n" % "!DATE")
  
  for i in range(system.NAtoms):
    
    tempStr =  ("%s%d" % (system.specieList[system.specie[i]], i+1))
    tempStr = "{:<7}".format(tempStr)
    
    fout.write("%7s %13.10f %13.10f %13.10f XXXX 1      xx      %2s %.4f\n" % (tempStr, 
      system.pos[3*i], system.pos[3*i+1], system.pos[3*i+2], 
      "{:<2}".format(system.specieList[system.specie[i]]), system.charge[i]))
  
  fout.write("end\n")
  fout.write("end\n")
  
  fout.close()
  
  return success, error

def writeGIN(system, outputFile, controlFile=None, outputXYZ=False):
  """
  Writes system as a GIN file.
  
  """
  
  error = ""
  success = True
  
  if (controlFile is None):
    masterGinFile = "Master.gin"
  else:
    masterGinFile = controlFile
  
  if system is None:
    success = False
    error = __name__ + ": no data to write"
    
    return success, error
  
  if (not os.path.isfile(masterGinFile)):
    success = False
    error = __name__ + ": could not locate Master.gin file"

    return success, error
    
  masterSystem, error, header, footer = readSystemFromFileGIN(masterGinFile, outputMode=True)
  
  if (masterSystem is None):
    success = False
    return success, error
  
  try:
    fout = open(outputFile, "w")
  except:
    success = False
    error = __name__ + ": Cannot open: " + filePath
     
    return success, error
  
  fout.write(header)

  for i in range(system.NAtoms):
    tempStr = system.specieList[system.specie[i]]
    
    gulpSpecies = masterSystem.gulpSpecies[tempStr]
    
    gulpSpeciesArr = gulpSpecies.split(";")
    
    for j in range(len(gulpSpeciesArr)):
      if len(gulpSpeciesArr[j].strip()) > 0:
        
        gulpSpeciesArr2 = gulpSpeciesArr[j].split(",")
                
        type = gulpSpeciesArr2[0].strip()
        
        try:
          charge = float(gulpSpeciesArr2[1].strip())
        except:
          charge = 0.0

        fout.write("%s %s %13.10f %13.10f %13.10f %.4f\n" % 
                  (tempStr, type, system.pos[3*i], system.pos[3*i+1], system.pos[3*i+2], charge)) 
     
  fout.write(footer)
  
  # adding output to xyz
  if outputXYZ:
    fout.write("\noutput xyz %s" % (outputFile[:-4]))

  fout.close()
    
  return success, error

def writeAimsGeometry(system, outputFile):
  """
  Writes system as a geometry file.
  
  """
  
  error = ""
  success = True
  
  if system is None:
    success = False
    error = __name__ + ": no data to write"
    
    return success, error
  
  try:
    fout = open(outputFile, "w")
  except:
    success = False
    error = __name__ + ": Cannot open: " + filePath
    
    return success, error
  
  for i in range(system.NAtoms):
    fout.write("atom %.10f %.10f %.10f %s\n" % (system.pos[3*i], system.pos[3*i+1], system.pos[3*i+2], 
      system.specieList[system.specie[i]]))
  
  fout.close()
  
  return success, error

def writeXYZ(system, outputFile):
  """
  Writes system as an XYZ file.
  
  """
  
  error = ""
  success = True
  
  if system is None:
    success = False
    error = __name__ + ": no data to write"
    
    return success, error
  
  if (len(system.specieList) < 1):
    success = False
    error = __name__ + ": no data to write"
    
    return success, error
  
  try:
    fout = open(outputFile, "w")
    
  except:
    success = False
    error = __name__ + ": Cannot open: " + outputFile
    
    return success, error
  
  fout.write("%d\n" % system.NAtoms)
  
  metaData = "SCF Done             %.10e;" % (system.totalEnergy)
    
  fout.write("%s\n" % (metaData))
  for i in range(system.NAtoms):
        
    fout.write("%s %.10f %.10f %.10f %.2f\n" % (system.specieList[system.specie[i]], 
      system.pos[3*i], system.pos[3*i+1], system.pos[3*i+2], system.charge[i]))
  
  fout.close()
    
  return success, error