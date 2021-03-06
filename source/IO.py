"""
Input/Output module.

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import copy
import os
import sys
import glob
import time
import numpy as np

# import System
# import Gulp
# import Atoms
import System
import Gulp
import Atoms
#import Fhiaims

const_file_ext_xyz = "xyz"
const_file_ext_out = "out"

def checkDirectory(dirPath, createMd=0):
  """
  Checks if directory exists
  """
  
  exists = os.path.exists(dirPath)
  
  if ((not exists) and (createMd)):
    try:
      os.makedirs(dirPath)
      exists = True
    except:
      exists = False
      
  return exists

def checkFile(filePath):
  """
  Checks if file exists.
  
  """
  
  return os.path.isfile(filePath)

def countLines(fileName):
  """
  Counts the number of lines in a file
  
  """

  success = False
  error = ""
  linesCount = 0
  
  if not os.path.isfile(fileName):
    error = "File [%s] doesn't exist." % (fileName)
    return success, error, linesCount
  
  try:
    f = open(fileName)
    
  except:
    error = "Cannot read file [%s]" % (fileName)
    return success, error, linesCount
  
  for line in f:
    linesCount += 1
  
  f.close()
  
  success = True
  return success, error, linesCount
  
def countMixAtoms(fileName):
  """
  Counts the number of atoms to be used in mixing.
  
  """
  
  success = False
  error = ""
  atomsCnt = 0
  
  if not os.path.isfile(fileName):
    error = "File [%s] doesn't exist." % (fileName)
    return success, error, atomsCnt
  
  try:
    f = open(fileName)
  except:
    error = "Cannot read file [%s]" % (fileName)
    return success, error, atomsCnt
    
  regionSt = False
  regionRead = True
  
  for line in f:
    
    if (("extra" in line) or ("species" in line)):
      regionSt = False
      regionRead = False
    
    elif regionRead and regionSt:
      if "#" not in line:
        atomsCnt += 1
    
    if (("fractional" in line) or ("cartesian" in line)):
      regionSt = True
    
  f.close()
  
  #print "Counting the number of atoms (including vacancies): %d" % (atomsCnt)
  
  success = True
  return success, error, atomsCnt

def get_unique_systems_hashkeys(systems_list, hashkeyRadius=None):
  """
  Evaluates the hashkeys
  
  """
  
  system_list_len = len(systems_list)
  
  unique_hashkeys = []
  unique_systems = []
  unique_cnt = 0
  
  hashkeys_count = np.zeros(system_list_len, np.int16)
  
  system_cnt = 0
  for system in systems_list:
    
    systems_list[system_cnt].calculateHashkey(hashkeyRadius=hashkeyRadius)
    hashkey = systems_list[system_cnt].hashkey
      
    try:
      hashkey_idx = unique_hashkeys.index(hashkey)
    except:
      hashkey_idx = -1
    
    # storing duplicate numbers
    if hashkey_idx != -1:
      hashkeys_count[hashkey_idx] += 1
      
      prevEnergy = unique_systems[hashkey_idx].totalEnergy
      currentEnergy = systems_list[system_cnt].totalEnergy
      
      if currentEnergy < prevEnergy:
        unique_systems[hashkey_idx] = copy.deepcopy(systems_list[system_cnt])
        
    # saving only unique: 
    else:
      unique_cnt += 1
      
      unique_hashkeys.append(hashkey)
      unique_systems.append(systems_list[system_cnt])
    
    system_cnt += 1
    if (system_cnt % 100 == 0): print ("Getting the hashkeys %d/%d" % (system_cnt, system_list_len))
  
  for i in range(unique_cnt):
    unique_systems[i].hashkey_duplicate_cnt = hashkeys_count[i]
  
  return unique_systems
  
def get_file_list(extension="*"):
  """
  Returns a list of files with a specific extension
  
  """
  
  file_list = glob.glob("*.%s" % (extension))
  
  return file_list

def get_file_list_recursive_simple(extension="*", file_list=[]):
  """
  Returns a list of files with a specific extension by analysing directories recursively
  
  """

  cwd = os.getcwd()
  
  for root, dirs, files in os.walk("./"):
    
    dir_path_add = os.path.join(cwd, root[2:])
    os.chdir(dir_path_add)
    
    # looking for files in a directory
    file_list_dir = glob.glob("*.%s" % (extension))
    
    # adding directory to the path 
    for i in range(len(file_list_dir)):
      file_list_dir[i] = os.path.join(dir_path_add, file_list_dir[i])
       
    # joining two lists, python makes it so easy :)
    file_list += file_list_dir
    
    os.chdir(cwd)
        
  return file_list

def get_file_list_recursive(extension="*", file_list=[], dir_path=None, recurs_iter=0, recurs_max=8):
  """
  Returns a list of files with a specific extension by analysing directories recursively
  
  """

  cwd = os.getcwd()
  if dir_path is not None:
    os.chdir(dir_path)
  
  dir_path_add = os.getcwd()
  
  for _, dirs, files in os.walk("./"):
    
    if len(dirs) > 0:
      for dir_recur in dirs:
        file_list_dir = get_file_list_recursive(extension=extension, file_list=file_list, 
                                      dir_path=dir_recur, recurs_iter=recurs_iter+1, recurs_max=recurs_max)
      
    else:
      file_list_dir = glob.glob("*.%s" % (extension))
      
      # adding directory to the path 
      for i in range(len(file_list_dir)):
        file_list_dir[i] = os.path.join(dir_path_add, file_list_dir[i])
        
      # joining two lists, python makes it so easy :)
      file_list += file_list_dir
        
  os.chdir(cwd)
    
  return file_list

def get_file_list_pre_sub(prefix="*", subfix="*"):
  """
  Returns a list of files with specified prefix and/or subfix
  
  """
  
  file_list = glob.glob("%s*.%s" % (prefix, subfix))
  
  return file_list

def get_dir_list(file_name=None):
  """
  Get a list of directories which contain file_name file.
  
  """
  
  dir_list = []
  
  for root, dirs, files in os.walk("./"):
    
    if file_name is not None:
      if file_name in files:
        
        dir_list.append(root)
    else:
      dir_list.append(root)
        
  return dir_list
  
def lookForFiles(extension):
  """
  Looks for files with a specific extension and returns the last one (according to the file system)
  
  """
  
  csvFile = None
  
  for file in glob.glob("*.%s" % (extension)):
    csvFile = file
  
  return csvFile

def removeFile(filePath):
  """
  Deletes a file.
  
  """
  
  success = True
  error = ""
  
  try:
    os.unlink(filePath)
    
  except:
    error = "Cannot delete file: %s" % (filePath)
    success = False
  
  return success, error

def read_in_systems(systems_paths_list):
  """
  Reads in systems form a list of paths and returns a system list
  
  """
  
  systems = []
  systems_paths_cnt = len(systems_paths_list)
  
  systems_paths_iter = 1
  for system_file in systems_paths_list:
    
    try:
      system, error = readSystemFromFile(system_file)
    except:
      system = None
      error = "cannot read in the file"
    
    if system is not None:
      systems.append(system)
    else:
      print ("error: %s" % (error))
     
    if (systems_paths_iter % 1000 == 0):
      print ("Reading %d/%d" % (systems_paths_iter, systems_paths_cnt))
    
    systems_paths_iter += 1
    
  return systems

def readGulpOutputPolymerGulpOutput(fileName):
  """
  Reads input information from a gulp ouput of a polymer simulation
  
  """

  polymer = None 
  
  # counting the number of cores
  success, error, coreCount = Gulp.readGulpOutputPolymerInputCountCoresShells(fileName)
  
  if not success:
    return success, error, polymer
  
  polymer = System.System(coreCount)
  
  # reading in the initial polymer coordinates and system parameters
  success, error = Gulp.readGulpOutputPolymerInput(polymer, fileName)
  
  if not success:
    return success, error, polymer
  
  polymer.iniPos = copy.deepcopy(polymer.pos)
  
  # reading the final coordinates and final cell parameter
  success, error = Gulp.readGulpOutputPolymerOutput(polymer, fileName)
  
  return success, error, polymer

def readSystemFromFile(file_name):
  """
  Read a system from a file
  
  """
  
  system = None
  error = ""
  
  if file_name.endswith(".gin"):
    system, error = readSystemFromFileGIN(file_name)

  elif file_name.endswith(".car"):
    system = readSystemFromFileCAR(file_name)
    
    if system is None:
      error = "System cannot be read from file: %s" % (file_name)
    
  elif file_name.endswith(".xyz"):
    system = readSystemFromFileXYZ(file_name)
    
    if system is None:
      error = "System cannot be read from file: %s" % (file_name)
  
  elif file_name.endswith(".out"):
    
    print ("Trying to read [%s]." % (file_name))
    
    cwd = os.getcwd()
        
    # changing into fhiaims simulation dir
    os.chdir(os.path.dirname(file_name))
    
    # getting fhiaims output file name
    aims_output_file = os.path.basename(file_name)
        
    success, error, system = Fhiaims._readAimsStructure(Fhiaims._const_geometry_in, 
                                                        aims_output_file, 
                                                        relaxed=False, eigenvalues=False)
    
    os.chdir(cwd)
      
  else:
    error = "Unidentified file format"
  
  return system, error

def readSystemFromFileARC(fileName):
    """
    Reads in the structure of a system from an ARC file.
    
    """
    
    system = None
    
    if not os.path.isfile(fileName):
        print ("File [%s] doesn't exist." % (fileName))
        return system
    
    try:
        f = open(fileName)
    except:
        print ("Cannot read file [%s]" % (fileName))
        return system
    
    i = 0
    totAtomCnt = 0
    
    for line in f:
      
      i += 1
      line = line.strip()
      
      if ("end" in line):
        break
      
      if (i > 5):
        totAtomCnt += 1
      
    f.close()
    
    system = System.System(totAtomCnt)
    
    # Adding system name
    system.name = os.path.splitext(os.path.basename(fileName))[0]
    
    f = open(fileName)
    
    i = 0
    atomCnt = 0
    
    for line in f:
      
      i += 1
      line = line.strip()
      
      if ("end" in line):
        break
      
      if (i > 5):
        array = line.split()
        
        sym = array[7].strip()
        
        if sym not in system.specieList:
          system.addSpecie(sym)
       
        specInd = system.specieIndex(sym)
         
        system.specieCount[specInd] += 1
         
        system.specie[atomCnt] = specInd
         
        for j in range(3):
          system.pos[atomCnt*3 + j] = float(array[j+1])
        
        system.charge[atomCnt] = float(array[8])
        atomCnt += 1
      
      # reading in lattice parameters
      elif (i == 5):
        array = line.split()
        
        # cell dimensions
        system.cellDims[0] = float(array[1])
        system.cellDims[1] = float(array[2])
        system.cellDims[2] = float(array[3])
        
        # cell angles
        system.cellAngles[0] = array[4]
        system.cellAngles[1] = array[5]
        system.cellAngles[2] = array[6]
        
    f.close()

    return system

def readSystemFromFileCAR(fileName):
    """
    Reads in the structure of a system from a CAR file.
    
    """
    
    system = None
    
    if not os.path.isfile(fileName):
        print ("File [%s] doesn't exist." % (fileName))
        return system
    
    try:
        f = open(fileName)
    except:
        print ("Cannot read file [%s]" % (fileName))
        return system
    
    i = 0
    totAtomCnt = 0
    
    for line in f:
      
      i += 1
      line = line.strip()
      
      if ("end" in line):
        break
      
      if (i >= 5):
        totAtomCnt += 1
      
    f.close()
    
    system = System.System(totAtomCnt)
    
    f = open(fileName)
    
    i = 0
    atomCnt = 0
    
    for line in f:
      
      i += 1
      line = line.strip()
      
      if ("end" in line):
        break
      
      if (i >= 5):
        array = line.split()
        
        sym = array[7].strip()
                 
        if sym not in system.specieList:
            system.addSpecie(sym)
       
        specInd = system.specieIndex(sym)
         
        system.specieCount[specInd] += 1
         
        system.specie[atomCnt] = specInd
         
        for j in range(3):
            system.pos[atomCnt*3 + j] = float(array[j+1])
                
        system.charge[atomCnt] = float(array[8])
        atomCnt += 1
    
    f.close()

    return system

def readSystemFromFileGIN(fileName, outputMode=False):
  """
  Reads in the structure of a system from a GIN file.
  
  """
  
  header = ""
  footer = ""
  
  system = None
  
  success, error, NAtoms = countMixAtoms(fileName)
    
  if not success:
    if outputMode:
      return system, error, header, footer
    
    else:
      return system, error
  
  if not os.path.isfile(fileName):
    error = "File [%s] doesn't exist." % (fileName)
    
    if outputMode:
      return system, error, header, footer
    
    else:
      return system, error
  
  try:
    f = open(fileName)
    
  except:
    error = "Cannot read file [%s]" % (fileName)
    
    if outputMode:
      return system, error, header, footer
    
    else:
      return system, error
    
  system = System.System(NAtoms)
    
  atomsCnt = 0
  regionSt = False
  regionRead = True
  regionReached = False
  
  speciesSt = False
  speciesEnd = False
  speciesCnt = 0
  
  getCellParams = False
    
  coordTypeRead = False
  
  for line in f:
    
    if outputMode:
      if (not regionRead and not regionSt):
        footer += line
        
      elif (regionRead and not regionSt):
        header += line
        
    if getCellParams:
      line = line.strip()
      array = line.split()
      
      system.cellDims[0] = array[0]
      system.cellDims[1] = array[1]
      system.cellDims[2] = array[2]
      
      system.cellAngles[0] = array[3]
      system.cellAngles[1] = array[4]
      system.cellAngles[2] = array[5]
      
      system.PBC[0] = 1
      system.PBC[1] = 1
      system.PBC[2] = 1
      
      #TODO: extend to not only cubic systems
      if ((system.cellAngles[0] != 90.0) or (system.cellAngles[1] != 90.0) or (system.cellAngles[2] != 90.0)):
        print (" : At the moment we can deal only with cubic cells. Please contact if you want to use different systems")
        
      else:
        sys.exit(__name__ +" : Cannot determine cell type")
      
      getCellParams = False
      
    elif (("extra" in line) or ("species" in line)):
      regionSt = False
      regionRead = False
      
      if ("species" in line):
        speciesSt = True
        
        array = line.split()

        speciesCnt = int(array[1])
        
    elif (speciesSt and (speciesCnt > 0)):
      # Reading in the species
      
      line = line.strip()
      array = line.split()
      
      if len(array) > 2:
        try:
          system.gulpSpecies[array[0].strip()]
          system.gulpSpecies[array[0].strip()] += array[1].strip() + "," + array[2].strip() + ";"
        except:
          system.gulpSpecies[array[0].strip()] = array[1].strip() + "," +  array[2].strip() + ";"
      else:
        try:
          system.gulpSpecies[array[0].strip()]
          system.gulpSpecies[array[0].strip()] += array[1].strip() + ";"
        except:
          system.gulpSpecies[array[0].strip()] = array[1].strip() + ";"
        
      speciesCnt -= 1
    
    elif regionRead and regionSt:

      line = line.strip()
      array = line.split()
      array_len = len(array)
      
      if "#" not in line:
                
        #print line
        # reading the specie
        sym = array[0].strip()
        
        if sym not in system.specieList:
          system.addSpecie(sym)
        
        specInd = system.specieIndex(sym)
        system.specieCount[specInd] += 1
        system.specie[atomsCnt] = specInd
        
        system.gulpAtomType[atomsCnt] = array[1][0]
        
        # positions
        for j in range(3):
          system.pos[atomsCnt*3 + j] = float(array[j+2])
        
        # charge 
        try:
          system.charge[atomsCnt] = float(array[5])
          
        except:
          system.charge[atomsCnt] = 0.0

        extraInfo = ""
        for j in range(5, array_len):
          extraInfo = "%s %s" % (extraInfo, array[j])
        
        system.gulpAtomExtraInfo[atomsCnt] = extraInfo
        
        atomsCnt += 1
      
    elif ("fractional" in line):
      
      regionSt = True
      regionReached = True

    elif ("cartesian" in line):
      
      regionSt = True
      regionReached = True
      
    if atomsCnt == NAtoms and atomsCnt > 0:
      regionSt = False
      regionRead = False

  f.close()
      
  if outputMode:
    return system, error, header, footer
  
  else:
    return system, error
  
def readSystemFromFileXYZ(fileName):
    """
    Reads in the structure of a system from an XYZ file.
    
    """
    
    system = None
    
    if not os.path.isfile(fileName):
        print ("File [%s] doesn't exist." % (fileName))
        return system
    
    try:
        f = open(fileName)
    except:
        print ("Cannot read file [%s]" % (fileName))
        return system
    
    line = f.readline().strip()
    
    NAtoms = int(line)
        
    system = System.System(NAtoms)
    
    # additional info
    line = f.readline().strip()
    line_array = line.split()
    line_array_len = len(line_array)

    # reading energy if it is a result from ...
    if "SCF Done" in line:
      energy = float(line_array[2].strip()[:-1])
      system.totalEnergy = energy
    
    elif "Generated by KLMC DM" in line:
      energy = float(line_array[4].strip())
      system.totalEnergy = energy
    
    elif "Generated by KLMC" in line:
      energy = float(line_array[3].strip())
      system.totalEnergy = energy
    
    elif (("energy=" in line) and (line_array_len == 3)):      
      energy = float(line_array[2].strip())
      system.totalEnergy = energy
    
    elif line_array_len == 3:      
      system.PBC = np.ones(3, np.int32)
      system.cellDims = np.array([float(line_array[0]), float(line_array[1]), float(line_array[2])], np.float64)
      system.cellAngles = np.array([90.0, 90.0, 90.0], np.float64)
      
    # atoms and their positions
    i = 0
    for line in f:
        array = line.strip().split()

        sym = array[0].strip()
        
        if sym not in system.specieList:
            system.addSpecie(sym)
        
        specInd = system.specieIndex(sym)
        
        system.specieCount[specInd] += 1
        
        system.specie[i] = specInd
        
        for j in range(3):
            system.pos[i*3 + j] = float(array[j+1])
        
        try:
            system.charge[i] = array[4]
        except:
            system.charge[i] = 0.0
        
        i += 1
        
        if i == NAtoms:
          break
    
    f.close()
    
    system.name = os.path.splitext(os.path.basename(fileName))[0]
    
    return system

def save_systems_to_xyz(systems_list, dir_path):
  """
  Saves systems into xyz files
  
  """
  
  rank = 1
  for system in systems_list:
    
    file_name = "n%02d_%03d_%s.xyz" % (system.NAtoms, rank, system.name)

    file_path = os.path.join(dir_path, file_name)
    success_, error_ = writeXYZ(system, file_path)
  
    rank += 1
    
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
    
#     gulpSpecies = masterSystem.gulpSpecies[tempStr]
#     
#     gulpSpeciesArr = gulpSpecies.split(";")
#     
#     for j in range(len(gulpSpeciesArr)):
#       if len(gulpSpeciesArr[j].strip()) > 0:
#         
#         gulpSpeciesArr2 = gulpSpeciesArr[j].split(",")
#                 
#         type = gulpSpeciesArr2[0].strip()
#         
#         try:
#           charge = float(gulpSpeciesArr2[1].strip())
#         except:
#           charge = 0.0

    fout.write("%s %s %13.10f %13.10f %13.10f %s\n" % 
              (tempStr, system.gulpAtomType[i], system.pos[3*i], system.pos[3*i+1], system.pos[3*i+2], system.gulpAtomExtraInfo[i])) 
     
  fout.write(footer)
  
  # adding output to xyz
  if outputXYZ:
    fout.write("\noutput xyz %s" % (outputFile[:-4]))

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
 
  if (system.PBC[0] != 0 or system.PBC[0] != 0 or system.PBC[0] != 0):
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

def writeXYZ(system, outputFile, scfDone=True):
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
  
  metaData = ""
  
  if system.cellDims[0] != 0.0 and system.cellDims[1] != 0.0 and system.cellDims[2] != 0.0:
    meta_data = "%.10f %.10f %.10f" % (system.cellDims[0], system.cellDims[1], system.cellDims[2])

  elif scfDone and system.totalEnergy != None and system.totalEnergy != 0.0:
    metaData = "SCF Done             %.10e;" % (system.totalEnergy)
  
  fout.write("%s\n" % (metaData))

  for i in range(system.NAtoms):
        
    fout.write("%s %.10f %.10f %.10f %.2f\n" % (system.specieList[system.specie[i]], 
      system.pos[3*i], system.pos[3*i+1], system.pos[3*i+2], system.charge[i]))
  
  fout.close()
    
  return success, error