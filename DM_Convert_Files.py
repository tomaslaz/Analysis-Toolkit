#!/usr/bin/env python

"""
A script to convert files. Works with XYZ, CAR, and GIN formats.

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import glob
import os
import sys
import numpy as np
from optparse import OptionParser

import Atoms

class Cluster(object):
  """
  A class to save the structure of a cluster.
  
  NAtoms: number of atoms in cluster (N)
  specie[N]: array of symbols of atoms
  pos[3N]: array of positions of atoms
  charge[N]: array of charges of atoms
  
  """
    
  def __init__(self, NAtoms):
      
    self.totalEnergy = 0.0
    self.NAtoms = NAtoms
    self.cellDims = np.zeros(3, np.float64)
    self.cellAngles = np.empty(3, np.float64)
    
    self.specie = np.empty(self.NAtoms, np.int32)
    self.pos = np.empty(3*self.NAtoms, np.float64)
    self.minPos = np.empty(3, np.float64)
    self.maxPos = np.empty(3, np.float64)
    
    self.charge = np.empty(self.NAtoms, np.float64)
    
    self.com = np.empty(3, np.float64)
    self.momentOfInertia = np.zeros([3, 3], np.float64)
    
    dt = np.dtype((str, 2))
    self.specieList = np.empty(0, dt)
    self.specieCount = np.empty(0, np.int32)
    
    self.PBC = np.zeros(3, np.int32)
    
    self.gulpSpecies = {}
        
  def addAtom(self, sym, pos, charge):
    """
    Add an atom to the cluster
    
    """
    
    if sym not in self.specieList:
        self.addSpecie(sym)
    
    specInd = self.specieIndex(sym)
    
    self.specieCount[specInd] += 1
    
    pos = np.asarray(pos, dtype=np.float64)
    
    self.specie = np.append(self.specie, np.int32(specInd))
    self.pos = np.append(self.pos, pos)
    self.charge = np.append(self.charge, charge)

    self.NAtoms += 1
    
  def calcCOM(self):
      
    """
    Calculates the centre of mass of a cluster
    
    """
    
    totMass = 0.0
    self.com[:] = 0.0
    
    for i in range(self.NAtoms):
      atomMass = Atoms.atomicMassAMU(self.specieList[self.specie[i]])
      totMass += atomMass
      
      for j in range(3):
        self.com[j] += atomMass * self.pos[3*i + j]

    self.com = self.com / totMass
  
  def calcMOI(self):
    """
    Calculates moment of inertia
    
    """
    
    moi = np.zeros(6, np.float64)

    
    self.momentOfInertia[:] = 0.0
    
    for i in range(self.NAtoms):
      atomMass = Atoms.atomicMassAMU(self.specieList[self.specie[i]])
        
      moi[0] += (self.pos[3*i+1]**2 + self.pos[3*i+2]**2) * atomMass
      moi[1] += (self.pos[3*i+0]**2 + self.pos[3*i+2]**2) * atomMass
      moi[2] += (self.pos[3*i+0]**2 + self.pos[3*i+1]**2) * atomMass
      moi[3] += -(self.pos[3*i+0] * self.pos[3*i+1]) * atomMass
      moi[4] += -(self.pos[3*i+0] * self.pos[3*i+2]) * atomMass
      moi[5] += -(self.pos[3*i+1] * self.pos[3*i+2]) * atomMass
    
    self.momentOfInertia[0][0] = moi[0]
    self.momentOfInertia[1][1] = moi[1]
    self.momentOfInertia[2][2] = moi[2]
    
    self.momentOfInertia[0][1] = moi[3]
    self.momentOfInertia[0][2] = moi[4]
    
    self.momentOfInertia[1][0] = moi[3]
    self.momentOfInertia[1][2] = moi[5]
    
    self.momentOfInertia[2][0] = moi[4]
    self.momentOfInertia[2][1] = moi[5]
  
  def rotateToMOI(self, basis):
    
    for i in range(self.NAtoms):
      for j in range(3):
        elSum = 0.0
        
        for k in range(3):
          elSum += self.pos[3*i+k] * basis[k][j]

        self.pos[3*i+j] = elSum
      
  def moveToCOM(self):
    """
    Centers the cluster on the centre of of mass.
    
    """
    
    for i in range(self.NAtoms):
      for j in range(3):
        self.pos[3*i + j] -= self.com[j]
      
  def removeAtom( self, index ):
    """
    Remove an atom from the structure
    
    """
    
    specInd = self.specie[index]
    self.specie = np.delete(self.specie, index)
    self.pos = np.delete(self.pos, [3*index,3*index+1,3*index+2])
    self.charge = np.delete(self.charge, index)

    self.NAtoms -= 1
    
    self.specieCount[specInd] -= 1
    if self.specieCount[specInd] == 0 and not self.specieListForced:
        self.removeSpecie(specInd)
  
  def removeSpecie(self, index):
    """
    Remove a specie from the specie list.
    
    """
    self.specieCount = np.delete(self.specieCount, index)
    self.specieList = np.delete(self.specieList, index)
    
    for i in xrange(self.NAtoms):
        if self.specie[i] > index:
            self.specie[i] -= 1

  def specieIndex(self, check):
    """
    Index of sym in specie list
    
    """
    
    count = 0
    index = -1
    for sym in self.specieList:
        if sym == check:
            index = count
            break
        
        count += 1
    
    return index 
  
  def addSpecie(self, sym, count=None):
    """
    Add specie to specie list
    
    """
            
    if sym in self.specieList:
        if count is not None:
            specInd = self.specieIndex(sym)
            self.specieCount[specInd] = count
        
        return
    
    if count is None:
        count = 0
    
    self.specieList = np.append(self.specieList, sym)
    self.specieCount = np.append(self.specieCount, np.int32(count))
    
  def minMaxPos(self, PBC):
      
    for i in xrange(3):
        if not PBC[i]:
            self.minPos[i] = self.pos[i::3].min()
            self.maxPos[i] = self.pos[i::3].max()

def cmdLineArgs():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog inputFile outputFile"
  
  parser = OptionParser(usage=usage)

  parser.disable_interspersed_args()
      
  (options, args) = parser.parse_args()

  if (len(args) != 2):
      parser.error("incorrect number of arguments")

  return options, args

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
    
      atomsCnt += 1
    
    if (("fractional" in line) or ("cartesian" in line)):
      regionSt = True
    
  f.close()
  
  print "Counting the number of atoms (including vacancies): %d" % (atomsCnt)
  
  success = True
  return success, error, atomsCnt

def readClusterFromFileCAR(fileName):
    """
    Reads in the structure of a cluster from a CAR file.
    
    """
    
    cluster = None
    
    if not os.path.isfile(fileName):
        print "File [%s] doesn't exist." % (fileName)
        return cluster
    
    try:
        f = open(fileName)
    except:
        print "Cannot read file [%s]" % (fileName)
        return cluster
    
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
    
    cluster = Cluster(totAtomCnt)
    
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
         
        if sym not in cluster.specieList:
            cluster.addSpecie(sym)
       
        specInd = cluster.specieIndex(sym)
         
        cluster.specieCount[specInd] += 1
         
        cluster.specie[atomCnt] = specInd
         
        for j in range(3):
            cluster.pos[atomCnt*3 + j] = float(array[j+1])
        
        cluster.charge[atomCnt] = float(array[8])
        atomCnt += 1
    
    f.close()

    return cluster

def readClusterFromFileGIN(fileName, outputMode=False):
  """
  Reads in the structure of a cluster from a GIN file.
  
  """
  
  header = ""
  footer = ""
  
  cluster = None
  
  success, error, NAtoms = countMixAtoms(fileName)
  
  if not success:
    if outputMode:
      return cluster, error, header, footer
    
    else:
      return cluster, error
  
  if not os.path.isfile(fileName):
    error = "File [%s] doesn't exist." % (fileName)
    
    if outputMode:
      return cluster, error, header, footer
    
    else:
      return cluster, error
  
  try:
    f = open(fileName)
    
  except:
    error = "Cannot read file [%s]" % (fileName)
    
    if outputMode:
      return cluster, error, header, footer
    
    else:
      return cluster, error
    
  cluster = Cluster(NAtoms)
  
  atomsCnt = 0
  regionSt = False
  regionRead = True
  regionReached = False
  
  speciesSt = False
  speciesEnd = False
  speciesCnt = 0
  
  getCellParams = False
  
  for line in f:
    
    if outputMode:
      if (not regionRead and not regionSt):
        footer += line
        
      elif (regionRead and not regionSt):
        header += line
      
    if getCellParams:
      line = line.strip()
      array = line.split()
      
      cluster.cellDims[0] = array[0]
      cluster.cellDims[1] = array[1]
      cluster.cellDims[2] = array[2]
      
      cluster.cellAngles[0] = array[3]
      cluster.cellAngles[1] = array[4]
      cluster.cellAngles[2] = array[5]
      
      cluster.PBC[0] = 1
      cluster.PBC[1] = 1
      cluster.PBC[2] = 1
      
      #TODO: extend to not only cubic clusters
      if ((cluster.cellAngles[0] <> 90.0) or (cluster.cellAngles[1] <> 90.0) or (cluster.cellAngles[2] <> 90.0)):
        print " : At the moment we can deal only with cubic cells. Please contact if you want to use different clusters"
        
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

      try:
        cluster.gulpSpecies[array[0].strip()]
        cluster.gulpSpecies[array[0].strip()] += array[1].strip() + "," + array[2].strip() + ";"
      except:
        cluster.gulpSpecies[array[0].strip()] = array[1].strip() + "," +  array[2].strip() + ";"
        
      speciesCnt -= 1
    
    elif regionRead and regionSt:

      line = line.strip()
      array = line.split()
      
      # reading the specie
      sym = array[0].strip()
      
      if sym not in cluster.specieList:
        cluster.addSpecie(sym)
      
      specInd = cluster.specieIndex(sym)
      cluster.specieCount[specInd] += 1
      cluster.specie[atomsCnt] = specInd
      
      # positions
      for j in range(3):
        cluster.pos[atomsCnt*3 + j] = float(array[j+2])
      
      # charge 
      try:
        cluster.charge[atomsCnt] = float(array[5])
        
      except:
        cluster.charge[atomsCnt] = 0.0
      
      atomsCnt += 1
      
    elif ("fractional" in line):
      
      regionSt = True
      regionReached = True

    elif ("cartesian" in line):
      
      regionSt = True
      regionReached = True
      
    
    if atomsCnt == NAtoms:
      regionSt = False
      regionRead = False

  f.close()
      
  if outputMode:
    return cluster, error, header, footer
  
  else:
    return cluster, error
  
def readClusterFromFileXYZ(fileName):
    """
    Reads in the structure of a cluster from an XYZ file.
    
    """
    
    cluster = None
    
    if not os.path.isfile(fileName):
        print "File [%s] doesn't exist." % (fileName)
        return cluster
    
    try:
        f = open(fileName)
    except:
        print "Cannot read file [%s]" % (fileName)
        return cluster
    
    line = f.readline().strip()
    
    NAtoms = int(line)
        
    cluster = Cluster(NAtoms)
    
    # additional info
    line = f.readline().strip()
    
    # size ?
    array = line.split()
#     cluster.cellDims[0] = float(array[0])
#     cluster.cellDims[1] = float(array[1])
#     cluster.cellDims[2] = float(array[2])

    # atoms and their positions
    i = 0
    for line in f:
        array = line.strip().split()

        sym = array[0].strip()
        
        if sym not in cluster.specieList:
            cluster.addSpecie(sym)
        
        specInd = cluster.specieIndex(sym)
        
        cluster.specieCount[specInd] += 1
        
        cluster.specie[i] = specInd
        
        for j in range(3):
            cluster.pos[i*3 + j] = float(array[j+1])
        
        try:
            cluster.charge[i] = array[4]
        except:
            cluster.charge[i] = 0.0
        
        i += 1
        
        if i == NAtoms:
          break
    
    f.close()

    return cluster

def writeCAR(cluster, outputFile):
  """
  Writes cluster as a CAR file.
  
  """
  
  error = ""
  success = True
  
  if cluster is None:
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
 
  if (cluster.PBC[0] <> 0 or cluster.PBC[0] <> 0 or cluster.PBC[0] <> 0):
    success = False
    error = __name__ + ": PBC are not implemented for CAR"
  
    return success, error
  
  fout.write("%s\n" % "PBC=OFF")
    
  fout.write("\n")
  fout.write("%s\n" % "!DATE")
  
  for i in range(cluster.NAtoms):
    
    tempStr =  ("%s%d" % (cluster.specieList[cluster.specie[i]], i+1))
    tempStr = "{:<7}".format(tempStr)
    
    fout.write("%7s %13.10f %13.10f %13.10f XXXX 1      xx      %2s %.4f\n" % (tempStr, 
      cluster.pos[3*i], cluster.pos[3*i+1], cluster.pos[3*i+2], 
      "{:<2}".format(cluster.specieList[cluster.specie[i]]), cluster.charge[i]))
  
  fout.write("end\n")
  fout.write("end\n")
  
  fout.close()
  
  return success, error

def writeGIN(cluster, outputFile):
  """
  Writes cluster as a GIN file.
  
  """
  
  error = ""
  success = True
  
  masterGinFile = "Master.gin"
  
  if cluster is None:
    success = False
    error = __name__ + ": no data to write"
    
    return success, error
  
  if (not os.path.isfile(masterGinFile)):
    success = False
    error = __name__ + ": could not locate Master.gin file"

    return success, error
    
  masterCluster, error, header, footer = readClusterFromFileGIN(masterGinFile, outputMode=True)
  
  if (masterCluster is None):
    success = False
    return success, error
  
  try:
    fout = open(outputFile, "w")
  except:
    success = False
    error = __name__ + ": Cannot open: " + filePath
     
    return success, error
  
  fout.write(header)

  for i in range(cluster.NAtoms):
    tempStr = cluster.specieList[cluster.specie[i]]
    
    gulpSpecies = masterCluster.gulpSpecies[tempStr]
    
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
                  (tempStr, type, cluster.pos[3*i], cluster.pos[3*i+1], cluster.pos[3*i+2], charge)) 
     
  fout.write(footer)
  
  fout.close()
  
  return success, error

def writeXYZ(cluster, outputFile):
  """
  Writes cluster as an XYZ file.
  
  """
  
  error = ""
  success = True
  
  if cluster is None:
    success = False
    error = __name__ + ": no data to write"
    
    return success, error
  
  try:
    fout = open(outputFile, "w")
    
  except:
    success = False
    error = __name__ + ": Cannot open: " + filePath
    
    return success, error
  
  fout.write("%d\n" % cluster.NAtoms)
  
  metaData = "SCF Done             %.10e;" % (cluster.totalEnergy)
    
  fout.write("%s\n" % (metaData))
  for i in range(cluster.NAtoms):
    
    fout.write("%s %.10f %.10f %.10f %.2f\n" % (cluster.specieList[cluster.specie[i]], 
      cluster.pos[3*i], cluster.pos[3*i+1], cluster.pos[3*i+2], cluster.charge[i]))
  
  fout.close()
  
  return success, error

def convertFile(cluster, outFile):
  """
  Saves a cluster in to a file according to the file format
  
  """
  error = ""
  success = True
  
  if outFile.endswith(".xyz"):
    success, error = writeXYZ(cluster, outFile)
    
  elif outFile.endswith(".car"):
    success, error = writeCAR(cluster, outFile)
    
  elif outFile.endswith(".gin"):
    success, error = writeGIN(cluster, outFile)
    
  else:
    error = "Undefined output format"
    success = False
  
  return success, error
    
if __name__ == "__main__":
  
  ok = 1
  error = ""
  
  _, args = cmdLineArgs()
  
  inFileName = args[0]
  outFileName = args[1]
  
  if inFileName[-3:] == outFileName[-3:]:
    sys.exit("Formats must differ!")
  
  if inFileName == ".gin":
    fileList = glob.glob("*.gin")
      
    for fileName in fileList:
      cluster = readClusterFromFileGIN(fileName)
        
      outputFile = fileName[:-3] + outFileName[-3:]
      
      ok, error = convertFile(cluster, outputFile)
  
  elif inFileName == ".car":
    fileList = glob.glob("*.car")
    
    for fileName in fileList:
      
      cluster = readClusterFromFileCAR(fileName)
       
      outputFile = fileName[:-3] + outFileName[-3:]
      
      ok, error = convertFile(cluster, outputFile)
  
  elif inFileName == ".xyz":
    fileList = glob.glob("*.xyz")
    
    for fileName in fileList:
      cluster = readClusterFromFileXYZ(fileName)
        
      outputFile = fileName[:-3] + outFileName[-3:]
      
      ok, error = convertFile(cluster, outputFile)
    
  else:
    
    if inFileName.endswith(".xyz"):
      cluster = readClusterFromFileXYZ(inFileName)
      
    elif inFileName.endswith(".car"):
      cluster = readClusterFromFileCAR(inFileName)
    
    elif inFileName.endswith(".gin"):
      cluster, error = readClusterFromFileGIN(inFileName)
    
    else:
      print "Unrecognised input file ", inFileName, " format"
      sys.exit()
      
    ok, error = convertFile(cluster, outFileName)
    
    if ok:
      print "Finished!"
    
    else:
      print "Error: ", error
      