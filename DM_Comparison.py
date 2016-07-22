#!/usr/bin/env python

"""
A script to compare unique structures in terms of their energy ranking between different levels of theory (IPs (GULP) and DFT (FHIaims))

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines

import source.IO as IO

inputTypeKLMC_GA = "KLMC_GA"
inputTypeKLMC_DM = "KLMC_DM"
inputTypeKLMC_DM_TOP = "KLMC_DM_TOP"

klmcDMEnergiesFile = "energies"
klmcDMHashkeysFile = "hashkeys"
klmcDMStatsFile = "statistics"

klmcRunDir = "run/"
klmcTopDir = "top_structures/"
klmcProdStatsFile = "prodStatistics.csv"
klmcOutFile = "KLMC.out"
klmcLogFile = "KLMC.log"

# Input 
theoryTypeCnt = 5

theory1Label = 'IP (cores)'
theory1Path = '/Volumes/DATA/ZnO/90_GA_10x_Results/n12/run/run/200/'
theory1Limit = 50
theory1Type = inputTypeKLMC_GA
theory1Format = ""

theory2Label = 'IP (shells)'
theory2Path = '/Volumes/DATA/ZnO/91_Reoptimise_with_shells/n12/'
theory2Limit = 50
theory2Type = inputTypeKLMC_DM
theory2Format = ""

def lookForIniStructureKLMC(klmcOut, structureName):
  """
  Scans klmc outfile, expecting it to be in the DM mode and tries to locate the initial file name
  
  """
  
  prevRank = -1
  
  # TODO: Names can also start with B, C, D, E!
  structureNameReplaced = structureName.replace("A", "X")
  
  try:
    f = open(klmcOut, "r")
   
  except:
    return prevRank
  
  strExpr = " as %s " % (structureNameReplaced)

  for line in f: 
    if strExpr in line:
      startExpr = "from restart/"
      
      fileNameStart = line.find(startExpr)
      fileNameEnd = line.find(strExpr)
      
      fileName = line[fileNameStart+len(startExpr):fileNameEnd]
      
      # TODO: This might not be allways the case!
      prevRank = int(fileName.split("_")[0])

  f.close()
  
  return prevRank

def lookForStructureRankKLMC(filePath, structureName, structureHashkey):
  """
  Scans through a statistics file to find the rank of a structure
  
  """
  
  rank = -1
  lineCnt = 0
  
  # TODO: Names can also start with B, C, D, E!
  if (structureName is not None):
    structureNameReplaced = structureName.replace("A", "X")
  
  try:
    f = open(filePath, "r")
  except:
    return rank
  
  for line in f:
    if lineCnt > 0:
      line = line.strip()

      array = line.split()
      name = array[1].strip()
      rankFile = int(array[0].strip())
      
      if (structureHashkey is not None):
        hashkey = array[2].strip()
      else:
        hashkey = None
      
      if ((structureName is not None) and (name == structureName)):
        rank = rankFile
        return rank
      
      elif ((structureHashkey is not None) and (hashkey == structureHashkey)):
        rank = rankFile
        return rank
    
    lineCnt += 1
  
  return rank
  
def plot(structCnt1, energyArr1, hashkeys1, label1, prevRank1, currRank1, 
         structCnt2, energyArr2, hashkeys2, label2, prevRank2, currRank2):
  """
  
  """
  
  plt.xkcd()

  fig, axes = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=False)
  
  # Plot 1
  x1 = np.zeros(structCnt1, np.int16)
  y1 = energyArr1[:structCnt1]
  axes[0].scatter(x1, y1, s=500)
  
  my_xticks = [label1]
  plt.sca(axes[0])
  plt.xticks(x1, my_xticks)
  
  # Plot 2
  x2 = np.zeros(structCnt2, np.int16)
  y2 = energyArr2[:structCnt2]
  
  axes[1].scatter(x2, y2, s=500)
  
  my_xticks = [label2]
  plt.sca(axes[1])
  plt.xticks(x2, my_xticks)
  
  # Joining plot 1 - plot 2
  transFigure = fig.transFigure.inverted()
  
  i = 0
  
  coord1 = transFigure.transform(axes[0].transData.transform([x1[i], y1[i]]))
  coord2 = transFigure.transform(axes[1].transData.transform([x2[i], y2[i]]))
  
  line = lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),
                                 transform=fig.transFigure)
  fig.lines = line,

  
  
  plt.show()

def readDMFiles(structureLimit):
  """
  Extracts energies and hashkeys from DM simulations
  
  """
  
  success = True
  structCnt = 0
  lineCnt = 0
  
  energyArr = np.zeros(structureLimit, np.float64)
  hashkeys = []
  
  prevRank = np.zeros(structureLimit, np.int16)
  currRank = np.zeros(structureLimit, np.int16)
  hashkeyRank = np.zeros(structureLimit, np.int16)
  
  # reading the statistics file
  try:
    f = open("%s/%s" % (klmcRunDir, klmcProdStatsFile), "r")
   
  except:
    success = False
    return success, structCnt, energyArr, hashkeys, prevRank, currRank
  
  for line in f:
    line = line.strip()
    
    if lineCnt > 0:
      array = line.split(",")
      
      name = array[1].strip()
      hashkey = array[2].strip()
      edfn = int(array[3])
      energy = float(array[4])
      
      prevRankEl = lookForIniStructureKLMC(klmcLogFile, name)
      currRankEl = lookForStructureRankKLMC("%s/%s" % (klmcTopDir, klmcDMStatsFile), name, None)
      hashkeyRankEl = lookForStructureRankKLMC("%s/%s" % (klmcTopDir, klmcDMHashkeysFile), None, hashkey)
            
      print name, hashkey, edfn, energy, prevRank, currRank, hashkeyRank
    
      if edfn > 0 and energy != 0.0 and structCnt < structureLimit:
          
        energyArr[structCnt] = energy
        hashkeys.append(hashkey)
          
        prevRank[structCnt] = prevRankEl
        currRank[structCnt] = currRankEl
        hashkeyRank[structCnt] = hashkeyRankEl
          
        structCnt += 1

    lineCnt += 1
      
  f.close()
  
  # reading the KLMC output file to get the names
  
  
  sys.exit("BYE BYE") 
  
  return success, structCnt, energyArr, hashkeys, prevRank, currRank, hashkeyRank

def readDMTopFiles(energiesFile, hashkeysFile, structureLimit):
  """
  Extracts energies and hashkeys from DM simulations
  
  """
  
  success = True
  structCnt = 0
  lineCnt = 0
  
  energyArr = np.zeros(structureLimit, np.float64)
  hashkeys = []
  
  prevRank = np.zeros(structureLimit, np.int16)
  currRank = np.zeros(structureLimit, np.int16)
  
  # Reading the energies file
  
  try:
    f = open(energiesFile, "r")
  
  except:
    success = False
    return success, structCnt, energyArr, hashkeys, prevRank, currRank
  
  for line in f:
    line = line.strip()
    
    array = line.split()
    energy = float(array[1])
    
    if structCnt < structureLimit:
      energyArr[structCnt] = energy
      
      structCnt += 1
      
  f.close()
  
  # Reading the hashkeys file
  
  try:
    f = open(hashkeysFile, "r")
  
  except:
    success = False
    return success, structCnt, energyArr, hashkeys, prevRank, currRank
  
  structCnt2 = 0
  
  for line in f:
    
    if lineCnt > 0:
      line = line.strip()
      
      array = line.split()
      hashkey = array[2]
      
      if ((structCnt2 < structureLimit) and (structCnt2 < structCnt)):
        structCnt2 += 1
      
    lineCnt += 1

  f.close()
  
  return success, structCnt, energyArr, hashkeys, prevRank, currRank

def readGAStatsFile(filePath, structureLimit):
  """
  Extracts energies and hashkeys from the GA statistics file assuming that the entries are ordered 
  with respect to the energy.
  
  """
  
  success = True
  lineCnt = 0
  structCnt = 0
  
  energyArr = np.zeros(structureLimit, np.float64)
  hashkeys = []
  
  prevRank = np.zeros(structureLimit, np.int16)
  currRank = np.zeros(structureLimit, np.int16)
  
  try:
    f = open(filePath, "r")
  
  except:
    success = False
    return success, structCnt, energyArr, hashkeys, prevRank, currRank

  for line in f:

    if lineCnt > 0:
      line = line.strip()
      array = line.split(",")
      
      hashkey = array[2].strip()
      edfn = int(array[3])
      status = int(array[4])
      energy = float(array[5])

      if edfn > 0 and status == 1 and energy != 0.0 and structCnt < structureLimit:
        
        energyArr[structCnt] = energy
        hashkeys.append(hashkey)
        
        prevRank[structCnt] = -1
        currRank[structCnt] = structCnt
        
        structCnt += 1
        
    lineCnt += 1

  f.close()
  
  return success, structCnt, energyArr, hashkeys, prevRank, currRank

def readStructuresData(dirPath, limit, type, format):
  """
  
  """
  
  structures = None
  success = False
  error = ""
  energyArr = None
  hashkeys = None
  prevRank = None
  currRank = None
  
  cwd = os.getcwd()
  
  if not IO.checkDirectory(dirPath):
    error = "Directory does not exist: %s" % (dirPath)
    return success, error, structures, energyArr, hashkeys, prevRank, currRank
  
  os.chdir(dirPath)
  
  if (type == inputTypeKLMC_GA):
    
    # look for csv file with the statistics
    csvFile = IO.lookForFiles("csv")
    
    if ((csvFile is None) or (not IO.checkFile(csvFile))):
      error = "File does not exist: %s" % (csvFile)
      return success, error, structures, energyArr, hashkeys, prevRank, currRank
      
    success, structCnt, energyArr, hashkeys, prevRank, currRank = readGAStatsFile(csvFile, limit)
  
  elif (type == inputTypeKLMC_DM):
    
    # look for energies and hashkeys files
#     if ((not IO.checkFile(klmcDMEnergiesFile)) or (not IO.checkFile(klmcDMHashkeysFile))):
#       error = "File(s) does not exist: %s %s" % (klmcDMEnergiesFile, klmcDMHashkeysFile)
#       return success, error, structures, energyArr, hashkeys, prevRank, currRank
    
    success, structCnt, energyArr, hashkeys, prevRank, currRank = readDMFiles(limit)
         
  else:
    error = "Type '%s' cannot be used at the moment " % (type)
  
  os.chdir(cwd)
  
  return success, error, structCnt, energyArr, hashkeys, prevRank, currRank
  
if __name__ == "__main__":
  
  # reading the data1
  success, error, structCnt1, energyArr1, hashkeys1, prevRank1, currRank1 = readStructuresData(theory1Path, theory1Limit, theory1Type, theory1Format)
  
  if not success:
    print error
  
  # reading the data2
  success, error, structCnt2, energyArr2, hashkeys2, prevRank2, currRank2 = readStructuresData(theory2Path, theory2Limit, theory2Type, theory2Format)
  
  if not success:
    print error
  
  # plotting all the data
 # plot(structCnt1, energyArr1, hashkeys1, theory1Label, prevRank1, currRank1,
 #      structCnt2, energyArr2, hashkeys2, theory2Label, prevRank2, currRank2)
  
  print "Done!"