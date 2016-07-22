#!/usr/bin/env python

"""
A script to compare unique structures in terms of their energy ranking between different levels of theory (IPs (GULP) and DFT (FHIaims))

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines

import source.IO as IO

inputTypeKLMC_GA = "KLMC_GA"
inputTypeKLMC_DM_TOP = "KLMC_DM_TOP"
klmcDMEnergiesFile = "energies"
klmcDMHashkeysFile = "hashkeys"

# Input 
theoryTypeCnt = 5

theory1Label = 'IP (cores)'
theory1Path = '/Volumes/DATA/ZnO/90_GA_10x_Results/n12/run/run/200/'
theory1Limit = 50
theory1Type = inputTypeKLMC_GA
theory1Format = ""

theory2Label = 'IP (shells)'
theory2Path = '/Volumes/DATA/ZnO/91_Reoptimise_with_shells/n12/top_structures/'
theory2Limit = 50
theory2Type = inputTypeKLMC_DM_TOP
theory2Format = ""

def plot(structCnt1, energyArr1, hashkeys1, label1, 
         structCnt2, energyArr2, hashkeys2, label2):
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

def readDMTopFiles(energiesFile, hashkeysFile, structureLimit):
  """
  Extracts energies and hashkeys from DM simulations
  
  """
  
  success = True
  structCnt = 0
  lineCnt = 0
  
  energyArr = np.zeros(structureLimit, np.float64)
  hashkeys = []
  
  # Reading the energies file
  
  try:
    f = open(energiesFile, "r")
  
  except:
    success = False
    return success, structCnt, energyArr, hashkeys
  
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
    return success, structCnt, energyArr, hashkeys
  
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
  
  return success, structCnt, energyArr, hashkeys

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

  try:
    f = open(filePath, "r")
  
  except:
    success = False
    return success, structCnt, energyArr, hashkeys

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
        
        structCnt += 1
        
    lineCnt += 1

  f.close()
  
  return success, structCnt, energyArr, hashkeys

def readStructuresData(dirPath, limit, type, format):
  
  """
  
  
  """
  
  structures = None
  success = False
  error = ""
  
  cwd = os.getcwd()
  
  if not IO.checkDirectory(dirPath):
    error = "Directory does not exist: %s" % (dirPath)
    return success, error, structures
  
  os.chdir(dirPath)
  
  if (type == inputTypeKLMC_GA):
    
    # look for csv file with the statistics
    csvFile = IO.lookForFiles("csv")
    
    if ((csvFile is None) or (not IO.checkFile(csvFile))):
      error = "File does not exist: %s" % (csvFile)
      return success, error, structures
      
    success, structCnt, energyArr, hashkeys = readGAStatsFile(csvFile, limit)
  
  elif (type == inputTypeKLMC_DM_TOP):
    
    # look for energies and hashkeys files
    if ((not IO.checkFile(klmcDMEnergiesFile)) or (not IO.checkFile(klmcDMHashkeysFile))):
      error = "File(s) does not exist: %s %s" % (klmcDMEnergiesFile, klmcDMHashkeysFile)
      return success, error, structures
    
    success, structCnt, energyArr, hashkeys = readDMTopFiles(klmcDMEnergiesFile, klmcDMHashkeysFile, limit)
         
  else:
    error = "Type '%s' cannot be used at the moment " % (type)
  
  os.chdir(cwd)
  
  return success, error, structCnt, energyArr, hashkeys
  
if __name__ == "__main__":
  
  # reading the data1
  success, error, structCnt1, energyArr1, hashkeys1 = readStructuresData(theory1Path, theory1Limit, theory1Type, theory1Format)
  
  if not success:
    print error
  
  # reading the data2
  success, error, structCnt2, energyArr2, hashkeys2 = readStructuresData(theory2Path, theory2Limit, theory2Type, theory2Format)
  
  if not success:
    print error
  
  # plotting all the data
  plot(structCnt1, energyArr1, hashkeys1, theory1Label, structCnt2, energyArr2, hashkeys2, theory2Label)
  
  print "Done!"