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

_inputTypeKLMC_GA = "KLMC_GA"
_inputTypeKLMC_DM = "KLMC_DM"
_inputTypeJobArr = "JOBARR_DM"

_klmcDMEnergiesFile = "energies"
_klmcDMHashkeysFile = "hashkeys"
_klmcDMStatsFile = "statistics"

_klmcRunDir = "run/"
_klmcTopDir = "top_structures/"
_klmcProdStatsFile = "prodStatistics.csv"
_klmcOutFile = "KLMC.out"
_klmcLogFile = "KLMC.log"

_jobArrOutDir = "output/"
_jobArrUniqueDir = "output/unique/"
_jobArrStatsFile = "Stats.csv"

# TODO: Maybe this should be defined in a more general way, input file?
# Input 
theoryTypeCnt = 4

theory1Label = 'IP (cores)'
theory1Path = '/Volumes/DATA/ZnO/90_GA_10x_Results/n12/run/run/200/'
theory1Limit = 50
theory1Type = _inputTypeKLMC_GA
theory1Format = ""

theory2Label = 'IP (shells)'
theory2Path = '/Volumes/DATA/ZnO/91_Reoptimise_with_shells/n12/'
theory2Limit = 50
theory2Type = _inputTypeKLMC_DM
theory2Format = ""

theory3Label = 'PBESol (light)'
theory3Path = '/Volumes/DATA/ZnO/94_Aims_Archer/1_Light_12_32/2_Sims/n12/'
theory3Limit = 50
theory3Type = _inputTypeJobArr
theory3Format = ""

theory4Label = 'PBESol (tight)'
theory4Path = '/Volumes/DATA/ZnO/94_Aims_Archer/2_Tight_12_32/2_Sims/n12/'
theory4Limit = 50
theory4Type = _inputTypeJobArr
theory4Format = ""

class SimulationResult(object):
  """
  A class to save the simulation data
    
  """
  
  success   = False
  error     = None
  simType   = None
  structCnt = None
  energyArr = None
  hashkeys  = None
  prevRank  = None
  currRank  = None
  currRankHash = None
  
  def __init__(self, simType, structCntLimit):
    """
    Constructor
    
    """
    
    self.simType = simType
    
    self.structCnt = 0
    
    self.success   = False
    self.error     = None
    
    self.energyArr = np.zeros(structCntLimit, np.float64)
    self.hashkeys  = []
    self.prevRank  = np.zeros(structCntLimit, np.int16)
    self.currRank  = np.zeros(structCntLimit, np.int16)
    self.currRankHash = np.zeros(structCntLimit, np.int16)
  
  def appendStructure(self, energy, hashkey, prevRankEl, currRankEl, hashkeyRankEl):
    """
    Appends a structure to the results 
    
    """

    self.energyArr[self.structCnt] = energy
    self.hashkeys.append(hashkey)      
    self.prevRank[self.structCnt] = int(prevRankEl)
    self.currRank[self.structCnt] = int(currRankEl)
    self.currRankHash[self.structCnt] = int(hashkeyRankEl)
    
    self.structCnt += 1
    
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

def lookForStructureRankJobArr(filePath, structureName, structureHashkey):
  """
  Scans through a top structures file to find the rank of a structure
  
  """
  
  rank = -1
  lineCnt = 0
  lookForHashkey = False
  
  # TODO: Names can also start with B, C, D, E!
  if (structureName is not None):
    structureNameReplaced = structureName.replace("A", "X")
  
  if (structureHashkey is not None):
    lookForHashkey = True
  
  try:
    f = open(filePath, "r")
    
  except:
    return rank
  
  for line in f:
    if lineCnt > 0:
      line = line.strip()

      array = line.split(",")
      name = array[0].strip()
      hashkey = array[2].strip()
          
      if lookForHashkey:
        if hashkey == structureHashkey:
          rank = lineCnt
          break
        
      elif name == structureName:
        rank = lineCnt
        break
      
    lineCnt += 1
  
  return rank

def plot(mSimulationResult1, label1, mSimulationResult2, label2, mSimulationResult3, label3, 
         mSimulationResult4, label4):
  """
  Main plotting routine.
  
  TODO: Needs to be rewritten in a more general way
  
  """
  
  plt.xkcd()

  fig, axes = plt.subplots(nrows=1, ncols=4, sharex=False, sharey=False)
  
  # Plot 1
  x1 = np.zeros(mSimulationResult1.structCnt, np.int16)
  y1 = mSimulationResult1.energyArr[:mSimulationResult1.structCnt]
  axes[0].scatter(x1, y1, s=500)
  
  my_xticks = [label1]
  plt.sca(axes[0])
  plt.xticks(x1, my_xticks)
  
  # Plot 2
  x2 = np.zeros(mSimulationResult2.structCnt, np.int16)
  y2 = mSimulationResult2.energyArr[:mSimulationResult2.structCnt]
  axes[1].scatter(x2, y2, s=500)
  
  my_xticks = [label2]
  plt.sca(axes[1])
  plt.xticks(x2, my_xticks)
  
  # Plot 3
  x3 = np.zeros(mSimulationResult3.structCnt, np.int16)
  y3 = mSimulationResult3.energyArr[:mSimulationResult3.structCnt]
  
  axes[2].scatter(x3, y3, s=500)
  axes[2].ticklabel_format(style='plain')
  
  my_xticks = [label3]
  plt.sca(axes[2])
  plt.xticks(x3, my_xticks)
  
  # Plot 4
  x4 = np.zeros(mSimulationResult4.structCnt, np.int16)
  y4 = mSimulationResult4.energyArr[:mSimulationResult4.structCnt]
  
  axes[3].scatter(x4, y4, s=500)
  axes[3].ticklabel_format(style='plain')
  
  my_xticks = [label4]
  plt.sca(axes[3])
  plt.xticks(x4, my_xticks)
    
  # Joining plot 1 - plot 2
  transFigure = fig.transFigure.inverted()
  
  for i in range(mSimulationResult1.structCnt):
    
    currRank1 = mSimulationResult1.currRank[i]
    nextRank1 = -1
    
    for j in range(mSimulationResult2.structCnt):
      
      prevRank2 = mSimulationResult2.prevRank[j]
      
      if prevRank2 == currRank1:
        
        currRank2 = mSimulationResult2.currRank[j]
        currHashRank2 = mSimulationResult2.currRankHash[j]
        
        if currRank2 > 0:
          nextRank1 = currRank2
        else:
          nextRank1 = currHashRank2
        break
      
    if currRank1 > 0 and nextRank1 > 0:
      
      currPos1 = currRank1 - 1
      currPos2 = nextRank1 - 1
      
      coord1 = transFigure.transform(axes[0].transData.transform([x1[currPos1], y1[currPos1]]))
       
      coord2 = transFigure.transform(axes[1].transData.transform([x2[currPos1], y2[currPos1]]))
 
      line = lines.Line2D((coord1[0], coord2[0]),(coord1[1], coord2[1]), transform=fig.transFigure)
           
      fig.lines.append(line)
    
  # Joining plot 2 - plot 3  
  for i in range(mSimulationResult3.structCnt):
    
    prevRank3 = mSimulationResult3.prevRank[i]
    currRank3 = mSimulationResult3.currRank[i]
    hashRank3 = mSimulationResult3.currRankHash[i]
    
    if currRank2 != -1:
      rank3 = currRank3
    else:
      rank3 = hashRank3
    
    prevPos = prevRank3-1
    
    coord2 = transFigure.transform(axes[1].transData.transform([x2[prevPos], y2[prevPos]]))
      
    coord3 = transFigure.transform(axes[2].transData.transform([x3[i], y3[i]]))

    line = lines.Line2D((coord2[0], coord3[0]),(coord2[1], coord3[1]), transform=fig.transFigure)
          
    fig.lines.append(line)
    
  # Joining plot 3 - plot 4
  for i in range(mSimulationResult4.structCnt):
    
    prevRank4 = mSimulationResult4.prevRank[i]
    currRank4 = mSimulationResult4.currRank[i]
    hashRank4 = mSimulationResult4.currRankHash[i]
     
    prevIniRank4 = -1
        
    for j in range(mSimulationResult3.structCnt):
      
      if mSimulationResult3.currRank[j] != -1:
        currRank3 = mSimulationResult3.currRank[j]
        
      else:
        currRank3 = mSimulationResult3.currRankHash[j]
      
      if currRank3 == prevRank4:

        prevIniRank4 = mSimulationResult3.prevRank[j]
        break
        
    if currRank4 != -1:
      rank4 = currRank4
      
    else:
      rank4 = hashRank4
    
    prevPos = prevIniRank4-1
    currPos = rank4-1
        
    coord3 = transFigure.transform(axes[2].transData.transform([x3[j], y3[j]]))
       
    coord4 = transFigure.transform(axes[3].transData.transform([x4[currPos], y4[currPos]]))
 
    line = lines.Line2D((coord3[0], coord4[0]), (coord3[1], coord4[1]), transform=fig.transFigure)
           
    fig.lines.append(line)
    
  fig.savefig('Comparison.png')


def readDMFiles(structureLimit):
  """
  Extracts energies and hashkeys from DM simulations
  
  """
  
  mSimulationResult = SimulationResult(_inputTypeKLMC_DM, structureLimit)
  
  structCnt = 0
  lineCnt = 0
    
  # reading the statistics file
  try:
    f = open("%s/%s" % (_klmcRunDir, _klmcProdStatsFile), "r")
   
  except:
    return mSimulationResult
  
  for line in f:
    line = line.strip()
    
    if lineCnt > 0:
      array = line.split(",")
      
      # Extracting the data
      name = array[1].strip()
      hashkey = array[2].strip()
      edfn = int(array[3])
      energy = float(array[4])
      
      # Getting the ranking
      prevRankEl = lookForIniStructureKLMC(_klmcLogFile, name)
      
      currRankEl = lookForStructureRankKLMC("%s/%s" % (_klmcTopDir, _klmcDMStatsFile), name, None)
      
      hashkeyRankEl = lookForStructureRankKLMC("%s/%s" % (_klmcTopDir, _klmcDMHashkeysFile), None, hashkey)
            
      if edfn > 0 and energy != 0.0 and structCnt < structureLimit:
                
        # Appending the result
        mSimulationResult.appendStructure(energy, hashkey, prevRankEl, currRankEl, hashkeyRankEl)
          
        structCnt += 1

    lineCnt += 1
      
  f.close()
      
  return mSimulationResult

def readGAStatsFile(filePath, structureLimit):
  """
  Extracts energies and hashkeys from the GA statistics file assuming that the entries are ordered 
  with respect to the energy.
  
  """
  structCnt = 0
  lineCnt = 0
  
  mSimulationResult = SimulationResult(_inputTypeKLMC_GA, structureLimit)
    
  try:
    f = open(filePath, "r")
  
  except:
    return mSimulationResult

  for line in f:

    if lineCnt > 0:
      line = line.strip()
      array = line.split(",")
      
      hashkey = array[2].strip()
      edfn = int(array[3])
      status = int(array[4])
      energy = float(array[5])
        
      if edfn > 0 and status == 1 and energy != 0.0 and structCnt < structureLimit:
        structCnt += 1
         
        prevRankEl = -1
        currRankEl = structCnt
        hashkeyRankEl = -1
        
        mSimulationResult.appendStructure(energy, hashkey, prevRankEl, currRankEl, hashkeyRankEl)
      
        
        
    lineCnt += 1

  f.close()
  
  return mSimulationResult

def readJobArrFiles(structureLimit):
  """
  Extracts energies and hashkeys from job array simulations
  
  """
  
  mSimulationResult = SimulationResult(_inputTypeKLMC_DM, structureLimit)
  
  structCnt = 0
  lineCnt = 0
    
  # reading the statistics file
  try:    
    f = open("%s/%s" % (_jobArrOutDir, _jobArrStatsFile), "r")
   
  except:
    return mSimulationResult
    
  for line in f:
    line = line.strip()
        
    if lineCnt > 0:
      array = line.split(",")
             
      # Extracting the data
      name = array[0].strip()
      energy = float(array[1])
      hashkey = array[2].strip()
      
      # extracting previous rank from name 
      nameSplit = name.split("_")
            
      if (len(nameSplit) > 1):
        
        try:
          prevRankEl = int(nameSplit[1])
        except:
          prevRankEl = -1
        
      else:
        try:
          prevRankEl = int(name[1:])
        except:
          prevRankEl = -1
      
      currRankEl = lookForStructureRankJobArr("%s/%s" % (_jobArrUniqueDir, _jobArrStatsFile), name, None)
      hashkeyRankEl = lookForStructureRankJobArr("%s/%s" % (_jobArrUniqueDir, _jobArrStatsFile), name, hashkey)
             
      if energy != 0.0 and structCnt < structureLimit:       
        # Appending the result
        mSimulationResult.appendStructure(energy, hashkey, prevRankEl, currRankEl, hashkeyRankEl)
            
        structCnt += 1
            
    lineCnt += 1
      
  f.close()
    
  return mSimulationResult

def readStructuresData(dirPath, limit, type, format):
  """
  Depending on the input type, reads in the structural data
  
  """
  
  success = False
  error = ""

  mSimulationResult = None
  
  cwd = os.getcwd()
  
  if not IO.checkDirectory(dirPath):
    error = "Directory does not exist: %s" % (dirPath)
    
    return success, error, mSimulationResult
  
  os.chdir(dirPath)
  
  if (type == _inputTypeKLMC_GA):
    
    # look for csv file with the statistics
    csvFile = IO.lookForFiles("csv")
    
    if ((csvFile is None) or (not IO.checkFile(csvFile))):
      error = "File does not exist: %s" % (csvFile)
      
      return success, error, mSimulationResult
      
    mSimulationResult = readGAStatsFile(csvFile, limit)
  
  elif (type == _inputTypeKLMC_DM):
    
    mSimulationResult = readDMFiles(limit)
  
  elif (type == _inputTypeJobArr):
  
    mSimulationResult = readJobArrFiles(limit)
    
  else:
    error = "Type '%s' cannot be used at the moment " % (type)
  
  os.chdir(cwd)
  
  return success, error, mSimulationResult
  
if __name__ == "__main__":
  
  # reading the data1
  success, error, mSimulationResult1 = readStructuresData(theory1Path, theory1Limit, theory1Type, theory1Format)
  
  if not success:
    print error
  
  # reading the data2
  success, error, mSimulationResult2 = readStructuresData(theory2Path, theory2Limit, theory2Type, theory2Format)
  
  if not success:
    print error
    
  # reading the data3
  success, error, mSimulationResult3 = readStructuresData(theory3Path, theory3Limit, theory3Type, theory3Format)
  
  if not success:
    print error
    
  # reading the data4
  success, error, mSimulationResult4 = readStructuresData(theory4Path, theory4Limit, theory4Type, theory4Format)
    
  if not success:
    print error
  
  # plotting all the data
  plot(mSimulationResult1, theory1Label, mSimulationResult2, theory2Label, mSimulationResult3, theory3Label,
       mSimulationResult4, theory4Label)
  
  print "Done!"