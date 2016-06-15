"""
A script to find the GA iteration number on which the GM has been found

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import os, sys

_runKLMCDir = "run"

_statsPrefix = "gaStatistics"
_statsSubfix = ".csv"
_tarSubfix = ".tar.gz"

def findLowest(stepNo):
  """
  Looks for the lowest energy structure and returns its info from the statistics file.
  
  """
  
  # unzip the tar file
  cmdLine = "gzip -d < %s | tar xf -" % ("%d%s" % (stepNo, _tarSubfix))
  os.system(cmdLine)
  
  # reading the statistics file
  filePath = os.path.join(_runKLMCDir, "%d" % (stepNo), "%s%d%s" % (_statsPrefix, stepNo, _statsSubfix))
  found, rank, hashkey, energy, origin, gaGenNo, parent1, parent2 = readGAStatsFileTop(filePath)
  
  # cleaning up
  cmdLine = "rm -rf run"
  os.system(cmdLine)
  
  return found, rank, hashkey, energy, origin, gaGenNo, parent1, parent2

def getLatestGAStep(dirPath=None):
  """
  Gets the last GA iteration number
  
  """
  
  filesList = []
  gaItNoList = []
  
  if (dirPath is not None):
    cwd = os.getcwd()
    os.chdir(dirPath)
    
  else:
    cwd = None
  
  maxGAIter = -1
  
  subfixLen = len(_tarSubfix)

  for root, _, files in os.walk("./"):
    for fileName in files:
      if fileName.endswith(_tarSubfix):
        fileNameLen = len(fileName)
        gaItNo = int(fileName[:fileNameLen-subfixLen])
        
        if (gaItNo > maxGAIter):
          maxGAIter = gaItNo  
  
  if (dirPath is not None):
    os.chdir(cwd)
    
  return maxGAIter

def readGAStatsFileTop(filePath):
  """
  Reads in the GA iteration statistics file.
  
  """
  
  found = False
  rank = 0
  hashkey = ""
  energy = 0.0
  origin  = ""
  gaGenNo = -1
  parent1 = ""
  parent2 = ""
  
  lineCnt = 0
  
  try:
    f = open(filePath, "r")
    
  except:
    sys.exit("Error: cannot open: " + filePath)
  
  for line in f:
    if lineCnt == 1:
      line = line.strip()
      array = line.split(",")
      
      found = True
      rank = lineCnt
      hashkey = array[2]
      energy = float(array[5])
      parent1 = array[8]
      parent2 = array[9]
      origin = array[14]
      gaGenNo = int(array[15])

      break
    
    lineCnt += 1
    
  f.close()
  
  return found, rank, hashkey, energy, origin, gaGenNo, parent1, parent2

if __name__ == "__main__":

  lastGA = getLatestGAStep(_runKLMCDir)
  
  cwd = os.getcwd()
  os.chdir(_runKLMCDir)
  
  prevIter = lastGA
  
  _, _, masterHashkey, prevEnergy, _, prevGAGenNo, _, _ = findLowest(prevIter)

  iter = prevIter/2
  
  continueSearch = True
  while continueSearch:
    
    _, _, hashkey, energy, _, GAGenNo, _, _ = findLowest(iter)
        
    checkIterationCpy = iter
    
    if hashkey == masterHashkey:
      iter = iter / 2
      prevIter = checkIterationCpy 
      
    else:
      iter = (iter + prevIter) / 2
    
    if (prevIter - iter) < 2:
            
      continueSearch = False
  
  _, _, prevHashkey, prevEnergy, _, prevGAGenNo, _, _ = findLowest(prevIter)
  
  print prevIter, prevHashkey, prevEnergy, prevGAGenNo
  
  os.chdir(cwd)
