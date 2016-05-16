"""
Finds parents, grandparents, etc for a specific GA iteration.

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import os
import subprocess
import sys
from optparse import OptionParser

# 
GARecursiveDepth = 10

# Settings below should not be changed if you are not sure what they are for
#---------------------------------------------------------------------------
_analyzeKLMCDir = "analysis"
_runKLMCDir = "run"
_popKLMCDir = "POP"

_statsPrefix = "gaStatistics"
_statsSubfix = ".csv"
_tarSubfix = ".tar.gz"

def cmdLineArgs():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog"
  
  parser = OptionParser(usage=usage)

  parser.disable_interspersed_args()
  
  parser.add_option("-s", "--step", dest="step", default=None, type="int", help="")
  parser.add_option("-d", "--depth", dest="depth", default=None, type="int", help="")
  
  (options, args) = parser.parse_args()

  if (len(args) != 0):
    parser.error("No arguments required")

  return options, args
    
def analyze(latestGAStep, initiateMode=False, depth=0, lookFor=None, lookForString=None):
  """
  Performs the family tree search.
  
  """
  
  continueAnalysing = True
  
  _nextGAStep = latestGAStep
  _depth = depth+1
  _initiateMode = initiateMode
  _lookFor = lookFor
  _lookForString = lookForString
  
  if _depth > GARecursiveDepth:
    return
  
  while True:
    found = False
    
    if _initiateMode:
      _initiateMode = False
      found, rank, hashkey, energy, origin, gaGenNo, parent1, parent2 = findLowest(_nextGAStep)
         
    else:
      # looking for a hashkey in the statistics file
      if _lookFor == "#":
        found, rank, hashkey, energy, origin, gaGenNo, parent1, parent2 = findHashkey(_nextGAStep, _lookForString)
      
      # looking for file in POP folder
      elif _lookFor == "f":
        found, rank, hashkey, energy, origin, gaGenNo, parent1, parent2 = findFile(_nextGAStep, _lookForString)
        
      else:
        print "Analysing: %d"% (_nextGAStep), "could not understand what should I look for"
      
      _lookFor = None
      _lookForString = None
    
    if not found:
      print "Analysing: %d"% (_nextGAStep), "could not locate next step"
      break
    
    else:
      print "Analysing: %d"% (_nextGAStep), depth, found, rank, hashkey, energy, origin, gaGenNo, parent1, parent2
    
    if (origin == "REPOPM"):
      if gaGenNo < _nextGAStep:
        _nextGAStep = gaGenNo
      else:
        _nextGAStep -= 1
        
      _lookFor = "#"
      _lookForString = hashkey
    
    elif ((origin == "MUTATE") or (origin == "CROSSO")):
      if gaGenNo < _nextGAStep:
        _nextGAStep = gaGenNo - 1
      else:
        _nextGAStep -= 1
        
      _lookFor = "f"
      
      _lookForString = parent1
      analyze(_nextGAStep, depth=_depth, lookFor=_lookFor, lookForString=_lookForString)
      
      _lookForString = parent2
      analyze(_nextGAStep, depth=_depth, lookFor=_lookFor, lookForString=_lookForString)
      
      break
    
    elif (origin == "MUTCRS"):
      if gaGenNo < _nextGAStep:
        _nextGAStep = gaGenNo
      else:
        _nextGAStep -= 1
        
      _lookFor = "#"
      _lookForString = hashkey
      
    else:
      print "Analysing: %d"% (_nextGAStep), "could not understand the origin of the structure"
    
def findFile(stepNo, fileName):
  """
  Looks for a specific file and returns its info from the statistics file.
  
  """
  
  found = False
  rank = -1
  hashkey = ""
  origin  = ""
  gaGenNo = -1
  parent1 = ""
  parent2 = ""
  
  # unzip the tar file
  cmdLine = "gzip -d < %s | tar xf -" % ("%d%s" % (stepNo, _tarSubfix))
  os.system(cmdLine)
  
  cwd = os.getcwd()
  os.chdir(os.path.join(_runKLMCDir, "%d" % (stepNo), _popKLMCDir))
  
  cmdLine = "ls *%s_1.xyz" % (fileName.strip())  
  output, stderr, status = runSubProcess(cmdLine)
  
  os.chdir(cwd)
  
  try:
    rank = int(output.strip().split("_")[0])
    found = True
    
  except:
    rank = -1
  
  # lets read the data from the statistics file
  if found:    
    filePath = os.path.join(_runKLMCDir, "%d" % (stepNo), "%s%d%s" % (_statsPrefix, stepNo, _statsSubfix))
    
    found, hashkey, energy, origin, gaGenNo, parent1, parent2 = readGAStatsFileRank(filePath, rank)
      
  # cleaning up
  cmdLine = "rm -rf run"
  os.system(cmdLine)
  
  return found, rank, hashkey, energy, origin, gaGenNo, parent1, parent2
    
def findHashkey(stepNo, hashkey):
  """
  Looks for a specific hashkey and returns its info from the statistics file.
  
  """
  
  # unzip the tar file
  cmdLine = "gzip -d < %s | tar xf -" % ("%d%s" % (stepNo, _tarSubfix))
  os.system(cmdLine)
  
  # reading the statistics file
  filePath = os.path.join(_runKLMCDir, "%d" % (stepNo), "%s%d%s" % (_statsPrefix, stepNo, _statsSubfix))
  found, rank, energy, origin, gaGenNo, parent1, parent2 = readGAStatsFileHashkey(filePath, hashkey)
  
  # cleaning up
  cmdLine = "rm -rf run"
  os.system(cmdLine)
  
  return found, rank, hashkey, energy, origin, gaGenNo, parent1, parent2
  
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
  
def getLatestGaStep(dirPath=None):
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

def readGAStatsFileRank(filePath, rank):
  """
  Returns statistics from the statistics file with respect to the rank.
  
  """
    
  found = False
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

    if lineCnt == rank:
      found = True
      
      line = line.strip()
      array = line.split(",")

      hashkey = array[2]
      energy = float(array[5])
      parent1 = array[8]
      parent2 = array[9]
      origin = array[14]
      gaGenNo = int(array[15])

      break
      
    lineCnt += 1
    
  f.close()
  
  return found, hashkey, energy, origin, gaGenNo, parent1, parent2

def readGAStatsFileHashkey(filePath, hashkey):
  """
  Returns statistics from the statistics file with respect to the hashkey.
  
  """
  
  found = False
  rank = 0
  hashhey = ""
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
    if lineCnt > 0:
      line = line.strip()
      array = line.split(",")
            
      if array[2] == hashkey:
        found = True
        rank = lineCnt 
        energy = float(array[5])
        parent1 = array[8]
        parent2 = array[9]
        origin = array[14]
        gaGenNo = int(array[15])

        break
      
    lineCnt += 1
    
  f.close()
  
  return found, rank, energy, origin, gaGenNo, parent1, parent2

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

def runSubProcess(command, verbose=0):
  """
  Runs command using subprocess module.
    
  """
  
  if verbose:
      print command
  
  process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output, stderr = process.communicate()
  status = process.poll()
  
  return (output, stderr, status)

if __name__ == "__main__":
  
  options, args = cmdLineArgs()
  
  # setting the depth of recursion
  if options.depth is not None:
    GARecursiveDepth = options.depth
    
  if options.step is not None:
    latestGAStep = options.step
  else:
    # finds the latest GA step
    latestGAStep = getLatestGaStep(_runKLMCDir)

  cwd = os.getcwd()
  os.chdir(_runKLMCDir)
  
  # analyze
  analyze(latestGAStep, initiateMode=True)
  
  os.chdir(cwd)