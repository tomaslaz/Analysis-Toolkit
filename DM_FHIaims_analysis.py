#!/usr/bin/env python

"""
A script to analyse FHI-aims simulations

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

import copy
import os
import sys

import numpy as np

import source.Atoms as Atoms
import source.Fhiaims as FHIaims
import source.IO as IO

_fhiaimsGeometryFile = "geometry.in"
_fhiaimsOutFile = "FHIaims.out"
_outputDir = "output"
_topDir = "tops"

def getFileList(dirPath=None):
  """
  Checks the directory for fhiaims output files with ending "FHIaims.out" and returns a directory list
  
  """
    
  dirList = []

  for root, _, files in os.walk("./"):
    for fileName in files:
      if fileName == _fhiaimsOutFile:
        
        dirList.append(os.path.join(root[2:]))
    
  return dirList

def generateStatistics(systemlist):
  """
  Generates statistics about the FHI-aims simulations
    
  """
  
  noOfSystems = len(systemlist)
  sumOfCores = 0
  runTimes = np.zeros(noOfSystems, np.float64)
  
  f = open("Stats.csv", "w")
  
  f.write("%s,%s,%s,%s,%s,%s\n" % ("System", "Energy", "Hashkey", "Cores", "Time", "Tot.Time"))
  
  systemCnt = 0
  for system in systemlist:
    
    f.write("%s,%f,%s,%d,%f,%f\n" % (system.name, system.totalEnergy, system.hashkey,
                                     system.noOfcores, system.runTime, system.noOfcores*system.runTime))
    
    runTimes[systemCnt] = system.runTime
    
    sumOfCores += system.noOfcores
    systemCnt += 1
    
  f.close()
  
  avgCores = sumOfCores / noOfSystems
  
  print "Avg. number of cores: ", avgCores
  print "Avg. run time: ", np.average(runTimes)
  print "Min run time: ", np.min(runTimes)
  print "Max run time: ", np.max(runTimes)
  print "Stdev run time: ", np.std(runTimes)
  
def readFHIaimsSystems(fhiaimsDirs):
  """
  Reads in the FHI-aims systems.
  
  """
  
  systemsList = []
  
  for dirName in fhiaimsDirs:
    cwd = os.getcwd()
    
    os.chdir(dirName)

    systemName = dirName[len(_outputDir)+1:].strip()
    
    success, error, system = FHIaims._readAimsStructure(_fhiaimsGeometryFile, _fhiaimsOutFile)
    system.name = systemName
    
    systemsList.append(system)
    
    os.chdir(cwd)
  
  return systemsList

def saveFiles(systemsList):
  """
  Saves systems as xyz files
  
  """
  
  cwd = os.getcwd()
  os.system("rm -rf %s" % (_topDir))
  IO.checkDirectory(_topDir, True)
  os.chdir(_topDir)
  
  systemListLen = len(systemsList)
  
  for i in range(systemListLen):
    
    fileName = "%03d_%s.xyz" % (i+1, systemsList[i].name)
    
    IO.writeXYZ(systemsList[i], fileName)
    
    hashkeyRadius = Atoms.getRadius(systemsList[i]) + 0.4
        
    cmdLine = "python ~/git/hkg/hkg.py %s %f" % (fileName, hashkeyRadius)
    
    hashkey = os.popen(cmdLine).read().strip()
    
    systemsList[i].hashkey = copy.deepcopy(hashkey)
  
  os.chdir(cwd)

def sortSystems(systemsList):
  """
  Sorts systems according to their energy
  
  """
  
  systemListLen = len(systemsList)
  
  for i in range(systemListLen):
    for j in range(systemListLen):
      if systemsList[i].totalEnergy < systemsList[j].totalEnergy:
        tempSystem = copy.deepcopy(systemsList[i])
        
        systemsList[i] = copy.deepcopy(systemsList[j])
        systemsList[j] = copy.deepcopy(tempSystem)
  
if __name__ == "__main__":
  
  # looks for FHI-aims simulations in the output directory
  fhiaimsDirs = getFileList()
  
  # reads in the FHI-aims systems
  systems = readFHIaimsSystems(fhiaimsDirs)
  
  # sort the systems according to energy
  sortSystems(systems)
  
  # save files
  saveFiles(systems)
  
  # generate statistics
  generateStatistics(systems)
  
  
  print "Finished."
  