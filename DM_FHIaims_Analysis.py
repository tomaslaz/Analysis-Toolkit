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
from optparse import OptionParser

import source.Atoms as Atoms
import source.Fhiaims as FHIaims
import source.File as File
import source.IO as IO

_fhiaimsGeometryFile = "geometry.in"
_fhiaimsOutFile = "FHIaims.out"
_outputDir = "output"
_topDir = "tops"
_uniqueDir = "unique"

def cmd_line_args():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog "
  
  parser = OptionParser(usage=usage)
    
  parser.add_option("-r", "--relaxed", dest="relaxed", action="store_true", default=False, 
    help="Whether geometry relaxation was used.")
    
  parser.disable_interspersed_args()
      
  (options, args) = parser.parse_args()
  
  return options, args

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

def generateStatistics(systemlist, unique=False):
  """
  Generates statistics about the FHI-aims simulations
    
  """
    
  noOfSystems = len(systemlist)
  sumOfCores = 0
  runTimes = np.zeros(noOfSystems, np.float64)
  
  if not unique:
    f = open("Stats.csv", "w")
  else:
    f = open("%s/Stats.csv" % (_uniqueDir), "w")
  
  f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % ("System", "Energy", "Hashkey", "Cores", 
          "Time", "Tot.Time", "H-L", 
          "VBM", "VBMOcc", "VBMSpinChannel", 
          "CBM", "CBMOcc", "CBMSpinChannel", 
          "SpinN", "SpinS", "SpinJ","Size"))
  
  systemCnt = 0
  
  for system in systemlist:
    f.write("%s,%f,%s,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n" % (system.name, system.totalEnergy, system.hashkey,
                                     system.noOfcores, system.runTime, system.noOfcores*system.runTime,
                                     system.homo_lumo_gap, system.vbm, system.vbm_occ_num, system.vbm_spin_chan, 
                                     system.cbm, system.cbm_occ_num, system.cbm_spin_chan, 
                                     system.spin_N, system.spin_S, system.spin_J, float(system.NAtoms)))
    
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
  
def readFHIaimsSystems(fhiaimsDirs, relaxed=False):
  """
  Reads in the FHI-aims systems.
  
  """
  
  systemsList = []

  for dirName in fhiaimsDirs:
    cwd = os.getcwd()
    
    os.chdir(dirName)
    #systemName = dirName[len(_outputDir)+1:-2].strip()
    systemName = dirName.strip()
   
    success, error, system = FHIaims._readAimsStructure(_fhiaimsGeometryFile, _fhiaimsOutFile, relaxed=relaxed)
    
    if success:
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
  os.system("rm -rf %s" % (_uniqueDir))
  
  IO.checkDirectory(_uniqueDir, True)
  IO.checkDirectory(_topDir, True)
  os.chdir(_topDir)
  
  systemListLen = len(systemsList)
  
  uniquehashkeys = []
  uniqueSystems = []
  uniqueCnt = 0
   
  for i in range(systemListLen):
    
    nStr = "n%02d" % (systemsList[i].NAtoms)
    
    #fileName = "%s_%03d_%s.xyz" % (nStr, i+1, systemsList[i].name)
    fileName = "%s.xyz" % (systemsList[i].name)
        
    success_, error_ = File.writeXYZ(systemsList[i], fileName)
    
    hashkeyRadius = Atoms.getRadius(systemsList[i]) + 1.0

    cmdLine = "python ~/git/hkg/hkg.py %s %f" % (fileName, hashkeyRadius)
    
    hashkey = os.popen(cmdLine).read().strip()
    
    systemsList[i].hashkey = copy.deepcopy(hashkey)
    
    # saving only uniques:
    if not (hashkey in uniquehashkeys):
      uniqueCnt += 1
      uniquehashkeys.append(hashkey)
      
      uniqueSystems.append(systemsList[i])
      
      fileName2 = "%s_%03d_%s.xyz" % (nStr, uniqueCnt, systemsList[i].name)
      
      os.system("cp %s ../%s/%s" % (fileName, _uniqueDir, fileName2))
  
  os.chdir(cwd)
  
  return uniqueSystems

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
  
  # reading the command line arguments and options
  options, _ = cmd_line_args()
  
  # looks for FHI-aims simulations in the output directory
  fhiaimsDirs = getFileList()
    
  # reads in the FHI-aims systems
  systems = readFHIaimsSystems(fhiaimsDirs, options.relaxed)
  
  # sort the systems according to energy
  sortSystems(systems)
  
  # save files
  uniqueSystems = saveFiles(systems)
  sortSystems(uniqueSystems)
  
  # generate statistics
  generateStatistics(systems)
  
  # generating unique statistics
  generateStatistics(uniqueSystems, True)
  
  print "Finished."
  