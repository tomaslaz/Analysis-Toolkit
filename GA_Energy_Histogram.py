#!/usr/bin/env python

"""
A script to plot an energy distribution histogram from the last iteration of a KLMC GA simulation.

The script must be executed at the top level of the KLMC simulation directory

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import time

_energyRange = 2.0
_energyBinSize = 0.02

# Settings below should not be changed if you are not sure what they are for
_statsPrefix = "gaStatistics"
_statsSubfix = ".csv"

_runDirPath = "run/"
_zipSubFix = ".tar.gz"

_analyzeKLMCDir = "analysis"

def analyse(statsFile):
  """
  Reads in the statistics file, puts the data into bins and plots the histogram
  
  """
  
  # reads the energy values
  energyArr = readGAStatsFile(statsFile)
  
  # puts the data into bins
  binTheEnergies(energyArr)
  
def binTheEnergies(energyArr):
  """
  Generates bins and puts the energies into the bins
  
  """
  
  fig = plt.figure(figsize=(9, 6))
  ax = fig.add_subplot(1,1,1)
  
  minFreq = 0
  maxFreq = 0
  
  minEnergy = np.min(energyArr)
  maxEnergy = np.max(energyArr)
    
  if (maxEnergy - minEnergy > _energyRange):
    maxEnergy = minEnergy + _energyRange
  
  maxEnergy += _energyBinSize
  
  print "Bins are created within range: %f - %f with a stepsize %f" % (minEnergy, maxEnergy, _energyBinSize)
  
  bins = np.array(list(frange(minEnergy, maxEnergy, _energyBinSize)))
  
  inds = np.digitize(energyArr, bins)
  
  n, binsHist, patches = plt.hist(energyArr, bins=bins, normed=False, facecolor='r', alpha=0.75)
    
  maxFreq = np.max(n)

  plt.xlabel('Energy Bin (eV)', fontsize=18)
  plt.ylabel('Number of Configurations', fontsize=18)
  
  plt.grid(True)
  
  ax.set_xlim([minEnergy, maxEnergy])
  ax.xaxis.set_ticks(np.arange(minEnergy, maxEnergy, _energyBinSize*10))
  
  plt.subplots_adjust(left=0.08, bottom=0.11, top=0.98, right=0.98)
  
  fig.savefig('GA_Energy_Histogram.png')

def checkDirectory(dirPath, createMd=0):
  """
  Checks if the directory exists
  
  """
  
  exists = os.path.exists(dirPath)
  
  if ((not exists) and (createMd)):
    try:
      os.makedirs(dirPath)
      exists = True
      
    except:
      exists = False
      
  return exists

def getLastStatsFile():
  """
  Looks for zippped files and unzips them. Returns a list of directories.
  
  """
    
  cwd = os.getcwd()
  os.chdir(_runDirPath)
  
  maxStep = -1
  zipPath = None
  
  # finding the last iteration of the GA simulation
  for root, _, files in os.walk("./"):
    for fileName in files:
      if fileName.endswith(_zipSubFix):
        
        try:
          stepNo = int(fileName.split(".")[0])
        except:
          stepNo = -1
        
        if maxStep < stepNo:
          maxStep = stepNo
          zipPath = fileName
  
  if maxStep < 0:
    sys.exit("Error: could not find any data from a GA simulation")
  
  # extracting the zip file of the last iteration
  cmdLine = "gzip -d < %s | tar xf -" % (zipPath)
  os.system(cmdLine)
  
  # copy the analysis file from the zip file
  analysisFile = "%s%s%s" % (_statsPrefix, maxStep, _statsSubfix)
  analysisFilePath = os.path.join(cwd, _analyzeKLMCDir, analysisFile)
  
  cmdLine = "cp -rf %s %s" % (os.path.join(_runDirPath, str(maxStep), analysisFile), analysisFilePath)
  os.system(cmdLine)
  
  # clean up
  cmdLine = "rm -rf %s" % (_runDirPath)
  os.system(cmdLine)
  
  os.chdir(cwd)
  
  return analysisFilePath

def frange(a, b, step):
  """
  Range for float numbers
  
  """
  
  while a + sys.float_info.epsilon < b:
    yield a
    a += step

def main():
  timeStart = time.time()
  
  # does the analysis directory exist?
  checkDirectory(_analyzeKLMCDir, createMd=1)
  
  # finds the last iteration
  statsFile = getLastStatsFile()
  
  # reads in the statistics and plots the data
  analyse(statsFile)
  
  # Final message
  printFinalMessage(timeStart)
  
def printFinalMessage(timeStart=None):
  """
  Prints the author and the time taken to run the script
  
  """
  
  print
  
  if timeStart is not None:
    print "Finished analysing: in %f s." % (time.time() - timeStart)
  
  print
  print "-" * 30
   
  print "Author: Tomas Lazauskas, 2016"
  print "Website: www.lazauskas.net"
  print

def readGAStatsFile(filePath):
  """
  Reads in the GA iteration statistics file.
  
  """
  
  lineCnt = 0
  structCnt = 0
  
  lowEnergy = 0.0
  highEnergy = -1000000000000.0

  energiesTmp = np.zeros(10000, np.float64)
    
  statusUsed = 1
  
  try:
    f = open(filePath, "r")
    
  except:
    sys.exit("Error: cannot open: " + filePath)
  
  for line in f:

    if lineCnt > 2:
      line = line.strip()
      array = line.split(",")
      
      if statusUsed:
        edfn = int(array[3])
        status = int(array[4])
        energy = float(array[5])
        
      else:
        edfn = int(array[3])
        status = 1
        energy = float(array[4])
      
      if edfn > 0 and status == 1 and energy != 0.0:
        structCnt += 1
      
        energiesTmp[structCnt-1] = energy
        
        if energy < lowEnergy:
          lowEnergy = energy
          
        if energy > highEnergy:
          highEnergy = energy
        
    lineCnt += 1
    
  f.close()
  
  print "Lowest energy structure: %f" % (lowEnergy)
  print "Highest energy structure: %f" % (highEnergy)
  print "Total number of unique structures: %d" % (structCnt)
  
  energies = np.zeros(structCnt, np.float64)
  for i in range(structCnt): energies[i] = energiesTmp[i]

  return energies

if __name__ == "__main__":
  
  main()