"""
A script to plot the energy evolution of the lowest energy structures of a KLMC GA simulation.

The script must be executed at the top level of the KLMC simulation directory

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import sys, os, copy
import numpy as np
import matplotlib.pyplot as plt
import re
import time

# Number of structres to be presented on the graph
lowEnergyCnt = 20
# Energy range to be presented on the graph
energyRange = 10.0

# Settings below should not be changed if you are not sure what they are for
#---------------------------------------------------------------------------
_analyzeKLMCDir = "analysis"

_statsPrefix = "gaStatistics"
_statsSubfix = ".csv"

class GAIterStats:
  """
  Class to hold statistics from a GA iteration
  """

  def __init__(self, iter, structCnt, lowEnergy, avgEnergy, highEnergy, sigmaEnergy, lowEnergyArr):
    self.iter = iter
    self.structCnt = structCnt
    self.lowEnergy = lowEnergy
    self.avgEnergy = avgEnergy
    self.highEnergy = highEnergy
    self.sigmaEnergy = sigmaEnergy
    self.lowEnergyArr = lowEnergyArr

def analyseGARun(gaItNoList, filesList, unzipFlag=False):
  """
  Reads in the statistics files and generates statistics of the GA run.

  """

  totGaIter = len(gaItNoList)

  statistics = []
  prevIterStats = None

  for i in range(totGaIter):

    if i > 0:
      prevIterStats = gaIterStat

    filePath = getIterFile(i, totGaIter, gaItNoList, filesList)

    gaIterStat = readGAStatsFile(i, filePath, prevIterStats)

    statistics.append(gaIterStat)

  plotGASimStatistics(statistics)

  return statistics

def getIterFile(iter, totGaIter, gaItNoList, filesList):
  """
  Returns a filepath from filesList with respect to the generation number.

  """

  filePath = None

  for i in range(totGaIter):

    if gaItNoList[i] == iter:

      filePath = filesList[i]

      return filePath

def getFileList(dirPath=None):
  """
  Checks the directory for ga statistics files and makes a list of them

  """

  filesList = []
  gaItNoList = []

  if (dirPath is not None):
    cwd = os.getcwd()
    os.chdir(dirPath)

  else:
    cwd = None

  subfixLen = len(_statsSubfix)
  prefixLen = len(_statsPrefix)

  for root, _, files in os.walk("./"):

    for fileName in files:

      if fileName.startswith("ga") and fileName.endswith(".csv"):
        fileNameLen = len(fileName)

        gaItNo = int(fileName[prefixLen:fileNameLen-subfixLen])

        gaItNoList.append(gaItNo)

        filesList.append(os.path.join(_analyzeKLMCDir, fileName))

  if (dirPath is not None):
    os.chdir(cwd)

  return gaItNoList, filesList

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

def plotGASimStatistics(statistics):
  """
  Plots statistics regarding a GA simulation.

  """

  lenStats = len(statistics)

  maxy = -1000000000000000.0
  miny = maxy * -1

  x = np.zeros(lenStats, np.int16)
  s = np.zeros(lenStats, np.int16)
  a = np.zeros(lenStats, np.float64)
  l = np.zeros(lenStats, np.float64)
  e = np.zeros(lenStats, np.float64)

  la = np.zeros(lenStats, np.float64)
  lm = np.zeros([lenStats, lowEnergyCnt], np.float64)

  for i in range(lenStats):
    x[i] = statistics[i].iter
    s[i] = statistics[i].structCnt
    a[i] = statistics[i].avgEnergy
    l[i] = statistics[i].lowEnergy
    e[i] = statistics[i].sigmaEnergy

    la[i] = 0

    for j in range(lowEnergyCnt):
      lm[i, j] = statistics[i].lowEnergyArr[j]

      if lm[i, j] > maxy and lm[i, j] < 0.0:
        maxy = lm[i, j]

      if lm[i, j] < miny and lm[i, j] < 0.0:
        miny = lm[i, j]

      la[i] += statistics[i].lowEnergyArr[j]

    la[i] /= lowEnergyCnt

  fig = plt.figure(figsize=(9, 6))
  ax1 = fig.add_subplot(1,1,1)

  for j in range(lowEnergyCnt):
    ax1.scatter(x, lm[:, j], s=10, c='black', marker="s")

  ax1.plot(x, la, c='r', linewidth=3.0)

  ax1.set_xlabel('GA iteration', fontsize=18)
  ax1.set_ylabel('Energy (eV)', fontsize=18)

  ystart, yend = ax1.get_ylim()

  maxy = miny + energyRange
  print "Plotting energy range: ", miny, maxy

  ystart = miny
  yend   = maxy

  #lenStats = 2500

  ax1.set_xlim([0, lenStats])

  ax1.xaxis.set_ticks(np.arange(0, lenStats+1, int(round(lenStats/10))))

  ax1.set_ylim([ystart, yend])

  ax1.yaxis.set_ticks(np.arange(ystart, yend, energyRange/25.0))

  plt.subplots_adjust(left=0.13, bottom=0.11, top=0.98, right=0.97)

  plt.grid()

  fig.savefig('GA_top_analysis.png')

def readGAStatsFile(gaIter, filePath, prevStats):
  """
  Reads in the GA iteration statistics file.

  """

  lineCnt = 0
  structCnt = 0
  lowEnergy = 0.0
  avgEnergy = 0.0
  highEnergy = 0.0
  sigmaEnergy = 0.0

  energySum = 0.0
  energiesTmp = np.zeros(1000, np.float64)

  lowEnergyArr = np.zeros(lowEnergyCnt, np.float64)

  statusUsed = 1

  # if the stats file is missing we have to use stats from the previous step
  if filePath == None:
    if prevStats is not None:
      gaIterStat = GAIterStats(prevStats.iter, prevStats.structCnt, prevStats.lowEnergy, prevStats.avgEnergy,
                               prevStats.highEnergy, prevStats.sigmaEnergy, prevStats.lowEnergyArr)

    else:
      gaIterStat = GAIterStats(gaIter, structCnt, lowEnergy, avgEnergy, highEnergy, sigmaEnergy, lowEnergyArr)

    return gaIterStat

  try:
    f = open(filePath, "r")
  except:
    sys.exit(__name__ + ": error: cannot open <" + filePath)

  for line in f:

    if lineCnt > 0:
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

      if edfn > 0 and status == 1 and energy != 0.0 and structCnt < lowEnergyCnt:
        structCnt += 1

        energySum += energy

        energiesTmp[structCnt-1] = energy

        if energy < lowEnergy:
          lowEnergy = energy

        if energy > highEnergy:
          highEnergy = energy

    lineCnt += 1

  f.close()

  energies = np.zeros(structCnt, np.float64)

  for i in range(structCnt): energies[i] = energiesTmp[i]

  sigmaEnergy = np.std(energies)

  if structCnt > 0: avgEnergy = energySum / float(structCnt)

  # Preparing top x structures' energies
  # sorting according to the energy
  energiesRanked = np.sort(copy.deepcopy(energies))

 # print gaIter, energiesTmp

  noTopStructures = min(len(energiesRanked), lowEnergyCnt)

  lowEnergyArr[:noTopStructures] = copy.deepcopy(energiesRanked[:noTopStructures])

  gaIterStat = GAIterStats(gaIter, structCnt, lowEnergy, avgEnergy, highEnergy, sigmaEnergy, lowEnergyArr)

  return gaIterStat

def getZipFileList():
  """
  Looks for zippped files and unzips them. Returns a list of directories.

  """

  dirName = "run/"

  subFix = ".tar.gz"

  cwd = os.getcwd()
  os.chdir(dirName)

  for root, _, files in os.walk("./"):
    for fileName in files:
      if fileName.endswith(subFix):

        stepNo = fileName.split(".")[0]
        analysisFile = "%s%s%s" % (_statsPrefix, stepNo, _statsSubfix)
        analysisFilePath = os.path.join(cwd, _analyzeKLMCDir, analysisFile)

        if (not os.path.isfile(analysisFilePath)):

          cmdLine = "gzip -d < %s | tar xf -" % (fileName)

          os.system(cmdLine)

          cmdLine = "cp -rf %s %s" % (os.path.join(dirName, stepNo, analysisFile), analysisFilePath)

          os.system(cmdLine)

  cmdLine = "rm -rf %s" % (dirName)
  os.system(cmdLine)

  os.chdir(cwd)

if __name__ == "__main__":

  checkDirectory(_analyzeKLMCDir, createMd=1)

  unzipFlag = 1

  if unzipFlag:
    getZipFileList()

  gaItNoList, filesList = getFileList(dirPath=_analyzeKLMCDir)

  analyseGARun(gaItNoList, filesList)
