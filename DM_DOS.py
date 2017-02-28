#!/usr/bin/env python

"""
A script to plot DOS (integrated)

@author Tomas Lazauskas, David Mora Fonz, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import math
import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
import scipy.special
import sys

_delta = 0.01
_extraBins = 2

_sigma = 0.1
_kB = 0.0257/298.0

_colours = ['r', 'b', 'g', 'y', 'c', 'm', 'darkblue', 'sienna', 'indigo', 'orange', 'grey', 'brown']

def cmdLineArgs():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog inputFile"
  
  parser = OptionParser(usage=usage)

  parser.disable_interspersed_args()
  
  parser.add_option('-t', dest="temps", default=None, help="List of temperatures, separated by a comma (default t=0)")
    
  (options, args) = parser.parse_args()

  if (len(args) != 1):
    parser.error("incorrect number of arguments")

  return options, args

def getTheListOfTemps(tempString):
  """
  A function to process an argument string line to an array of temperatures
  
  """
  
  temps = None
  
  arrString = tempString.split(",")
  lenArrString = len(arrString)
  
  temps = np.zeros(lenArrString)
  
  for i in range(lenArrString):
    try:
      temp = float(arrString[i])
    except:
      sys.exit("Incorrect temperatures.")
    
    temps[i] = temp
  
  return temps

def plotDOS(energyBins, energyDOS, eMax):
  """
  Plots the DOS graph.
  
  """
  
  fig = plt.figure(figsize=(9, 6))
  ax1 = fig.add_subplot(1,1,1)
  plt.subplots_adjust(left=0.1, bottom=0.11, top=0.95, right=0.95)
  
  ax1.plot(energyBins, energyDOS, c='r', linewidth=2.0)
  
  plt.grid()
    
  stepSize = roundTo1St(eMax/10)
  
  ax1.xaxis.set_ticks(np.arange(0, eMax, stepSize))
  
  ax1.set_xlabel('Energy (eV)', fontsize=18)
  ax1.set_ylabel('DOS', fontsize=18)
  
  fig.savefig('DOS.png', dpi=300, bbox_inches='tight')

def plotDOSandIntegratedDOS(energyBins, energyDOS, tempArrs, noOfTemps, temps, eMax):
  """
  Plots DOS and integrated DOS on the same graph.
  
  """
  
  same = []
  labels = []
  
  fig = plt.figure(figsize=(9, 6))
  ax1 = fig.add_subplot(1,1,1)
  plt.subplots_adjust(left=0.1, bottom=0.11, top=0.95, right=0.90)
 
  ax1.set_ylabel('DOS', fontsize=18)
  
  # plotting DOS
  label = "DOS"
  series, = ax1.plot(energyBins, energyDOS, c=_colours[0], label=label, linewidth=3.0)
  
  labels.append(label)
  same.append(series)
 
  # plotting integrated DOS
  ax2 = ax1.twinx()
 
  for i in range(noOfTemps):
    label = "%d K" % (temps[i])
    series, = ax2.plot(energyBins, tempArrs[:, i], c=_colours[i+1], label=label, linewidth=2.0)
   
    labels.append(label)
    same.append(series)
   
  plt.grid()
   
  plt.legend(same, labels, loc=0)
   
  stepSize = roundTo1St(eMax/10)
   
  ax1.xaxis.set_ticks(np.arange(0, eMax, stepSize))
   
  ax1.set_xlabel('Energy (eV)', fontsize=18)
  ax2.set_ylabel('Integrated DOS', fontsize=18)
   
  fig.savefig('DOSandIntegratedDOS.png', dpi=300, bbox_inches='tight')

def plotIntegratedDOS(energyBins, tempArrs, noOfTemps, temps, eMax):
  """
  Plots the integrated DOS graph.
  
  """
  
  series = []
  labels = []
  
  fig = plt.figure(figsize=(9, 6))
  ax1 = fig.add_subplot(1,1,1)
  plt.subplots_adjust(left=0.1, bottom=0.11, top=0.95, right=0.95)

  for i in range(noOfTemps):
    label = "%d K" % (temps[i])
    serie, = ax1.plot(energyBins, tempArrs[:, i], c=_colours[i], label=label, linewidth=2.0)
  
    labels.append(label)
    series.append(serie)
  
  plt.grid()
  
  plt.legend(series, labels, loc=0, fontsize=18)
  
  stepSize = roundTo1St(eMax/10)
  
  ax1.xaxis.set_ticks(np.arange(0, eMax, stepSize))
  
  ax1.set_xlabel('Energy (eV)', fontsize=18)
  ax1.set_ylabel('Integrated DOS', fontsize=18)
  
  fig.savefig('Integrated_DOS.png', dpi=300, bbox_inches='tight')

def roundTo1St(x):
  """
  A function which rounds to the first significant number
  
  """
  
  return round(x, -int(math.floor(math.log10(abs(x)))))

def runDOS(energies, temps):
  """
  Calculates and plots DOS
  
  """
  
  # getting the number of temperatures
  noOfTemps = len(temps)
  
  # pushing by eMin
  energies -= energies.min()
  
  # getting the unique list of energies
  energiesUnique = np.unique(energies)
  
  # get min and max 
  eMin = energies.min()
  eMax = energies.max() + _extraBins*_delta
  
  # number of bins
  M = int((eMax - eMin)/_delta)
  
  # creating energy bins
  energyBins = np.arange(eMin, np.around([M * _delta], decimals=4), _delta)
  
  # preparing a DOS array
  energyDOS = np.zeros(M, dtype=np.float32)
  
  # creating temperature arrays
  tempArrs = np.zeros([M, noOfTemps], dtype=np.float32)
  
  # calculating DOS and integrated DOS with respect to the temperatures
  for i in range(M):
    energyDOS[i] = np.sum((1/(_sigma*np.pi**0.5)) * np.exp(-(energyBins[i] - energies)**2 / _sigma**2))
    
    if noOfTemps > 0:
      
      # going through the list of temperatures
      for j in range(noOfTemps):
        temp = temps[j]
         
        # calculating the integrated DOS       
        tempArrs[i][j] = ((1.0/(_sigma*np.pi**0.5)) * 
                          np.sum(((np.pi**0.5/(2*_sigma**-1)) * 
                                  (scipy.special.erf(1/_sigma * (energyBins[i] - energiesUnique)) - 
                                   scipy.special.erf(-np.infty))) * 
                                 np.exp(-(energiesUnique)/(_kB*temp) )))
  
  # printing DOS graph
  plotDOS(energyBins, energyDOS, eMax)    

  if noOfTemps > 0:
    # printing integrated DOS graph
    plotIntegratedDOS(energyBins, tempArrs, noOfTemps, temps, eMax)
    
    # printing DOS and integrated DOS
    plotDOSandIntegratedDOS(energyBins, energyDOS, tempArrs, noOfTemps, temps, eMax)

if __name__ == "__main__":
  
  options, args = cmdLineArgs()
  
  # getting the temperatures
  if options.temps is not None:
    temps = getTheListOfTemps(options.temps)
  
  else:
    temps = []
    
  # reading the energies from a file
  energies = np.loadtxt(args[0])
  
  runDOS(energies, temps)
  
  print "Finished."
  