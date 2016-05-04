#!/usr/bin/env python

"""
A script to calculate the radial distribution function

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

import copy
import math
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser

import source.Atoms as Atoms
import source.IO as IO

_constRdfCutOff = 15.0
_constRdfStepSize = 0.2

class Pair(object):
  """
  A class to save the rdf data of a pair
  
  """
  def __init__(self, el1, el2, maxRdfDist):
    """
    Constructor
    
    """
    
    self.el1 = el1
    self.el2 = el2
    self.pairName = "%s-%s" % (self.el1, self.el2)
    self.ndist = np.zeros(maxRdfDist, np.int32)
    self.gr = np.zeros(maxRdfDist, np.float64)

class RDF(object):
  """
  A RDF class
  
  """
  
  def __init__(self, system):
    """
    Constructor
    
    """
    
    self.__rdfCutOff = _constRdfCutOff
    self.__rdfCutOffSq = self.__rdfCutOff**2
    self.__rdfStepSize = _constRdfStepSize
    
    
    self.system = system
    
    self.maxRdfDist = getRdfBoxNumber(_constRdfCutOff, _constRdfStepSize) + 1
    self.NAtoms = self.system.NAtoms
    self.specieList = self.system.specieList
    self.specieListLen = len(self.specieList)
    
    self.NPairs = 0
    self.pairs, self.pairsIdxs = self.__createPairs()
    
    self.rdist = self.__rdfBoxes()
    
    self.__contstrains = False
    
    self._atomsToAnalyse = None
    self._atomsToAnalyseCnt = 0
    
    self.__findAtomsToAnalyse()
    
    self.__generateColours()
    
  def _calcRDF(self):
    """
    Estimates RDF for all possible pairs
    
    """
    
    atomsCount = self._atomsToAnalyseCnt
    
    for i in range(atomsCount):

      atomIdx = self._atomsToAnalyse[i]
            
      atomASym = self.specieList[self.system.specie[atomIdx]]
            
      # Finding the neighbouring atoms 
      neighboursCnt, neighboursArr, neighboursDistArr = self.system.findNN(atomIdx, self.__rdfCutOffSq)
            
      for j in range(neighboursCnt):
        atomBIdx = neighboursArr[j]
        atomBSym = self.specieList[self.system.specie[atomBIdx]]
        dist = neighboursDistArr[j]
        
        boxNum = getRdfBoxNumber(dist, self.__rdfStepSize)
        
        pairIdx = self.__getPairIdx(atomASym, atomBSym)
        self.pairs[pairIdx].ndist[boxNum] += 1
                  
    volume = 1.0
    grConst = 4.0 * math.pi * (atomsCount / volume)
    
    for i in range(self.NPairs):
      for j in range(self.maxRdfDist):
        distSq = self.rdist[j]**2
        
        self.pairs[i].gr[j] = self.pairs[i].ndist[j] / (grConst * distSq * self.__rdfStepSize)
        
        # averaging the values
        self.pairs[i].gr[j] /= atomsCount
        
  def __createPairs(self):
    """
    Creates an array of all possible pairs
    
    """
    
    pairs = []
    pairsIdx = {}
    
    self.NPairs = self.specieListLen + int(math.floor(((self.specieListLen-1)*self.specieListLen)/2))
    pairsCnt = 0
    
    for i in range(self.specieListLen-1):
      for j in xrange(i+1, self.specieListLen):
        
        pair = Pair(self.specieList[i], self.specieList[j], self.maxRdfDist)
        pairName =  "%s-%s" % (self.specieList[i], self.specieList[j])
        pairsIdx[pairName] = pairsCnt
        pairsCnt += 1
        
        pairs.append(pair)
        
    for i in range(self.specieListLen):
      
      pair = Pair(self.specieList[i], self.specieList[i], self.maxRdfDist)
      pairName =  "%s-%s" % (self.specieList[i], self.specieList[i])
      pairsIdx[pairName] = pairsCnt
      pairsCnt += 1
      
      pairs.append(pair)
            
    return pairs, pairsIdx
  
  def __findAtomsToAnalyse(self):
    """
    Generates a list of atoms which are within constrains
    
    """
    
    cntrRad = 15.0
    cntrRadSq = cntrRad**2
    
    atomCnt = 0
    atomsList = []
    
    if self.__contstrains:      
      for i in range(self.NAtoms):
        posx = self.system.pos[3*i+0]
        posy = self.system.pos[3*i+1]
        posz = self.system.pos[3*i+2]
        
        distSq = posx**2 + posy**2 + posz**2
          
        if distSq <= cntrRadSq:
          atomCnt += 1
          atomsList.append(i)

      self._atomsToAnalyseCnt = atomCnt
      self._atomsToAnalyse = np.array(atomsList, np.int32)
      
    else:
      self._atomsToAnalyse = np.array(range(self.NAtoms), np.int32)
      self._atomsToAnalyseCnt = self.NAtoms
  
  def __getPairIdx(self, atomASym, atomBSym):
    """
    Returns a pair-idx for an atom pair
    
    """

    try:
      pair = "%s-%s" % (atomASym, atomBSym)
      pairIdx = self.pairsIdxs[pair]
      
    except:
      pair = "%s-%s" % (atomBSym, atomASym)
      pairIdx = self.pairsIdxs[pair]
        
    return pairIdx
  
  def __generateColours(self):
    """
    Generate a list of colours for plotting
    
    """
    
    self._colours = ['b', 'g', 'y', 'c', 'm', 'r', 'darkblue', 'sienna', 'indigo', 'orange', 'grey', 'brown']

  def _plotRDF(self):
    """
    Plotting RDF for all pairs
    
    """
    
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1,1,1)
  
    plt.subplots_adjust(left=0.1, bottom=0.11, top=0.95, right=0.95)
    plt.grid()
    
    labels = []
    series = []
    
    for i in range(self.NPairs):   
      label = self.pairs[i].pairName
      
      if label not in ['Si-C', 'Ga-N']:
        continue
      
      serie, = ax.plot(self.rdist, self.pairs[i].gr, c=self._colours[i], label=label, linewidth=2.0)
      
      labels.append(label)
      series.append(serie)
      
    plt.legend(series, labels, loc=1)
    
    ax.set_yscale('log')
    
    ax.set_xlabel(r'r($\AA$)', fontsize=18)
    ax.set_ylabel(r'g(r)', fontsize=18)
    
    ax.set_ylim([10**-7, 10**-3])
    
    ax.set_xlim([0, self.__rdfCutOff])
    ax.xaxis.set_ticks(np.arange(0, self.__rdfCutOff+1, 1.0))
    
    fig.savefig('%s.png' % (self.system.name))
        
  def __rdfBoxes(self):
    """
    Creates an array with RDF box lenghts
    
    """
    
    rdfBoxArr = np.zeros(self.maxRdfDist, np.float64)
    
    for i in range(self.maxRdfDist):
      rdfBoxArr[i] = (i+1) * _constRdfStepSize
            
    return rdfBoxArr

def cmdLineArgs():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog inputFile"
  
  parser = OptionParser(usage=usage)

  parser.disable_interspersed_args()
      
  (options, args) = parser.parse_args()

  if (len(args) != 1):
    parser.error("incorrect number of arguments")

  return options, args
    
def getRdfBoxNumber(distance, boxSize):
  """
  A function to get the box number for the RDF function
  
  """
  
  return int(math.floor(distance/boxSize))

if __name__ == "__main__":
  
  _, args = cmdLineArgs()
  
  fileName = args[0]
  
  print "Reading: ", fileName
  
  system = IO.readSystemFromFileXYZ(fileName)
  
  system.calcCOG()
  system.moveToCOG()
    
  systemRDF = RDF(system)
  
  systemRDF._calcRDF()
  
  systemRDF._plotRDF()
  
  print "Finished."