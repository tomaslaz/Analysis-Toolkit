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
  
  def __init__(self, system, rdfCutOff, rdfStepSize, sigma, pairs):
    """
    Constructor
    
    """
    
    self.__rdfCutOff = rdfCutOff
    self.__rdfCutOffSq = self.__rdfCutOff**2
    self.__rdfStepSize = rdfStepSize
    self.__sigma = sigma
    self.__setPairs = pairs
    
    self.system = system
    
    self.maxRdfDist = getRdfBoxNumber(self.__rdfCutOff, self.__rdfStepSize) + 1
    self.NAtoms = self.system.NAtoms
    self.specieList = self.system.specieList
    self.specieListLen = len(self.specieList)
    
    self.NPairs = 0
    self.pairs, self.pairsIdxs = self.__createPairs(self.__setPairs)
    
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
    
    # finding the nearest neighbours
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
        
        if pairIdx is not None:
          self.pairs[pairIdx].ndist[boxNum] += 1
                
    if self.system.PBC[0] and self.system.PBC[0] and self.system.PBC[0]:
      volume = self.system.cellDims[0] * self.system.cellDims[1] * self.system.cellDims[2]
    else:  
      volume = 1.0
    
    rho = atomsCount / volume
    grConst = 4.0 * math.pi * rho * np.sqrt(2*np.pi)

    for i in range(self.NPairs):      
      for j in range(self.maxRdfDist):
        
        gr = 0.0
        for k in range(self.maxRdfDist):
          
          try:       
            gr += ((1.0 / (grConst * self.__sigma * np.power(self.rdist[k], 2.0))) * 
                   self.pairs[i].ndist[k] * 
                   (np.exp(-np.power(self.rdist[j] - self.rdist[k], 2.0) / (2.0 * np.power(self.__sigma, 2.)))))
          except:
            gr += 0.0
                
        self.pairs[i].gr[j] = gr

        # averaging the values
        if atomsCount > 0:
          self.pairs[i].gr[j] /= atomsCount
        
  def __createPairs(self, setPairs):
    """
    Creates an array of all possible pairs
    
    """
    
    setPairsArr = setPairs.split(",")
    setPairsArrLen = len(setPairsArr)
        
    pairs = []
    pairsIdx = {}
    
    pairsCnt = 0
    
    if setPairsArrLen > 0:
      for k in range(setPairsArrLen):
        for i in range(self.specieListLen-1):
          for j in xrange(i+1, self.specieListLen):
                        
            pairName =  "%s-%s" % (self.specieList[i], self.specieList[j])
            pairNameInv = "%s-%s" % (self.specieList[j], self.specieList[i])
            
            if pairName == setPairsArr[k] or pairNameInv == setPairsArr[k]:
              pairsIdx[pairName] = pairsCnt
              pairsCnt += 1
                
              pair = Pair(self.specieList[i], self.specieList[j], self.maxRdfDist)
              pairs.append(pair)
    
    else:
      for i in range(self.specieListLen-1):
        for j in xrange(i+1, self.specieListLen):
                      
          pairName =  "%s-%s" % (self.specieList[i], self.specieList[j])
          pairNameInv = "%s-%s" % (self.specieList[j], self.specieList[i])
          
          pairsIdx[pairName] = pairsCnt
          pairsCnt += 1
            
          pair = Pair(self.specieList[i], self.specieList[j], self.maxRdfDist)
          pairs.append(pair)
                      
    for i in range(self.specieListLen):
      
      pairName =  "%s-%s" % (self.specieList[i], self.specieList[i])
      
      foundPair = False
      
      if setPairsArrLen > 0:
        for k in range(setPairsArrLen):
          if pairName == setPairsArr[k]:
            foundPair = True
            break
      
      if not foundPair:
        continue
      
      pairsIdx[pairName] = pairsCnt
      pairsCnt += 1
          
      pair = Pair(self.specieList[i], self.specieList[i], self.maxRdfDist)
      pairs.append(pair)
    
    if pairsCnt < 1:
      sys.exit("Could not match the pairs in the system.")
    
    pairs.sort()
    
    # self.NPairs = self.specieListLen + int(math.floor(((self.specieListLen-1)*self.specieListLen)/2))
    self.NPairs = pairsCnt 
    
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
      try:
        pair = "%s-%s" % (atomBSym, atomASym)
        pairIdx = self.pairsIdxs[pair]
      except:
        pairIdx = None
      
    return pairIdx
  
  def __generateColours(self):
    """
    Generate a list of colours for plotting
    
    """
    
    self._colours = ['b', 'r', 'g', 'y', 'c', 'm', 'darkblue', 'sienna', 'indigo', 'orange', 'grey', 'brown']

  def _plotRDF(self):
    """
    Plotting RDF for all pairs
    
    """
    
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1,1,1)
  
    plt.subplots_adjust(left=0.02, bottom=0.13, top=0.95, right=0.95)
    plt.grid()
    
    labels = []
    series = []
    
    for i in range(self.NPairs):   
      label = self.pairs[i].pairName
            
      serie, = ax.plot(self.rdist, self.pairs[i].gr, c=self._colours[i], label=label, linewidth=2.5)
      
      labels.append(label)
      series.append(serie)
      
    plt.legend(series, labels, loc=1, fontsize=22)
  
    ax.set_xlabel(r'r($\AA$)', fontsize=22)
    ax.set_xlim([0, self.__rdfCutOff])
    ax.xaxis.set_ticks(np.arange(0, self.__rdfCutOff+1, 1.0))
    
    #ax.set_yscale('log')
    #ax.set_ylabel(r'g(r)', fontsize=18)
    #ax.set_ylim([10**-7, 10**-3])
    ax.yaxis.set_ticks([])
    
    fig.savefig('%s.png' % (self.system.name), dpi=300, bbox_inches='tight')
        
  def __rdfBoxes(self):
    """
    Creates an array with RDF box lenghts
    
    """
    
    rdfBoxArr = np.zeros(self.maxRdfDist, np.float64)
    
    for i in range(self.maxRdfDist):
      rdfBoxArr[i] = (i+1) * self.__rdfStepSize
            
    return rdfBoxArr

def cmdLineArgs():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog [options] inputFile"
  
  parser = OptionParser(usage=usage)

  parser.add_option("-r", "--radius", dest="rdfCutOff", default=10.0, type="float",
    help="RDF cut-off radius. Default = 10.0 A")
  
  parser.add_option("-x", "--step", dest="rdfCStepsize", default=0.01, type="float",
    help="RDF stepsize. Default = 0.01 A")
  
  parser.add_option("-s", "--sigma", dest="gausSigma", default=0.1, type="float",
    help="Sigma value for Gaussian smearing. Default = 0.1")
  
  parser.add_option("-p", "--pairs", dest="pairs", default="", type="string",
    help="A list of pairs for which RDF is plotted. Default = ''")

  parser.add_option("-c", "--cubic", dest="cubicPBC", default=False, action="store_true",
    help="Apply cubic periodicity. Default = False")

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
  
  options, args = cmdLineArgs()
  
  fileName = args[0]
  
  print "Reading: ", fileName
  
  system = IO.readSystemFromFileXYZ(fileName)
  
  # the system should not be moved if periodic boundaries are applied
  if not system.PBC[0] and not system.PBC[1] and not system.PBC[2]:
    system.calcCOG()
    system.moveToCOG()
    
  systemRDF = RDF(system, options.rdfCutOff, options.rdfCStepsize, options.gausSigma, options.pairs)
  
  systemRDF._calcRDF()
  
  systemRDF._plotRDF()
  
  print "Finished."