#!/usr/bin/env python

"""
A script to analyse systems (xyz format) in terms of avg. bond distance and coordination

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

import copy
import math
import os
import sys

import numpy as np
from optparse import OptionParser

import source.Atoms as Atoms
import source.IO as IO
import source.System as System
import source.Utilities as Utilities
            
class Connectivity(object):
  """
  A class to determine the connectivity of the system
  
  """
  
  def __init__(self, system, superCell, minRadius, maxRadius):
    """
    Initializes the class.
    
    """
    
    self._error = 0
    self.vertices = []
    
    self.nVert = system.NAtoms
    
    self.nTotVertColors = len(system.specieList)
    self.vertColors = copy.copy(system.specieList)
    
    self.__createVertices(system)
    
    self.__getConnectivity(system, superCell, minRadius, maxRadius)
    
  def __createVertices(self, system):
    """
    Creates vertices according to every atom in the system
    
    """
    
    for i in range(system.NAtoms):
      
      vertex = Vertex(i, system.specie[i], self.nTotVertColors)
      
      self.vertices.append(vertex)
  
  def __getConnectivity(self, system, superCell, minRadius, maxRadius):
    """
    Determines the connectivity of a system
    
    """
    
    minRadiusSq = minRadius * minRadius
    maxRadiusSq = maxRadius * maxRadius
     
    for i in range(system.NAtoms-1):
      for j in range(i+1, system.NAtoms):

        distSq = Utilities.distanceSq(system.pos[3*i+0], system.pos[3*i+1], system.pos[3*i+2],
                                      system.pos[3*j+0], system.pos[3*j+1], system.pos[3*j+2])
        
        if distSq > minRadiusSq and distSq < maxRadiusSq:
          
          iAtomColor = system.specie[i]
          jAtomColor = system.specie[j]
          
          self.__updateVertConnect(i, jAtomColor, distSq)
          self.__updateVertConnect(j, iAtomColor, distSq)
  
    if superCell is not None:
      
      self._errMsg = "PBC are not implemented yet"
      self._error = True
  
  def __updateVertConnect(self, vertIdx, connColor, distSq):
    """
    If two vertices are connected - update their information
    
    """
    
    vertex = self.vertices[vertIdx]
    
    if vertex.vertAtomIdx != vertIdx:
      
      self._errMsg = ": indexes must match!"
      self._error = True
    
    vertex.vertConnectColors[connColor] += 1
    vertex.vertConnectColorsDist[connColor] += math.sqrt(distSq)
  
  def avgCoordNumber(self):
    """
    Estimates the average coordination number of a vertex
    
    """
    
    sumCoord = 0.0
    
    for i in range(self.nVert):
      vertex = self.vertices[i]
            
      vertCoord = vertex.getCoordNumber()
  
      sumCoord += vertCoord
    
    if (self.nVert > 0):
      avgCoord = sumCoord / self.nVert
    else:
      avgCoord = 0.0
    
    return avgCoord
  
  def avgCoordNumberBySpecie(self):
    """
    Estimates the average coordination according to the colour (atom type)
    
    """
  
    specieCoord = {}
    
    for i in range(self.nTotVertColors):
      
      vertCnt = 0.0
      vertSumCoor = 0.0
      
      vertColor = self.vertColors[i]
      
      for j in range(self.nVert):
        vertex = self.vertices[j]
        
        if (self.vertColors[vertex.vertColorIdx] == vertColor):
          vertCnt += 1.0
          vertSumCoor += vertex.getCoordNumber()
        
        if vertCnt > 0.0:
          avgCoord = vertSumCoor / vertCnt
        else:
          avgCoord = 0.0
        
        specieCoord[vertColor] = avgCoord
      
    return specieCoord
  
  def avgBondDistance(self):
    """
    Estimates average bond distances between all possible pairs
    
    """
    
    bondDist = {}
    bondCnt = {}
    
    # creating dictionaries
    for i in range(self.nTotVertColors):
      for j in range(i, self.nTotVertColors):
        bondDist["%s-%s" % (self.vertColors[i], self.vertColors[j])] = 0.0 
        bondCnt["%s-%s" % (self.vertColors[i], self.vertColors[j])] = 0.0
    
    # summing the bond distances and counting them 
    for i in range(self.nVert):
      vertex = self.vertices[i]
      
      for j in range(self.nTotVertColors):
        try:
          bondDist["%s-%s" % (self.vertColors[j], self.vertColors[vertex.vertColorIdx])] += vertex.vertConnectColorsDist[j]
          bondCnt["%s-%s" % (self.vertColors[j], self.vertColors[vertex.vertColorIdx])] += vertex.vertConnectColors[j]
          
        except:
          bondDist["%s-%s" % (self.vertColors[vertex.vertColorIdx], self.vertColors[j])] += vertex.vertConnectColorsDist[j]
          bondCnt["%s-%s" % (self.vertColors[vertex.vertColorIdx], self.vertColors[j])] += vertex.vertConnectColors[j]
    
    # calculating averages
    for i in range(self.nTotVertColors):
      for j in range(i, self.nTotVertColors):

        if bondDist["%s-%s" % (self.vertColors[i], self.vertColors[j])] != 0.0:
          bondDist["%s-%s" % (self.vertColors[i], self.vertColors[j])] /= bondCnt["%s-%s" % (self.vertColors[i], self.vertColors[j])]
        
    return bondDist

class Vertex(object):
  """
  A class to save the data about a vertex
  
  """
  
  def __init__(self, vertAtomIdx, vertColorIdx, totAtomSpecieCount):
    """
    Initiates the vertex class
    
    """
    
    self._error = False
    self._errMsg = ""
    
    self.vertAtomIdx = vertAtomIdx
    self.vertColorIdx = vertColorIdx
    self.vertConnectColors = np.zeros(totAtomSpecieCount, np.int32)
    self.vertConnectColorsDist = np.zeros(totAtomSpecieCount, np.float64)
  
  def getCoordNumber(self):
    """
    Gets the coordination number of an atom
    
    """
    
    return np.sum(self.vertConnectColors)

def cmdLineArgs():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog inputFile minRadius maxRadius"
  
  parser = OptionParser(usage=usage)

  parser.disable_interspersed_args()
      
  (options, args) = parser.parse_args()

  if (len(args) != 3):
    parser.error("incorrect number of arguments")

  return options, args

if __name__ == "__main__":
  
  _, args = cmdLineArgs()
  
  fileName = args[0]
  minRadius = float(args[1])
  maxRadius = float(args[2])
  
  print "Reading: ", fileName
  print "Analysis radius: ", minRadius, "-", maxRadius
  
  system = IO.readSystemFromFileXYZ(fileName)
  
  connO = Connectivity(system, None, minRadius, maxRadius)
  
  if connO._error:
    print connO._error, connO._errMsg
  
  print "avgCoordNumber:", connO.avgCoordNumber()
  print "avgCoordNumberBySpecie:", connO.avgCoordNumberBySpecie()
  print "avgBondDistance:", connO.avgBondDistance()
  
  print "Finished."