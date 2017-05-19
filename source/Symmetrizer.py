#!/usr/bin/env python

#import jpype
import os
import jpype
from jpype import *

import IO
import Utilities

# import sys
#import constants

__pointGroupFound = "+ POINT GROUP FOUND:"
__pointGroupStLine = "1."

__tolerance = "0.001"
__verbose = "0"



class Symmetrizer:
  """
  A class to control JAVA environment and launch symmetry identification code.
  
  """
  
  def __del__(self):
    """
    Destructor
    """

    self.__shutdownJVM()
  
  def __init__(self, options=None, args=None, standAlone=True):
    """
    Constructor
    """
    
    self._globJVMOn = False
    self._globJVMSymm = None
    
    randomName = utilities.getRandomName(8)
    
    self._globJVMOutPath = "jvm%s.out" % (randomName)
    self._globJVMErrPath = "jvm%s.err" % (randomName)
    
    self._options = options
    self._args = args
    
    if standAlone:
      self._input = self._args[0]
    
      if self._options.aPointG: 
        self.__startJVM()
    
    else:
      self._input = None
      
      self.__startJVM()
  
  def __startJVM(self):
    """
    Initializes Java Virtual Machine and loads symmetrizer
    """
    
    jarpath = os.path.join(constants._sourceDir, constants._symmetrizerPath)
    
    try:
      jpype.startJVM(jpype.getDefaultJVMPath(), "-ea", "-Djava.class.path=%s" % (jarpath))
      
      self._globJVMOn = True
      self._globJVMSymm = JPackage('net.webmo.symmetry').Main
      messages.log(__name__, "JVM initialised!", verbose=1)
      
      #TODO: In future this should be done via streams
      fs = jpype.JClass("java.io.File")
      ps = jpype.JClass("java.io.PrintStream")
      jpype.java.lang.System.setOut(ps(fs(self._globJVMOutPath)))
      jpype.java.lang.System.setErr(ps(fs(self._globJVMErrPath)))

    except:
      self._globJVMOn = False
      self._globJVMSymm = None
      messages.warning(__name__, "Cannot start JVM!", verbose=0)
          
  def __shutdownJVM(self):
    """
    Shuts down the Java Virtual Machine
    """
    
    if self._globJVMOn:
      jpype.shutdownJVM()
      
      self._globJVMOn = False
      self._globJVMSymm = None
      messages.log(__name__, "JVM shut down!", verbose=1)
      
      self._clearJVMOutputFiles()
  
  def _clearJVMOutputFiles(self):
    """
    Deletes JVM output files.
    """
    
    try:
      io.removeFile(self._globJVMOutPath)
      io.removeFile(self._globJVMErrPath)
      
    except:
      messages.warning(__name__, "Cannot remove JVM output files!", verbose=2)
  
     
  def _getSymmetry(self, clusterO):
    """
    Evaluates the symmetry (point group) of a cluster.
    
    """
    
    return analysis.getSymmetry(self, clusterO)
      
  @timer.timeit
  def run(self):
    """
    The main run routine
    """
        
    self._clearJVMOutputFiles()

def getSymmetry(cluster):
  """
  Runs symmetrizer and tries to get the symmetry from the output.
  
  """
  
  success = True
  error = ""
  symmetry = ""
  
  if not cluspy._globJVMOn:
    success = False
    error = "JVM is has not been initialized"
    
    return success, error, symmetry
  
  # Saves a cluster in to an xyz file.
  tempClusterFile = "%s.xyz" % utilities.getRandomName(10)
  
  success, error = io.writeXYZ(cluster, tempClusterFile)
  
  if not success:
    return success, error, symmetry
  
  # Runs the symmetrizer to get the point group
  try:
    cluspy._globJVMSymm.main(["-v", __verbose, "-t", __tolerance, tempClusterFile])
    
  except:
    #TODO: catch the Symmetrizer error
    success = False
    error = "error while running Symmetrizer"
    
    return success, error, symmetry
  
  io.removeFile(tempClusterFile)
  
  #TODO: In future this should be done via streams
  success, error, symmetry = _getSymmetryFromFile(cluspy._globJVMOutPath, cluspy._globJVMErrPath)
  
  return success, error, symmetry

def _getSymmetryFromFile(outFilePath, errFilePath):
  """
  Reads the symmetry value from the symmetrizer output. At the end clears the files.
  
  """
    
  success = True
  error = ""
  symmetry = ""
  
  # reading the symmetrizer output and looking for symmetry
  if not os.path.isfile(outFilePath):
    print "File [%s] doesn't exist." % (outFilePath)
    return success, error, symmetry
  
  try:
    f = open(outFilePath)
  except:
    print "Cannot read file [%s]" % (outFilePath)
    return success, error, symmetry
  
  i = 0
  
  startLooking = 0
  symmetryFound = 0
  
  for line in f:
    
    i += 1
    line = line.strip()
    
    if not startLooking and __pointGroupFound in line:
      startLooking = 1
    
    if startLooking:
      if __pointGroupStLine in line:
        
        lineArray = line.split()

        if len(lineArray) > 2:
          symArray = lineArray[1].split(":")
          
          symmetry = symArray[0]
          symmetryFound = 1
          
        break
    
  f.close()
  
  if not symmetryFound:
    success = False
    error = "Could not find the symmetry from symmetrizer"
  
  # Removing the output
  try:
    open(outFilePath, 'w').close()
    open(errFilePath, 'w').close()
    
  except:
    f = open(outFilePath, 'r+')
    f.truncate()
    
    f = open(errFilePath, 'r+')
    f.truncate()

  return success, error, symmetry
  