#!/usr/bin/env python

#import jpype
import os
import jpype
from jpype import *

import Constants
import IO
import Messages
import timer
import Utilities

# import sys
#import constants

_pointGroupFound = "+ POINT GROUP FOUND:"
_pointGroupStLine = "1."

_tolerance = "0.01"
_verbose = "0"

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
    
    self.cwd = os.getcwd()
    
    self._globJVMOn = False
    self._globJVMSymm = None
    
    randomName = Utilities.get_random_name(8)
    
    self._globJVMOutPath =  os.path.join(self.cwd, "jvm%s.out" % (randomName))
    self._globJVMErrPath =  os.path.join(self.cwd, "jvm%s.err" % (randomName))
    
    self._options = options
    self._args = args
    
    self._input = None
    
    self.__startJVM()
  
  def __startJVM(self):
    """
    Initializes Java Virtual Machine and loads symmetrizer
    """
    
    jarpath = os.path.join(Constants._sourceDir, Constants._symmetrizerPath)
    
    try:
      jpype.startJVM(jpype.getDefaultJVMPath(), "-ea", "-Djava.class.path=%s" % (jarpath))
      
      self._globJVMOn = True
      self._globJVMSymm = JPackage('net.webmo.symmetry').Main
      Messages.log(__name__, "JVM initialised!", verbose=1)
            
      #TODO: In future this should be done via streams
      fs = jpype.JClass("java.io.File")
      ps = jpype.JClass("java.io.PrintStream")
      jpype.java.lang.System.setOut(ps(fs(self._globJVMOutPath)))
      jpype.java.lang.System.setErr(ps(fs(self._globJVMErrPath)))

    except:
      self._globJVMOn = False
      self._globJVMSymm = None
      Messages.warning(__name__, "Cannot start JVM!", verbose=0)
          
  def __shutdownJVM(self):
    """
    Shuts down the Java Virtual Machine
    """
    
    if self._globJVMOn:
      jpype.shutdownJVM()
      
      self._globJVMOn = False
      self._globJVMSymm = None
      Messages.log(__name__, "JVM shut down!", verbose=1)
      
      self._deleteJVMOutputFiles()
  
  def _deleteJVMOutputFiles(self):
    """
    Deletes JVM output files.
    """

    try:
      IO.removeFile(self._globJVMOutPath)
      IO.removeFile(self._globJVMErrPath)
      
    except:
      Messages.warning(__name__, "Cannot remove JVM output files!", verbose=2)
  
  def _clearJVMOutputFiles(self):
    """
    Clears JVM output files.
    """

    try:
      open(self._globJVMOutPath, 'w').close()
      open(self._globJVMErrPath, 'w').close()
      
    except:
      Messages.warning(__name__, "Cannot remove JVM output files!", verbose=2)
  
  def getSymmetry(self, clusterO):
    """
    Evaluates the symmetry (point group) of a cluster.
    
    """
    
    success = True
    error = ""
    symmetry = ""
    
    if not self._globJVMOn:
      success = False
      error = "JVM is has not been initialized"
      
      return success, error, symmetry
    
    # Saves a cluster in to an xyz file.
    cwd = os.getcwd()
    tempClusterFile = "%s/%s.xyz" % (cwd, Utilities.get_random_name(10))
    
    success, error = IO.writeXYZ(clusterO, tempClusterFile, scfDone=False)
    
    if not success:
      return success, error, symmetry
    
    # Runs the symmetrizer to get the point group
    try:
      self._globJVMSymm.main(["-v", _verbose, "-t", _tolerance, tempClusterFile])
      
    except:
      #TODO: catch the Symmetrizer error
      success = False
      error = "error while running Symmetrizer"
       
      return success, error, symmetry
    
    IO.removeFile(tempClusterFile)
    
    #TODO: In future this should be done via streams
    success, error, symmetry = _getSymmetryFromFile(self._globJVMOutPath, self._globJVMErrPath)
    
    self._clearJVMOutputFiles()
    
    return  success, error, symmetry   

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
    
    if not startLooking and _pointGroupFound in line:
      startLooking = 1
    
    if startLooking:
      if _pointGroupStLine in line:
        
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
  