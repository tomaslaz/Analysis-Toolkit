#!/usr/bin/env python
"""
A script to slowly push a system alongside the eigenvector direction.

@author Tomas Lazauskas, 2016
@web lazauskas.net
@email tomas[a]lazauskas.net
"""

import numpy as np
import os
import subprocess

_constStepFrom =  0.0
_constStepTo   =  1.0
_constStepStep =  0.01
_constEigenVectorNumber = 1

_constPathToAims = "/Users/Tomas/Software/fhi-aims.160328_2/bin/aims.160328_1.serial.x"
_constPathToAims = "mpirun -np 16 /home/tomas/bin/aims.071914_2.mpi.x"

_constControlFile = "control.in"
_constGeometryFile = "geometry.in"
_constEigenFile = "0.001.xyz"
_constFHIaimsOutFile = "fhiaims.out"

class Geometry():
  """
  A geometry class 
  
  """
  
  def __init__(self, atomCnt, pos, types):
    """
    Constructor.
    
    """
    
    self.atomCnt = atomCnt
    self.pos = pos
    self.types = types
    self.eigenVector = None
    
    self.stepSize = []
    self.stepEnergy = []
  
  def setEigenVector(self, eigenVector):
    """
    Sets the eigenvector.
    
    """
    
    self.eigenVector = eigenVector
  
  def performAnalysis(self):
    """
    Performs analysis within a certain step range
    
    """
    
    step = _constStepFrom
    while step < _constStepTo+_constStepStep/2:
      
      print "Analysing ", step
      
      self.pushGeometry(step)
      
      step += _constStepStep
  
  def printResults(self):
    """
    Prints the analysis results to the stdout
    
    """
    
    for i in range(len(self.stepSize)):
      print i, self.stepSize[i], self.stepEnergy[i]
    
  def pushGeometry(self, stepSize):
    """
    Pushes the structure in the direction of the eigenvector and evaluates the energy
    
    """
    
    newPos = self.pos + stepSize*self.eigenVector
    
    dirName = "step_%f" % (stepSize)
    
    os.system("rm -rf %s" % dirName)
    os.system("mkdir %s" % dirName)
    
    lcwd = os.getcwd()
    os.chdir(dirName)
        
    self.saveGeometryFile(newPos, self.types, self.atomCnt)
    os.system("cp ../%s ." % (_constControlFile))
    
    runAims()
    energy = readAimsEnergy()
    
    self.stepSize.append(stepSize)
    self.stepEnergy.append(energy)
    
    os.chdir(lcwd)    
  
  def saveGeometryFile(self, pos, types, atomCnt):
    """
    Saves new structure into a file
    
    """
    
    try:
      f = open(_constGeometryFile, "w")
      
    except:
      print "Cannot read file [%s]" % (_constGeometryFile)
    
    for i in range(atomCnt):
      f.write("%s %.10f %.10f %.10f %s\n" % ("atom", pos[i][0], pos[i][1], pos[i][2], types[i]))
      
    f.close()
  
def readinEigenVector(atomCnt, eigenVectorNo):
  """
  Reads in eigevectors from frequencies calculations
  
  """
  
  try:
    f = open(_constEigenFile)
  except:
    print "Cannot read file [%s]" % (_constEigenFile)
  
  # selecting the lines corresponding to the eigenvector
  fromLine = (eigenVectorNo - 1) * (atomCnt + 2) + 2
  toLine =    eigenVectorNo      * (atomCnt + 2)
  
  eigenVector = np.zeros((atomCnt, 3), np.float)
  
  i = 0
  ii = 0
  for line in f:
        
    if i >= fromLine and i < toLine:
      line = line.strip()
      array = line.split()
      
      # first three numbers are coordinates, reading in the next three
      for j in range(3):
        eigenVector[ii][j] = array[4+j]
      
      ii += 1
    i += 1
  
  f.close()
  
  return eigenVector
  
def readinGeometry():
  """
  Reads in geometry.in file and returns a geometry object
  
  """
  
  try:
    f = open(_constGeometryFile)
  except:
    print "Cannot read file [%s]" % (_constGeometryFile)
  
  # getting the number of atoms
  atomCnt = 0  
  for line in f:
    
    line = line.strip()
    array = line.split()
    
    if len(array) == 5:
      atomCnt += 1
  
  # set the pointer to the beginning
  f.seek(0, 0)

  pos = np.zeros((atomCnt, 3), np.float)
  
  dt = np.dtype((str, 2))
  types = np.empty(atomCnt, dt)
  
  # reading the values in
  i = 0  
  for line in f:
    
    line = line.strip()
    array = line.split()
    
    if len(array) == 5:
      for j in range(3):
        pos[i][j] = array[j+1]
      
      types[i] = array[4]
      
      i += 1

  f.close()
    
  return Geometry(atomCnt, pos, types)

def readAimsEnergy():
  """
  Reads the system's energy from the fhiaims output file
  
  """
  
  energy = None

  try:
    fin = open(_constFHIaimsOutFile, "r")
    
  except:    
    return energy
  
  for line in fin:
    
    fields = line.strip().split()
    
    if ((len(fields) > 5) and (' '.join(fields[1:4]) == "Total energy uncorrected")):
      energy = float(fields[5])
            
  fin.close()
  
  return energy

def runAims():
  """
  Executes FHI-aims in the current directory
  
  """
  
  out = open(_constFHIaimsOutFile, "w")
  null = open('/dev/null', "r")
  
  cmd = _constPathToAims
  
  p = subprocess.Popen(cmd, shell=True, stdin=null, stdout=out, stderr=subprocess.STDOUT) 
  p.communicate()

  status = p.poll() 

  out.close()
  null.close()
  
if __name__ == "__main__":
    
  # reads in geometry.in
  geometry = readinGeometry()
  
  # reads in the eigenvector
  geometry.setEigenVector(readinEigenVector(geometry.atomCnt, _constEigenVectorNumber))
  
  # perform analysis by pushing the geometry in the direction of the eigenvector
  geometry.performAnalysis()
  
  # prints the result to the standard Output
  geometry.printResults()
  
  print "Finished."