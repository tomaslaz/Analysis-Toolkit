#!/usr/bin/env python

"""
A script to convert files. Works with XYZ, CAR, and GIN formats.

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import glob
import os
import sys
import numpy as np
from optparse import OptionParser

import source.Atoms as Atoms
import source.IO as IO

def cmdLineArgs():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog inputFile outputFile"
  
  parser = OptionParser(usage=usage)

  parser.disable_interspersed_args()
      
  (options, args) = parser.parse_args()

  if (len(args) < 2 and len(args) > 3):
    parser.error("incorrect number of arguments")

  return options, args

def convertFile(cluster, outFile, controlFile=None):
  """
  Saves a cluster in to a file according to the file format
  
  """
  error = ""
  success = True
  
  if outFile.endswith(".xyz"):
    success, error = IO.writeXYZ(cluster, outFile)
    
  elif outFile.endswith(".car"):
    success, error = IO.writeCAR(cluster, outFile)
    
  elif outFile.endswith(".gin"):
    success, error = IO.writeGIN(cluster, outFile, controlFile=controlFile)
    
  else:
    error = "Undefined output format"
    success = False
    
  return success, error
    
if __name__ == "__main__":
  
  ok = 1
  error = ""
  
  _, args = cmdLineArgs()
  
  inFileName = args[0]
  outFileName = args[1]
  
  if (len(args) > 2):
    controlFile = args[2]
  else:
    controlFile = None
    
  if inFileName[-3:] == outFileName[-3:]:
    sys.exit("Formats must differ!")
  
  if inFileName == ".gin":
    fileList = glob.glob("*.gin")
      
    for fileName in fileList:
      cluster, error = IO.readSystemFromFileGIN(fileName)
        
      outputFile = fileName[:-3] + outFileName[-3:]
      
      ok, error = convertFile(cluster, outputFile, controlFile)
  
  elif inFileName == ".car":
    fileList = glob.glob("*.car")
    
    for fileName in fileList:
      
      cluster = IO.readSystemFromFileCAR(fileName)
       
      outputFile = fileName[:-3] + outFileName[-3:]
      
      ok, error = convertFile(cluster, outputFile, controlFile)
  
  elif inFileName == ".xyz":
    fileList = glob.glob("*.xyz")
    
    for fileName in fileList:
      cluster = IO.readSystemFromFileXYZ(fileName)
        
      outputFile = fileName[:-3] + outFileName[-3:]
      
      ok, error = convertFile(cluster, outputFile, controlFile)
    
  else:
    
    if inFileName.endswith(".xyz"):
      cluster = IO.readSystemFromFileXYZ(inFileName)
      
    elif inFileName.endswith(".car"):
      cluster = IO.readSystemFromFileCAR(inFileName)
    
    elif inFileName.endswith(".gin"):
      cluster, error = IO.readSystemFromFileGIN(inFileName)
    
    else:
      print "Unrecognised input file ", inFileName, " format"
      sys.exit()
      
    ok, error = convertFile(cluster, outFileName, controlFile)
    
    if ok:
      print "Finished!"
    
    else:
      print "Error: ", error
      