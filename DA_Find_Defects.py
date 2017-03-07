"""
A script to identify defects by comparing initial and final systems

@author Tomas Lazauskas, 2017
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

from optparse import OptionParser
import numpy as np

import source.IO as IO

def check_systems(input_system, final_system):
  """
  Checks whether the analysis can be performed.
  
  """
  
  success = True
  error = ""
  
  # checking whether the cell dimensions match
  if not np.array_equal(input_system.cellDims, final_system.cellDims):
    success = False
    error = "Input and final systems do not have matching cell dimensions"
    
    return success, error
  
  # checking if the number of atoms matches
  if input_system.NAtoms != final_system.NAtoms:
    success = False
    error = "Input and final systems do not have identical number of atoms"
  
  return success, error

def cmdLineArgs():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog inputFile outputFile radius"
  
  parser = OptionParser(usage=usage)

  parser.disable_interspersed_args()
      
  (options, args) = parser.parse_args()

  if (len(args) != 3):
    parser.error("incorrect number of arguments")

  return options, args

def read_systems(input_file, final_file):
  """
  Reads in input and final systems
  
  """
  
  success = True
  error = ""
  
  input_system = None
  final_system = None
  
  # reading in the initial system
  input_system, error = IO.readSystemFromFile(input_file)
  
  if len(error) > 0:
    success = False
  
  if success:
    # reading in the final system
    final_system, error = IO.readSystemFromFile(final_file)
    
    if len(error) > 0:
      success = False
  
  return success, error, input_system, final_system

if __name__ == "__main__":
  
  # command line arguments and options
  options, args = cmdLineArgs()
  
  # reading in the systems
  success, error, input_system, final_system = read_systems(args[0], args[1])
  
  if success:
    success, error = check_systems(input_system, final_system)
  
  if success:
    vac_radius = np.float(args[2])
  
  if success:
    print "Finished!"
  else:
    print "ERROR: %s" % (error)
   