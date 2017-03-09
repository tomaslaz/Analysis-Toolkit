"""
Utilities module.

@author Tomas Lazauskas, 2017
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

import random
import string
import subprocess

def distanceSq(pos1x, pos1y, pos1z, pos2x, pos2y, pos2z):
  """
  Returns a squared distance between two positions in the Cartesian system 
  
  """
  
  return ((pos1x - pos2x)**2.0 + (pos1y - pos2y)**2.0 + (pos1z - pos2z)**2.0)

def atomicSeparation2(atomPos1, atomPos2, cellDims, PBC):
  """
  Return atomic separation squared with accounted periodic boundary conditions
  
  """
  
  rx = atomPos1[0] - atomPos2[0]
  ry = atomPos1[1] - atomPos2[1]
  rz = atomPos1[2] - atomPos2[2]
  
  if (PBC[0] == 1):
    rx = rx - round( rx / cellDims[0] ) * cellDims[0]

  if (PBC[1] == 1):
    ry = ry - round( ry / cellDims[1] ) * cellDims[1]

  if (PBC[2] == 1):
    rz = rz - round( rz / cellDims[2] ) * cellDims[2]

  sep2 = rx * rx + ry * ry + rz * rz;
      
  return sep2

def get_random_name(n=10):
  """
  Generates a n length random string
  
  """

  return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(n))

def run_sub_process(command, verbose=0):
  """
  Run command using subprocess module.
  Return tuple containing STDOUT, STDERR, STATUS
  Caller can decide what to do if status is true
  
  """
  
  if verbose:
    print command
  
  process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output, stderr = process.communicate()
  status = process.poll()
  
  return output, stderr, status