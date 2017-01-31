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