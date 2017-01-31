"""
A script to calculate surface energy for a cluster.

@author Tomas Lazauskas, 2017
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import os
import sys

from optparse import OptionParser

from source.Messages import log

# verbose level: 0 - off, 1 - on
_verbose = 1
# full path to the compiled adul v2
_adul_v2_path = "/Users/Tomas/Software/ARVO/arvo_c/arvo_c"

def cmd_line_args():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog input_file radius"
  
  parser = OptionParser(usage=usage)
  
  parser.disable_interspersed_args()
  
  (options, args) = parser.parse_args()
    
  if (len(args) != 2):
    parser.error("incorrect number of arguments")

  return options, args

if __name__ == "__main__":
  
  # reading in input arguments
  options, args = cmd_line_args()
  
  file_name = args[0]
  radius = args[1]
  
  message = "reading file: %s" % (file_name)
  log(__name__, message, verbose=_verbose)