#!/usr/bin/env python

"""
A script to convert coordinates.

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

from optparse import OptionParser

def cmd_line_args():
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

if __name__ == "__main__":
  
  ok = 1
  error = ""
  
  _, args = cmd_line_args()
  
  if ok:
    print "Finished!"
  
  else:
    print "Error: ", error