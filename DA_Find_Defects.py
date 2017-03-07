"""
A script to identify defects by comparing initial and final systems

@author Tomas Lazauskas, 2017
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

from optparse import OptionParser

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

if __name__ == "__main__":
  
  options, args = cmdLineArgs()