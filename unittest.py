"""
A script to perform unittests for the analysis tool

@author Tomas Lazauskas, 2017
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import sys

_available_tests = ["DM_Surface_Energy"]

if __name__ == "__main__":
  """
  Performs the unit tests.
  
  """

  if len(sys.argv) > 1:
    
    # performunittests for a single module
    if sys.argv[1] in _available_tests:
      
      print sys.argv
    
    # perform all unittests
    elif sys.argv[1] == "all":
      
      print sys.argv
    
    else:
      print "Usage: unittest.py [options]\n\nOptions are:\n    [module]"

      for i in range(len(_available_tests)):
        print "        ", _available_tests[i]
        
      print "    all"
      sys.exit(8)
    
  else:
    print "Usage: unittest.py [options]\n\nOptions are:\n    [module]"

    for i in range(len(_available_tests)):
      print "        ", _available_tests[i]
      
    print "    all"
    sys.exit(8)