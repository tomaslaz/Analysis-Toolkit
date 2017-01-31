"""
A script to perform unittests for the analysis tool

@author Tomas Lazauskas, 2017
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import os
import sys
import unittest

import source.IO as IO

_available_tests = ["DM_Surface_Energy"]

class Test_DM_Surface_Energy(unittest.TestCase):
  """
  DM_Surface_Energy unittest class
  
  """
  
  @classmethod
  def setUpClass(self):
    """
    Set ups the environment for testing
    
    """
    
    file_name = "unittests/DM_Surface_Energy/Ti_n13.xyz"
    
    self.cwd = os.getcwd()
    
    self.system = IO.readSystemFromFileXYZ(file_name)
  
  def test_surface_energy_calculation(self):
    """
    Testing coordinate transformation
    """
      
    self.assertEqual(1, 1)

def perform_unit_tests(analysis_tool):
  """
  Executes unit tests for the specified analysis_tool
  """  
  
  test_class_name = "Test_" + analysis_tool.strip()
    
  try:
    suite = unittest.TestLoader().loadTestsFromTestCase(eval(test_class_name))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
  
  except:
    return False

if __name__ == "__main__":
  """
  Performs the unit tests.
  
  """

  if len(sys.argv) > 1:
    
    # performunittests for a single module
    if sys.argv[1] in _available_tests:
      perform_unit_tests(sys.argv[1])
    
    # perform all unittests
    elif sys.argv[1] == "all":
      
      for tool in _available_tests:
        perform_unit_tests(tool)
    
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