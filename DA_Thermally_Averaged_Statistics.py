"""
A script to calculate thermally averaged statistics of a property

@author Tomas Lazauskas, 2017
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

from optparse import OptionParser
import math
import numpy as np

# Analysis toolkit modules
import source.Constants as Constants
import source.IO as IO
import DM_DOS as DOS

def calculate_thermaly_averaged_statistics(data_array, temps, unique_mode):
  """
  Calculates canonical potential
  
  """
  
  energies = data_array[:, 0]
  property = data_array[:, 1]
  
  occurences = None
  if unique_mode:
    occurences = data_array[:, 2]
  
  e_min = np.min(energies)
  e_cnt = len(energies)
  
  print "-" * 100
  
  for temparature in temps:
     
    kT = Constants.kB * temparature
     
    sum = 0.0
    sum2 = 0.0
    
    for i in range(e_cnt):
      e_diff = energies[i] - e_min
      
      expr = math.exp(-1.0*(e_diff) / (kT))
      
      if unique_mode:
        expr *= occurences[i]
       
      sum2 += expr
     
      sum += property[i] * expr
              
    thermally_averaged = sum / sum2
    
    print "Thermally averaged (%8.2f K) value: %.4f" % (temparature, thermally_averaged)
       
  print "-" * 100

def cmd_line_args():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog [options] input_file"
  
  parser = OptionParser(usage=usage)
  
  parser.add_option('-t', dest="temps", default="293", help="List of temperatures, separated by a comma (default t=293)")
  
  parser.add_option("-u", dest="unique", default=False, action="store_true",
    help="A flag to say whether the input file contains only unique values. Default = False")
  
  parser.disable_interspersed_args()
      
  (options, args) = parser.parse_args()

  if (len(args) != 1):
    parser.error("incorrect number of arguments")

  return options, args

if __name__ == "__main__":
  
  success = True
  error = ""
  
  # command line arguments and options
  options, args = cmd_line_args()
  
  # getting the temperatures
  if options.temps is not None:
    temps = DOS.getTheListOfTemps(options.temps)
  
  else:
    temps = [0]
  
  data_file = args[0]
  
  if not IO.checkFile(data_file):
    success = False
    error = "Cannot find file: %s" % (data_file)
  
  if success:
    
    # reading in the data file
    data_array = np.loadtxt(data_file, delimiter=',')
    
    # calculating the statistics
    calculate_thermaly_averaged_statistics(data_array, temps, options.unique)
      
  if success:  
    print "Finished!"
    
  else:
    print "ERROR: %s" % (error)