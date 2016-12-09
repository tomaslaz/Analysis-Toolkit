"""
A script to prepare FHI-aims directories for spin polarisation analysis

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import os
from optparse import OptionParser

import source.IO as IO

# a directory where systems in the xyz format and control.in file should be put
_input_directory = "input"
# input file file extension
_input_extension = "xyz"
# a directory to save the prepared input files for the simulations
_output_directory = "output"
# fhi-aims control file name
_aims_control = "control.in"
# fhi-aims geometry file name
_aims_geometry = "geometry.in"

def cmd_line_args():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog "
  
  parser = OptionParser(usage=usage)
  
  parser.add_option("-a", "--spinfr", dest="spinfr", default=0.0, type="float",
    help="Initial spin values from")
  
  parser.add_option("-b", "--spinto", dest="spinto", default=0.0, type="float",
    help="Initial spinn values to")
      
  parser.disable_interspersed_args()
      
  (options, args) = parser.parse_args()
  
  
  return options, args

def read_systems(input_dir):
  """
  Reading the input files as system objects
  
  """
  
  systems_list = []
  
  cwd = os.getcwd()
  os.chdir(_input_directory)
  
  files_list = IO.get_file_list(input_dir)
  
  for file_name in files_list:
    system = IO.readSystemFromFileXYZ(file_name)
    
    systems_list.append(system)
  
  os.chdir(cwd)
  
  return systems_list

def write_systems(systems_list):
  """
  Writing systems as geometry.in files
  
  """
  
  for system in systems_list:
    _, _ = IO.writeAimsGeometry(system, _aims_geometry)
  
  
if __name__ == "__main__":
  options, args = cmd_line_args()
  
  prepare_directories(options, args)
  
#   systems_list = read_systems(_input_extension)
# 
#   write_systems(systems_list)