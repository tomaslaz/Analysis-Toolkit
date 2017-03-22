"""
Finds unique structures by examining the working directory (recursively) and reading in systems with the predefined extension. Then examines them using the hashkey and creates a directory with all the unique systems and a statistics file.

@author Tomas Lazauskas, 2017
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

from optparse import OptionParser

# Analysis toolkit modules
import source.IO as IO
import source.Utilities as Utilities

def cmd_line_args():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog extension"
  
  parser = OptionParser(usage=usage)
  
  parser.disable_interspersed_args()
  
  (options, args) = parser.parse_args()
    
  if (len(args) != 1):
    parser.error("incorrect number of arguments")

  return options, args

if __name__ == "__main__":
  
  # reading in input arguments
  options, args = cmd_line_args()
  
  # finding the structures (paths)
  systems_files_list = IO.get_file_list_recursive(args[0])
  
  # reading in the systems
  systems = IO.read_in_systems(systems_files_list)
  
  # getting the hashkeys
  unique_systems = IO.get_unique_systems_hashkeys(systems)
  
  # finding unique systems
  
  
  
  # sorting the systems
  #Utilities.sort_systems(systems)
  
  # 