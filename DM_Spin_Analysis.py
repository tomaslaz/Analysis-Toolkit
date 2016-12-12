"""
A script to prepare FHI-aims directories for spin polarisation analysis

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import os
from optparse import OptionParser

import source.IO as IO
import source.Messages as Messages

# a directory where systems in the xyz format and control.in file should be put
_input_directory = "input"
# input file file extension
_input_extension = "xyz"
# a directory to save the prepared input files for the simulations
_output_directory = "output"
# output directory prefix
_output_prefix = "spin_ini_"
# fhi-aims control file name
_aims_control = "control.in"
# fhi-aims geometry file name
_aims_geometry = "geometry.in"
# default_initial_moment keyword in the control.in file
_aims_keyword_def_ini_moment = "default_initial_moment"

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

def prepare_control_file(ini_spin_value):
  """
  Copies the control file and adjusts the "default_initial_moment"
  
  """
  
  aims_control_temp = "_%s" % (_aims_control)
  
  f_cntlr_in = open(_aims_control, "r")
  f_cntlr_out = open(aims_control_temp, "w")
    
  for line in f_cntlr_in:
    
    # change the default initial spin value
    if _aims_keyword_def_ini_moment in line:
      f_cntlr_out.write("%s %s\n" % (_aims_keyword_def_ini_moment, str(ini_spin_value)))
      
    # write the rest of the lines as they are
    else:
      f_cntlr_out.write(line)
  
  f_cntlr_in.close()
  f_cntlr_out.close()
  
  os.system("mv -f %s %s" % (aims_control_temp, _aims_control))
  
def prepare_directories(options, args):
  """
  The main method to prepare directories for the spin calculations
  
  """
  
  # reading the input systems
  systems_list = read_systems(_input_extension)
  
  # preparing directories with adjusted spin
  prepare_spins(systems_list, options)

def prepare_spins(systems_list, options):
  """
  Prepares simulation directories for the read systems
  
  """
  
  main_dir_path = os.getcwd()
  
  IO.checkDirectory(_output_directory, createMd=1)
  os.chdir(_output_directory)
  
  output_dir_path = os.getcwd()
  
  spins_from = options.spinfr
  spins_to = options.spinto
  
  for spin in range(int(spins_from), int(spins_to)+1):
    
    # creates directories for the systems
    spin_dir_name = "%s%s" % (_output_prefix, str(spin))
    IO.checkDirectory(spin_dir_name, createMd=1)
    
    os.chdir(spin_dir_name)
    
    # prepares control.in and geometry.in files
    write_systems(systems_list, spin, main_dir_path)
    
    os.chdir(output_dir_path)
      
  os.chdir(main_dir_path)
  
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

def write_systems(systems_list, spin_file, main_dir_path):
  """
  Writing systems as geometry.in files and copy a control file to the system directory
  
  """
  
  cwd = os.getcwd()
  
  for system in systems_list:
    system_name = system.name
    
    IO.checkDirectory(system_name, createMd=1)
    os.chdir(system_name)
    
    # writes system as a geometry.in file
    _, _ = IO.writeAimsGeometry(system, _aims_geometry)
    
    # copies the control.in file
    cmd_line = "cp %s/%s/%s . " % (main_dir_path, _input_directory, _aims_control)
    os.system(cmd_line)
    
    # adjusts the control.in file
    prepare_control_file(spin_file)    
    
    os.chdir(cwd)
    
if __name__ == "__main__":
  options, args = cmd_line_args()
  
  prepare_directories(options, args)
  
  Messages.printAuthor()
  