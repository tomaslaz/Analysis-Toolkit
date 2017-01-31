"""
A script to prepare/execute/analyse FHI-aims polarisation calculations

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import os
from optparse import OptionParser

import source.IO as IO
import source.Messages as Messages
from source.Messages import log

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
# fhiaims execution command (as an example)
_aims_exe_cmd = "source /opt/intel/composer_xe_2015/mkl/bin/mklvars.sh intel64; mpirun -n 8 /Users/Tomas/Software/fhi-aims.160328/bin/aims.160328_1.mpi.x > fhiaims.out"

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
  
  parser.add_option("-p", "--prepare", dest="prepare", action="store_true", default=False, 
    help="Prepares directories for the spin calculations")
  
  parser.add_option("-x", "--execute", dest="execute", action="store_true", default=False, 
    help="Executes the simulations")
  
#   parser.add_option("-l", "--analyse", dest="analyse", action="store_true", default=False, 
#     help="Analyses the simulations")
  
  parser.disable_interspersed_args()
      
  (options, args) = parser.parse_args()
  
  return options, args

def execute():
  """
  Executes the fhi-aims calculations
  
  """
  
  main_dir_path = os.getcwd()
  
  Messages.log(__name__, "Running the simulations in: %s" % (_output_directory))
  
  dir_list = IO.get_dir_list(_aims_geometry)
  
  cwd = os.getcwd()
  
  dir_tot_str = str(len(dir_list))
  dir_cnt = 1
  
  for dir_path in dir_list:
    
    Messages.log(__name__, "%s/%s Executing FHI-aims in: %s " % (str(dir_cnt), dir_tot_str, dir_path), 1)
    
    os.chdir(dir_path)
    
    os.system(_aims_exe_cmd)
    
    os.chdir(cwd)
    
    dir_cnt += 1
  
  Messages.log(__name__, "Finished executing FHI-aims calculations")
    
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
  
def prepare_directories(options):
  """
  The main method to prepare directories for the spin calculations
  
  """
  
  Messages.log(__name__, "Preparing directories for the spin calculations")
  
  # reading the input systems
  systems_list = read_systems(_input_extension)
  
  # preparing directories with adjusted spin
  prepare_spins(systems_list, options)

def prepare_spins(systems_list, options):
  """
  Prepares simulation directories for the read systems
  
  """
    
  main_dir_path = os.getcwd()
  
  Messages.log(__name__, "Preparing the simulation files. The files will be saved in: %s" % (_output_directory))
  
  IO.checkDirectory(_output_directory, createMd=1)
  os.chdir(_output_directory)
  
  output_dir_path = os.getcwd()
  
  spins_from = options.spinfr
  spins_to = options.spinto
  
  for spin in range(int(spins_from), int(spins_to)+1):
    
    Messages.log(__name__, "Preparing: %s %s" % (_aims_keyword_def_ini_moment, str(spin)), 1)
    
    # creates directories for the systems
    spin_dir_name = "%s%s" % (_output_prefix, str(spin))
    IO.checkDirectory(spin_dir_name, createMd=1)
    
    os.chdir(spin_dir_name)
    
    # prepares control.in and geometry.in files
    write_systems(systems_list, spin, main_dir_path)
    
    os.chdir(output_dir_path)
      
  os.chdir(main_dir_path)
  
def read_systems(file_extension):
  """
  Reading the input files as system objects
  
  """
  
  Messages.log(__name__, "Reading in the systems from: %s" % (_input_directory), 1)
  
  systems_list = []
  
  cwd = os.getcwd()
  os.chdir(_input_directory)
  
  files_list = IO.get_file_list(file_extension)
  
  for file_name in files_list:
    system = IO.readSystemFromFileXYZ(file_name)
    
    systems_list.append(system)
  
  os.chdir(cwd)
  
  Messages.log(__name__, "Read %s systems." % (str(len(systems_list))), 1)
  
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
  
  Messages.log(__name__, "Running..")
  
  # reading the command line arguments and options
  options, _ = cmd_line_args()
  
  # prepare the directories?
  if options.prepare:
    prepare_directories(options)
    
  # execute the simulations?
  if options.execute:
    execute()
  
  Messages.log(__name__, "Finished.")
  Messages.printAuthor()
  