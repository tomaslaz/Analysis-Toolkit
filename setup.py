"""
A script to compile the thirdparty software for the analysis tools

@author Tomas Lazauskas, 2017
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import os
import sys

c_flags_set = False
CC = ""
CCFLAGS = ""
CCLINKFLAGS = ""

# Operating system: MacOS
if os.uname()[0] == "Darwin":

    CC = "gcc "
    CCFLAGS = "-c -O3 -fPIC "
    CCLINKFLAGS = "-lm"  
    
    c_flags_set = True

###################################################
# NOTHING BELOW THIS LINE SHOULD NEED TO BE CHANGED
###################################################
if os.path.dirname(__file__):
    os.chdir(os.path.dirname(__file__))

if not c_flags_set:
  sys.exit("C compilation flags has not been set")

cwd = os.getcwd()

def clean():
  
  print "CLEANING"
  print "-------------"
  
  command = "rm -f *.pyc"
  run_command(command)
  
  os.chdir(cwd)
  
  os.chdir("thirdparty")
  os.chdir("arvo_c")
    
  command = "rm -f arvo_c"
  run_command(command)
    
def main():
  print "COMPILING"
  print "-------------"
  
  os.chdir(cwd)
  
  os.chdir("thirdparty")
  os.chdir("arvo_c")
  
  command = "%s %s -o %s %s %s " % (CC, "", "arvo_c", "arvo_c.c", "")
  run_command(command)
      
def run_command(command, check_status=True):
  
  print "%50s ::  %s" %  (os.getcwd()[-50:], command)
  
  status = os.system(command)
  
  if check_status and status:
    print "COMMAND FAILED"
    sys.exit(25)

if __name__ == '__main__':
  
  if len(sys.argv) > 1:
    if sys.argv[1] == "clean":
      clean()
      sys.exit(0)
    
    elif sys.argv[1] == "make":
      main()
      sys.exit(0)
      
    elif sys.argv[1] == "all":
      clean()
      main()
      sys.exit(0)
      
    else:
      print "Usage: setup.py [options]\n\nOptions are:\n    clean\n    make\n    all"
      sys.exit(8)
      
  else:
    print "Usage: setup.py [options]\n\nOptions are:\n    clean\n    make\n    all"
    sys.exit(8)
  
 