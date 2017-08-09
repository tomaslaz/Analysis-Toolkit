"""
A script to compile the thirdparty software for the analysis tools

@author Tomas Lazauskas, 2017
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import os
import sys
import glob

c_flags_set = False
CC = ""
CCFLAGS = ""
CCLINKFLAGS = ""
INCLUDES = ""

# Operating system: MacOS
if os.uname()[0] == "Darwin":

    CC = "gcc "
    CCFLAGS = "-c -O3 -fPIC "
    #CCLINKFLAGS = "-lm"  
    CCLINKFLAGS = "-bundle -flat_namespace -undefined suppress "  
    c_flags_set = True


# if os.uname()[0] == "Darwin":
#   # use gcc
#   CC = "gcc "
#   CCFLAGS = "-c -O3 -fPIC "
#   CCLINKFLAGS = "-bundle -flat_namespace -undefined suppress "  
#   INCLUDES = ""

# Linux box?
elif os.uname()[0] == "Linux":

  CC = "gcc "
  CCFLAGS = "-c -O3 -fPIC "
  CCLINKFLAGS = "-shared "  
  INCLUDES = ""
  
  c_flags_set = True

###################################################
# NOTHING BELOW THIS LINE SHOULD NEED TO BE CHANGED
###################################################
if os.path.dirname(__file__):
    os.chdir(os.path.dirname(__file__))

if not c_flags_set:
  sys.exit("C compilation flags has not been set")

cwd = os.getcwd()

os.chdir("source")
os.chdir("c_libs")

# MAKE SURE NO WRAP C LIBS LEFT OVER
wrap_c_files = glob.glob("*_wrap.c")
if len(wrap_c_files):
  for fn in wrap_c_files:
    os.unlink(fn)

CLIBSRC = [fn[:-2] + "c" for fn in glob.glob("*[!__init__].py") if not fn.endswith("_utils.py")]

DEPS = {}

SOBJ = {}
OBJ = {}
for foo in CLIBSRC:
  SOBJ[foo] = "_%s.so" % foo[:-2]
  OBJ[foo] = "%s.o" % foo[:-2]

EXTC = []
tmpc = glob.glob("*.c")
for foo in tmpc:
  match = False
  if foo in CLIBSRC:
    match = True
  
  if match:
    pass
  
  else:
    EXTC.append(foo)

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
  
  os.chdir(cwd)
  os.chdir("thirdparty")
  os.chdir("nauty25r9")
  
  command = "make clean"
  run_command(command)
    
  command =  "rm -f dreadnaut makefile nauty.h"
  run_command(command)
  
  os.chdir(cwd)
  os.chdir("source")
  os.chdir("c_libs")
  
  command = "rm -f *.pyc *.o *.so *_wrap.c "
  run_command(command)
  
def main():
  print "COMPILING"
  print "-------------"
  
  # c arvo
  os.chdir(cwd)
  
  os.chdir("thirdparty")
  os.chdir("arvo_c")
  
  command = "%s %s -o %s %s %s " % (CC, "", "arvo_c", "arvo_c.c", "")
  run_command(command)
  
  os.chdir(cwd)
  os.chdir("thirdparty")
  
  os.chdir("nauty25r9")
  print "-------------"
  
  command = "./configure && make nauty"
  run_command(command)
  
  print "-------------"
  
  # c libraries
  os.chdir(cwd)
  
  os.chdir("source")
  os.chdir("c_libs")
  
   # first compile non swig C libraries
  for source in EXTC:
    command = "%s %s %s %s" % (CC, CCFLAGS, INCLUDES, source)
    run_command(command)

  print "-------------"
  
  # compile non swig c files
  for src in CLIBSRC:
    # compile src file
    command = "%s %s %s" % (CC, CCFLAGS, src)
    run_command(command)
    
    print "-------------"
  
  # link non swig libraries
  for src in CLIBSRC:
    try:
      EXTO = DEPS[src]
        
    except KeyError:
      EXTO = ""
        
    command = "%s %s %s %s -o %s" % (CC, CCLINKFLAGS, OBJ[src], EXTO, SOBJ[src])
    run_command(command)
    
    print "-------------"
  
      
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
  
 