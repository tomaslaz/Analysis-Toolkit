"""
Utilities module.

@author Tomas Lazauskas, 2017
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

import copy
import os
import random
import string
import subprocess

try:
  import vtk
  vtk_imported = True
except:
  vtk_imported = False

_systems_stats_file = "Stats.csv"

def delaunay3DArea(system):
  """
  Estimates system's area by summing all the triangles from Delaunay's triangulation
  
  """
  
  delaunayAlpha = 3.0
  delaunayTolerance = 2.9
  
  clusterPoints = vtk.vtkPoints()
    
  for i in range(system.NAtoms):
    clusterPoints.InsertPoint(i, system.pos[i*3 + 0], system.pos[i*3 + 1], system.pos[i*3 + 2])
  
  polyCluster = vtk.vtkPolyData()
  polyCluster.SetPoints(clusterPoints)
  
  delaunayCluster = vtk.vtkDelaunay3D()
  delaunayCluster.SetInputData(polyCluster)
#   delaunayCluster.SetTolerance(delaunayTolerance)
  delaunayCluster.SetAlpha(delaunayAlpha)
  delaunayCluster.BoundingTriangulationOff()
  delaunayCluster.Update()
  
  clusterSurfaceFilter = vtk.vtkGeometryFilter()
  clusterSurfaceFilter.SetInputConnection(delaunayCluster.GetOutputPort())
  clusterSurfaceFilter.Update()
  
  mapClipper = vtk.vtkDataSetMapper()
  mapClipper.SetInputConnection(clusterSurfaceFilter.GetOutputPort())
  
  surfaceArea = 0.0
    
  clipperPointInput = mapClipper.GetInput()  
  for i in range(clipperPointInput.GetNumberOfCells()):
    p1 = clipperPointInput.GetCell(i).GetPoints().GetPoint(0)
    p2 = clipperPointInput.GetCell(i).GetPoints().GetPoint(1)
    p3 = clipperPointInput.GetCell(i).GetPoints().GetPoint(2)
    
    triangleArea = vtk.vtkTriangle.TriangleArea(p1, p2, p3)
    surfaceArea += triangleArea
   
  return surfaceArea

def distanceSq(pos1x, pos1y, pos1z, pos2x, pos2y, pos2z):
  """
  Returns a squared distance between two positions in the Cartesian system 
  
  """
  
  return ((pos1x - pos2x)**2.0 + (pos1y - pos2y)**2.0 + (pos1z - pos2z)**2.0)

def atomicSeparation2(atomPos1, atomPos2, cellDims, PBC):
  """
  Return atomic separation squared with accounted periodic boundary conditions
  
  """
  
  rx = atomPos1[0] - atomPos2[0]
  ry = atomPos1[1] - atomPos2[1]
  rz = atomPos1[2] - atomPos2[2]
  
  if (PBC[0] == 1):
    rx = rx - round( rx / cellDims[0] ) * cellDims[0]

  if (PBC[1] == 1):
    ry = ry - round( ry / cellDims[1] ) * cellDims[1]

  if (PBC[2] == 1):
    rz = rz - round( rz / cellDims[2] ) * cellDims[2]

  sep2 = rx * rx + ry * ry + rz * rz;
      
  return sep2

def get_random_name(n=10):
  """
  Generates a n length random string
  
  """

  return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(n))

def getPathToDreadnaut():
  """
  Return a path to dreadnaut executable.
  
  """
  #nauty25r9
  return os.path.join(os.path.abspath(os.path.dirname(__file__)+"/../"), "thirdparty", "nauty25r9", "dreadnaut")

def run_sub_process(command, verbose=0):
  """
  Run command using subprocess module.
  Return tuple containing STDOUT, STDERR, STATUS
  Caller can decide what to do if status is true
  
  """
  
  if verbose:
    print command
  
  process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output, stderr = process.communicate()
  status = process.poll()
  
  return output, stderr, status

def sort_systems(systems_list):
  """
  Sorts systems according to their energy
  
  """
  
  systems_list_len = len(systems_list)

  for i in range(systems_list_len):
    for j in range(systems_list_len):
      if systems_list[i].totalEnergy < systems_list[j].totalEnergy:
        temp = copy.deepcopy(systems_list[i])
        
        systems_list[i] = copy.deepcopy(systems_list[j])
        systems_list[j] = copy.deepcopy(temp)
    
    if (i % 100 == 0): print "Sorting %d/%d" % (i+1, systems_list_len)

def systems_statistics(systems_list, dir_path=None):
  """
  Generates statistics about the FHI-aims simulations
    
  """
  
  if dir_path is None:
    f = open(_systems_stats_file, "w")
  else:
    f = open(os.path.join(dir_path, _systems_stats_file), "w")
  
  f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % ("System", "Method", "Energy", "Hashkey", "Cores", 
          "Time", "Tot.Time", "H-L", 
          "VBM", "VBMOcc", "VBMSpinChannel", 
          "CBM", "CBMOcc", "CBMSpinChannel", 
          "SpinN", "SpinS", "SpinJ","Size",
          "DimX", "DimY", "DimZ"))
  
  for system in systems_list:
    f.write("%s,%s,%f,%s,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n" % (system.name, system.energyDefinition,
                                     system.totalEnergy, system.hashkey,
                                     system.noOfcores, system.runTime, system.noOfcores*system.runTime,
                                     system.homo_lumo_gap, system.vbm, system.vbm_occ_num, system.vbm_spin_chan, 
                                     system.cbm, system.cbm_occ_num, system.cbm_spin_chan, 
                                     system.spin_N, system.spin_S, system.spin_J, float(system.NAtoms),
                                     system.cellDims[0], system.cellDims[1], system.cellDims[2]))
  f.close()
  
def runDreadnaut(graphString):
   """
   Runs dreadnaut on the given graph
   
   """
   
   if graphString is None:
       return None
   
   dreadnaut = getPathToDreadnaut()

   command = "%s << %s" % (dreadnaut, graphString)

   output, stderr, status = run_sub_process(command)
   
   if status:
       print "WARNING: dreadnaut FAILED"
       print stderr
       return None
   
   array = output.strip().split("\n")
   line = array[-1].strip()
       
   if line[:1] == "[" and line[-1:] == "]":
       return modifyLineHashkey(line)
   
   else:
       print "WARNING: dreadnaut FAILED"
       return None
     
def modifyLineHashkey(line):
    """
    Modifies the hashkey line
    
    """
    
    line = string.replace(line, " ", "_")
    line = string.replace(line, "[", "")
    line = string.replace(line, "]", "")
    
    return line