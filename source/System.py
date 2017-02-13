"""
System module.

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

import copy
import math
import numpy as np
import os

import Utilities
from scipy.constants.constants import Rydberg

_const_zero_value = 0.0
_const_def_value = -9999999999.9
_const_path_to_arvo = "thirdparty/arvo_c/arvo_c"

class System(object):
  """
  A class to save the systems.
  
  NAtoms: number of atoms in system (N)
  specie[N]: array of symbols of atoms
  pos[3N]: array of positions of atoms
  charge[N]: array of charges of atoms
  
  """
    
  def __init__(self, NAtoms):
    self.name = ""
    self.hashkey = ""
    
    self.totalEnergy = 0.0
    self.NAtoms = NAtoms
    self.cellDims = np.zeros(3, np.float64)
    self.cellAngles = np.empty(3, np.float64)
    
    self.specie = np.empty(self.NAtoms, np.int32)
    self.pos = np.empty(3*self.NAtoms, np.float64)
    self.minPos = np.empty(3, np.float64)
    self.maxPos = np.empty(3, np.float64)
    
    self.charge = np.empty(self.NAtoms, np.float64)
    
    self.com = np.empty(3, np.float64)
    self.cog = np.empty(3, np.float64)
    self.momentOfInertia = np.zeros([3, 3], np.float64)
    
    dt = np.dtype((str, 5))
    self.specieList = np.empty(0, dt)
    self.specieCount = np.empty(0, np.int32)
    
    self.PBC = np.zeros(3, np.int32)
    
    self.gulpSpecies = {}
    self.gulpAtomType = ["" for x in range(self.NAtoms)]
    self.gulpAtomExtraInfo = ["" for x in range(self.NAtoms)]
    
    self.noOfcores = 0
    self.runTime = 0.0
    
    # spin polarization
    self.spin_N = _const_def_value
    self.spin_S = _const_def_value
    self.spin_J = _const_def_value
    
    # homo-lumo
    self.homo_lumo_gap = _const_def_value
    
    # homo
    self.vbm = _const_def_value 
    
    self.vbm_occ_num = _const_def_value
    self.vbm_spin_chan = _const_def_value
    
    # lumo
    self.cbm = _const_def_value
    
    self.cbm_occ_num = _const_def_value
    self.cbm_spin_chan = _const_def_value
    
    # surface energy
    self.surface_energy = _const_zero_value
    
    # geometrical properties with arvo
    self.arvo_calc = False
    self.arvo_volume = _const_zero_value
    self.arvo_area = _const_zero_value
    self.arvo_spheres = _const_zero_value
    
    # eigenvalues
    self.eigenvalues = None
    self.ev_dos_bins = None
    self.ev_dos = None
    
  def addAtom(self, sym, pos, charge):
    """
    Add an atom to the system
    
    """
    
    if sym not in self.specieList:
        self.addSpecie(sym)
    
    specInd = self.specieIndex(sym)
    
    self.specieCount[specInd] += 1
    
    pos = np.asarray(pos, dtype=np.float64)
    
    self.specie = np.append(self.specie, np.int32(specInd))
    self.pos = np.append(self.pos, pos)
    self.charge = np.append(self.charge, charge)

    self.NAtoms += 1
  
  def calcCOG(self):
      
    """
    Calculates the the geometric centre
    
    """
    
    totMass = 0.0
    self.cog[:] = 0.0
    
    for i in range(self.NAtoms):      
      for j in range(3):
        self.cog[j] += self.pos[3*i + j]

    self.cog = self.cog / self.NAtoms
    
  def calcCOM(self):
      
    """
    Calculates the centre of mass of a system
    
    """
    
    totMass = 0.0
    self.com[:] = 0.0
    
    for i in range(self.NAtoms):
      atomMass = Atoms.atomicMassAMU(self.specieList[self.specie[i]])
      totMass += atomMass
      
      for j in range(3):
        self.com[j] += atomMass * self.pos[3*i + j]

    self.com = self.com / totMass
  
  def calc_ev_dos(self, ev_from=-100, ev_to=20, delta=0.01, sigma=0.1):
    """
    Calculates dos of electronic eigenvalues    
    
    """
    
    success = True
    error = ""
    
    if self.eigenvalues is None:
      success = False
      error = "Eigenvalues were not found"
      
      return success, error
    
    _extraBins=2
    
    # get min and max 
    ev_dos_min = np.float128(ev_from)
    ev_dos_max = np.float128(ev_to)
        
    # number of bins
    ev_dos_n_bins = np.around(int((ev_dos_max - ev_dos_min) / delta) + _extraBins, decimals=0)
        
    # array to hold the ev bin values
    ev_dos_bins = np.arange(ev_dos_min, np.around(ev_dos_min + ev_dos_n_bins * delta, decimals=4), delta)

    # array for the dos values
    ev_dos = np.zeros(ev_dos_n_bins, dtype=np.float128) 
        
    # calculating DOS
    for i in range(ev_dos_n_bins):
      ev_dos[i] = np.sum((1/(sigma*np.pi**0.5)) * np.exp(-(ev_dos_bins[i] - self.eigenvalues)**2 / sigma**2))
    
    self.ev_dos_bins = copy.deepcopy(ev_dos_bins)
    self.ev_dos = copy.deepcopy(ev_dos)
     
    return success, error
  
  def calcMOI(self):
    """
    Calculates moment of inertia
    
    """
    
    moi = np.zeros(6, np.float64)

    
    self.momentOfInertia[:] = 0.0
    
    for i in range(self.NAtoms):
      atomMass = Atoms.atomicMassAMU(self.specieList[self.specie[i]])
        
      moi[0] += (self.pos[3*i+1]**2 + self.pos[3*i+2]**2) * atomMass
      moi[1] += (self.pos[3*i+0]**2 + self.pos[3*i+2]**2) * atomMass
      moi[2] += (self.pos[3*i+0]**2 + self.pos[3*i+1]**2) * atomMass
      moi[3] += -(self.pos[3*i+0] * self.pos[3*i+1]) * atomMass
      moi[4] += -(self.pos[3*i+0] * self.pos[3*i+2]) * atomMass
      moi[5] += -(self.pos[3*i+1] * self.pos[3*i+2]) * atomMass
    
    self.momentOfInertia[0][0] = moi[0]
    self.momentOfInertia[1][1] = moi[1]
    self.momentOfInertia[2][2] = moi[2]
    
    self.momentOfInertia[0][1] = moi[3]
    self.momentOfInertia[0][2] = moi[4]
    
    self.momentOfInertia[1][0] = moi[3]
    self.momentOfInertia[1][2] = moi[5]
    
    self.momentOfInertia[2][0] = moi[4]
    self.momentOfInertia[2][1] = moi[5]
  
  def calc_geo_measures(self, radius):
    """
    Calculates geometrical measures: volume, area
    
    """
    
    _temp_file = Utilities.get_random_name()
    _temp_file = "temp"
    
    # TODO: Need a way set the surface radius for a each atom. Probably add it as 
    #       as a new column in the atoms.in file.
    
    # prepare a temporary file
    File.writeATS(self, _temp_file, radius)
    
    command = "%s protein=%s" % ("./thirdparty/arvo_c/arvo_c", _temp_file)
    output, stderr, status = Utilities.run_sub_process(command)
    
    # if the execution of the was successful:
    if not status:
      
      output_array = output.split()
      
      self.arvo_calc = True
      self.arvo_volume = np.float64(output_array[1])
      self.arvo_area = np.float64(output_array[3])
      self.arvo_spheres = np.int16(output_array[6])
  
    os.unlink(_temp_file)
    
  def calc_surface_energy(self, radius):
    """
    Calculates surface energy using arvo_c thirdparty code
    
    """
    
    # calculate the geometrical measures
    self.calc_geo_measures(radius)
    
    if self.arvo_calc:
    
      self.surface_energy = self.totalEnergy / self.arvo_area
   
  def findNN(self, atomIdx, rdfCutOffSq):
    
    cntrx = self.pos[3*atomIdx+0]
    cntry = self.pos[3*atomIdx+1]
    cntrz = self.pos[3*atomIdx+2]
    
    neighboursCnt = 0
    neighboursArr = np.zeros(self.NAtoms, np.int32)
    neighboursDistArr = np.zeros(self.NAtoms, np.float64)
    
    xdim = self.cellDims[0]
    ydim = self.cellDims[1]
    zdim = self.cellDims[2]
    
    for i in range(self.NAtoms):
      if (i != atomIdx):
        
        # distances
        rx = cntrx - self.pos[3*i+0]
        ry = cntry - self.pos[3*i+1]
        rz = cntrz - self.pos[3*i+2]
        
        # applying cubic periodic boundary conditions
        if self.PBC[0]:
          prev_rx = rx
          rx = rx - np.round( rx / xdim ) * xdim
          
        if self.PBC[1]:
          prev_ry = ry
          ry = ry - np.round( ry / ydim ) * ydim

        if self.PBC[2]:
          prev_rz = rz
          rz = rz - np.round( rz / zdim ) * zdim
        
        distSq = (rx**2 + ry**2 + rz**2)
        
        if (distSq <= rdfCutOffSq):
          neighboursArr[neighboursCnt] = i
          neighboursDistArr[neighboursCnt] = math.sqrt(distSq)
          
          neighboursCnt += 1
    
    return neighboursCnt, neighboursArr, neighboursDistArr
  
  def rotateToMOI(self, basis):
    
    for i in range(self.NAtoms):
      for j in range(3):
        elSum = 0.0
        
        for k in range(3):
          elSum += self.pos[3*i+k] * basis[k][j]

        self.pos[3*i+j] = elSum
  
  def moveToCOG(self):
    """
    Centers the system on the centre of geometry.
    
    """
    
    for i in range(self.NAtoms):
      for j in range(3):
        self.pos[3*i + j] -= self.cog[j]
       
  def moveToCOM(self):
    """
    Centers the system on the centre of mass.
    
    """
    
    for i in range(self.NAtoms):
      for j in range(3):
        self.pos[3*i + j] -= self.com[j]
      
  def removeAtom( self, index ):
    """
    Remove an atom from the system
    
    """
    
    specInd = self.specie[index]
    self.specie = np.delete(self.specie, index)
    self.pos = np.delete(self.pos, [3*index,3*index+1,3*index+2])
    self.charge = np.delete(self.charge, index)

    self.NAtoms -= 1
    
    self.specieCount[specInd] -= 1
    if self.specieCount[specInd] == 0 and not self.specieListForced:
        self.removeSpecie(specInd)
  
  def removeSpecie(self, index):
    """
    Remove a specie from the specie list.
    
    """
    self.specieCount = np.delete(self.specieCount, index)
    self.specieList = np.delete(self.specieList, index)
    
    for i in xrange(self.NAtoms):
        if self.specie[i] > index:
            self.specie[i] -= 1

  def specieIndex(self, check):
    """
    Index of sym in specie list
    
    """
    
    count = 0
    index = -1
    for sym in self.specieList:
        if sym == check:
            index = count
            break
        
        count += 1
    
    return index 
  
  def addSpecie(self, sym, count=None):
    """
    Add specie to specie list
    
    """
            
    if sym in self.specieList:
        if count is not None:
            specInd = self.specieIndex(sym)
            self.specieCount[specInd] = count
        
        return
    
    if count is None:
        count = 0
    
    self.specieList = np.append(self.specieList, sym)
    self.specieCount = np.append(self.specieCount, np.int32(count))
    
  def minMaxPos(self, PBC):
      
    for i in xrange(3):
        if not PBC[i]:
            self.minPos[i] = self.pos[i::3].min()
            self.maxPos[i] = self.pos[i::3].max()