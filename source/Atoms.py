"""
Atoms module.

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com
"""

import sys, os
 
def atomicNumber(sym):
    global atomicNumberDict
     
    try:
        value = atomicNumberDict[sym]
    except:
        sys.exit(__name__+": ERROR: no atomic number for " + sym)
     
    return value
 
def atomicMassAMU(sym):
    global atomicMassDict
     
    try:
        value = atomicMassDict[sym]
    except:
        sys.exit(__name__+": ERROR: no atomic mass for " + sym)
     
    return value
 
def covalentRadius(sym):
    global covalentRadiusDict
     
    try:
        return covalentRadiusDict[sym]
    except:
        sys.exit(__name__+": ERROR: no covalent radius for " + sym)

def getRadius(clusterO, radiusType='ir'):
  
  tempDist = 0.0
  maxDist = 0.0
  
  for sp1 in clusterO.specieList:
    for sp2 in clusterO.specieList:
      
      if sp1 == sp2 and len(clusterO.specieList) > 1:
        continue
      
      if radiusType == 'cr':
        rad1 = covalentRadius(sp1)
        rad2 = covalentRadius(sp2)
        
      elif radiusType == 'ir':
        rad1 = ionicRadius(sp1)
        rad2 = ionicRadius(sp2)
        
      else:
        rad1 = 0.0
        rad2 = 0.0
        
      tempDist = rad1 + rad2
      
      if tempDist > maxDist:
        maxDist = tempDist
  
  return maxDist

def ionicRadius(sym):
    global ionicRadiusDict
     
    try:
        return ionicRadiusDict[sym]
    except:
        sys.exit(__name__+": ERROR: no ionic radius for " + sym)

def initialise():
    global atomicNumberDict, specieDict, atomicMassDict, covalentRadiusDict, ionicRadiusDict
    
    path = os.path.dirname(__file__)
    if len(path):
      fileName = os.path.join(path, 'atoms.in')
    else:
      fileName = 'atoms.in'
    
    if os.path.exists(fileName):
      try:
        f = open(fileName, "r" )
      except:
        sys.exit('error: could not open atoms file: ' + fileName)
    else:
      sys.exit('error: could not find atoms file: ' + fileName)
    
    atomicNumberDict = {}
    specieDict = {}
    atomicMassDict = {}
    covalentRadiusDict = {}
    ionicRadiusDict = {}

    count = 0
    for line in f:
      if count == 0:
        count += 1
      else:
    
        line = line.strip()
        array = line.split(',')
        
        key = array[1].strip()
        
        atomicNumberDict[key] = int(array[0])
        specieDict[count] = key
        atomicMassDict[key] = float(array[2])
        covalentRadiusDict[key] = float(array[3])
        ionicRadiusDict[key] = float(array[4])
  
        count += 1
        
    f.close()
    
if __name__ == '__main__':
    pass
else:
    initialise()



