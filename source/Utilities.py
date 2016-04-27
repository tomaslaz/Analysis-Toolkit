"""
Utilities module.

@author Tomas Lazauskas, 2016
@web www.lazauskas.net
@email tomas.lazauskas[a]gmail.com

"""

def distanceSq(pos1x, pos1y, pos1z, pos2x, pos2y, pos2z):
  """
  Returns a squared distance between two positions in the Cartesian system 
  
  """
  
  return ((pos1x - pos2x)**2.0 + (pos1y - pos2y)**2.0 + (pos1z - pos2z)**2.0)
