import os

from ctypes import CDLL, c_double, POINTER, c_int, c_char_p, c_char

from .numpy_utils import CPtrToDouble, CPtrToInt, CPtrToChar
from .numpy_utils import Allocator as alloc


################################################################################

# load lib
_lib = CDLL(os.path.join(os.path.dirname(__file__), "_defects.so"))

################################################################################

# defectVolumeSeparation prototype
_lib.defectVolumeSeparation.restype = c_double
_lib.defectVolumeSeparation.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_int), 
                                        c_int, POINTER(c_int), POINTER(c_double)]

################################################################################

# findDefects prototype
_lib.findDefects.restype = c_int
_lib.findDefects.argtypes = [c_int, c_int, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), 
                                    POINTER(c_int), c_int, POINTER(c_int), c_int, POINTER(c_int), c_int, POINTER(c_int), c_int, POINTER(c_char), 
                                    POINTER(c_int), POINTER(c_double), c_int, POINTER(c_char), POINTER(c_int), POINTER(c_double), POINTER(c_double), 
                                    POINTER(c_int), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int]

# findDefects
def findDefects(includeVacs, includeInts, includeAnts, defectList, NDefectsByType, vacancies, interstitials, antisites, onAntisites, 
                       inclSpec, exclSpecInput, exclSpecRef, NAtoms, specieList, specie, pos, refNAtoms, specieListRef, specieRef, refPos,
                       cellDims, PBC, vacancyRadius, inclusionRadius, minPos, maxPos, verboseLevel, debugDefects):
    """
    findDefects
    
    """
    return _lib.findDefects(includeVacs, includeInts, includeAnts, CPtrToInt(defectList), CPtrToInt(NDefectsByType), CPtrToInt(vacancies), 
                            CPtrToInt(interstitials), CPtrToInt(antisites), CPtrToInt(onAntisites), len(inclSpec), CPtrToInt(inclSpec), 
                            len(exclSpecInput), CPtrToInt(exclSpecInput), len(exclSpecRef), CPtrToInt(exclSpecRef), NAtoms, CPtrToChar(specieList), 
                            CPtrToInt(specie), CPtrToDouble(pos), refNAtoms, CPtrToChar(specieListRef), CPtrToInt(specieRef), CPtrToDouble(refPos),
                            CPtrToDouble(cellDims), CPtrToInt(PBC), vacancyRadius, inclusionRadius, CPtrToDouble(minPos), CPtrToDouble(maxPos), 
                            verboseLevel, debugDefects)