import os

# Boltzman constant
kB = 8.6173324 * 0.00001 # eV/K
planckConst = 4.13566766225 * (10**(-15))
reducedPlanckConst = 6.582119514 * (10**(-16))
avogadroConst = 6.022140857* (10**(23))
lightSpeedConts = 29979245800 # cm/s

verbosity = 0

# Environment
_sourceDir = os.path.dirname(os.path.abspath(__file__))

# Symmetry
_cluspyErrFile = "cluspy.err"
_symmetrizerPath = os.path.join(os.path.dirname(__file__), "/Users/Tomas/git/Analysis-Toolkit/thirdparty/symmetrizer/symmetrizer.jar")

# Scott key
_scottKeySeparator = "-"
