"""
Python script for processing disk models using CASA.

This script assumes a full-installation of CASA is available.
"""

import sys

"""
Convolution
"""
def convolve(argv,line):
    """
    Convolve a data cube with a simple beam

    Parameters
    ----------
    model   : path to the model
        The directory where the model resides.
    line    : line identifier
        The line identifier for in the ProDiMo standard
    """

    # Import fits image and convert to CASA format
    cube=model+"LINE_3D_"+line+".fits"
    print(cube)
    #importfits(fitsimage="/Users/bportilla/Documents/project2/ProDiMo_models/run07/LINE_3D_001.fits",imagename="LINE_3D_001.im")
    return None

#convolve()
