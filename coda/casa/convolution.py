#from coda.mcmax.header import *

"""
Python script for processing disk models using CASA.
"""

def convolve(cube,beam):
    """
    Convolves a data cube with a simple beam.

    Parameters
    ----------
    cube    : fits data cube
        The model data cube.
    beam    : list
        A list containing the major-axis, minor-axis, and the
        orientation of the beam. Units are arcsec, arcsec and degree,
        respectively.

    Example
    -------
    >>> convolve([0.15,0.10,45])
    """
    #from astropy import units as u
    import astropy

    bmaj=beam[0]#*u.arcsec
    bmin=beam[1]#*u.arcsec
    bang=beam[2]#*u.deg

    ia.fromfits(infile=cube)
    #ia.imhead()

    return None

convolve("/Users/bportilla/Documents/project2/ProDiMo_models/run07/LINE_3D_001.fits",[0.1,0.1,45])
