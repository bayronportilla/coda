"""

Python script for processing ProDiMo models using CASA.
--------------------------------------------------------------------------------
This script assumes a full-installation of CASA and it has only been tested with
CASA version 6.4.

How to use this module?
--------------------------------------------------------------------------------
There are several methods that can be used independently. For a flawless
procedure however, it is adviced to run the whole set of routines in the
following order:

    - convolve()
    - subtract_continuum()
    - create_moments()
    - convert_units()
    - convert_to_fits()

"""

import sys

def convolve():
    """
    Convolve a data cube with a simple beam

    Parameters
    ----------
    model   : path to the model - argv[1]
        The directory where the model resides.
    line    : line identifier - argv[2]
        The line identifier for in the ProDiMo standard
    """

    # Define variables
    #model=sys.argv[1]
    #line=sys.argv[2]
    model='/Users/bportilla/Documents/project2/ProDiMo_models/run07'
    line='001'
    cube=model+"/LINE_3D_"+line+".fits"

    # Import fits image and convert to CASA format
    #importfits(fitsimage=cube,imagename="LINE_3D_"+line+".im",overwrite=True)

    # Convolution
    #beam=["5.833333333333333e-05deg","4.9999999999999996e-05deg","73.1deg"]
    beam=["1.0278e-04deg","9.1666e-05deg","67.1deg"] # ---> Can this be in arcsec? (!)
    #ia=iatool()
    #ia.open(infile=cube) # Attaching image analysis tool to the specified cube
    #ia.close()
    ia.fromfits(infile=cube) # Converting fits image to CASA format.
    #imconv=ia.convolve2d(outfile="LINE_3D_"+line+".fits.conv",
    imconv=ia.convolve2d(outfile=cube+".conv",
    major=beam[0],minor=beam[1],pa=beam[2],
    overwrite=True)

    imconv.done()
    ia.done()
    return None

def subtract_continuum():
    """
    Subtract the continuum from a convolved data cube

    Parameters
    ----------

    """

    default('imcontsub')

    model='/Users/bportilla/Documents/project2/ProDiMo_models/run07'
    line='001'
    cube=model+"/LINE_3D_"+line+".fits.conv"

    # Get number of channels from convolved image
    Nchans=imhead(cube,mode='get',hdkey='shape')[2]

    # Define fraction of channels assumed to be line-free
    Fchansfree=3

    # Determine range of line-free channels
    # ----> From Chrsitian (!)
    nrange=max([1,int(Nchans/Fchansfree)])
    Rchans="0~"+str(nrange-1)+","+str(Nchans-nrange)+"~"+str(Nchans-1)

    # ----> Removing file if exists (!)

    # Subtracting continuum using imcontsub task
    imcontsub(imagename=cube,
    linefile=cube+".line",
    contfile=cube+".cont",
    fitorder=0,
    chans=Rchans)

    return None

def create_moments():
    """
    Create moment maps from a convolved and continuum subtracted
    data cube.

    Parameters
    ----------
    """

    default('immoments')

    model='/Users/bportilla/Documents/project2/ProDiMo_models/run07'
    line='001'
    cube=model+"/LINE_3D_"+line+".fits.conv.line"

    # Computing first moment of the convolved line cube
    immoments(imagename=cube,
    moments=[0],
    outfile=cube+".mom0")

    return None


def convert_units(model,lineID,nu0,bmaj,bmin):

    """

    Convert line cube units from Jy beam^-1 Km s^-1 to
    K Km s^-1.

    Parameters
    ----------
    model:  path to the ProDiMo model directory [str]
    lineID: identificator of the emission line fits file produced by ProDiMo [str]
    nu0:    rest frequency of the emission line [Hz]
    bmaj:   the beam's major axis [arcsec]
    bmin:   the beam's minor axis [arcsec]

    Output
    ------
    A CASA image for the integrated line emission in K*km/s

    """

    default('immath')

    #model='/Users/bportilla/Documents/project2/ProDiMo_models/run07'
    #line='001'
    cube=model+"/LINE_3D_"+lineID+".fits.conv.line.mom0"

    # Rest frequency in GhZ
    nu0=nu0/1e9

    # Rest frequency squared
    nu0_pow_2=nu0**2

    # From intensity to brightness temperature. The method is also used here:
    # https://casaguides.nrao.edu/index.php/VLA_CASA_Imaging-CASA5.0.0#Image_Conversion
    immath(imagename=cube,
    mode='evalexpr',
    expr=('1.222e6*IM0/%s/(%s*%s)'%(nu0_pow_2,bmaj,bmin)),
    outfile=cube+'.Tb')

    # Change unit in header
    imhead(imagename=cube+'.Tb',
    mode='put',
    hdkey='bunit',
    hdvalue='K.km/s')

    return None


def convert_to_fits():

    """

    Convert a .Tb CASA image to fits format to
    be analyzed with GoFish

    Parameters
    ----------

    """

    default('exportfits')

    model='/Users/bportilla/Documents/project2/ProDiMo_models/run07'
    line='001'
    cube=model+"/LINE_3D_"+line+".fits.conv.line.mom0.Tb"

    # Convert CASA image into a fits image
    exportfits(imagename=cube,
    fitsimage=cube+'.fits',
    overwrite=True)

    return None

#convolve()

#subtract_continuum()

#create_moments()


convert_units('/Users/bportilla/Documents/project2/ProDiMo_models/run07',
'001',
2.30537939e11,
0.37,
0.33)

#convert_to_fits()
