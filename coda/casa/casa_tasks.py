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


def convolve(model,lineID,bmaj,bmin,pa):

    """

    Convolve a data cube with a beam

    Parameters
    ----------
    model   : path to the ProDiMo model directory (str)
    lineID  : identificator of the emission line fits file produced by ProDiMo (str)
    bmaj    : the beam's major axis [arcsec] (float)
    bmin    : the beam's minor axis [arcsec] (float)
    pa      : the position angle of the beam [deg] (float)

    Output
    ------
    A CASA image with the convolved data cube.

    """

    # Printing info
    print("\n Convolving the model... \n")

    # Define the variable 'cube'
    cube=model+"/LINE_3D_"+lineID+".fits"

    # Convolution
    bmaj=str(bmaj)+"arcsec"
    bmin=str(bmin)+"arcsec"
    pa=str(pa)+"deg"

    # Converting fits image to CASA format.
    ia.fromfits(infile=cube)

    # Convolving image with beam
    imconv=ia.convolve2d(outfile=cube+".conv",
            major=bmaj,minor=bmin,pa=pa,
            overwrite=True)

    # Detaching tool from image
    imconv.done()
    ia.done()

    return None


def subtract_continuum(model,lineID):

    """

    Subtract the continuum emission from a convolved data cube.

    Parameters
    ----------
    model   : path to the ProDiMo model directory (str)
    lineID  : identificator of the emission line fits file produced by ProDiMo (str)

    Output
    ------
    Two CASA images with the continuum flux and a CASA image with the line flux.
    Both are data cubes with the same dimension; the cube for the continuum
    stores the same information in all the channels.

    """

    # Printing info
    print("\n Subtracting the continuum... \n")

    default('imcontsub')

    # Define 'cube' variable
    cube=model+"/LINE_3D_"+lineID+".fits.conv"

    # Get total number of channels from convolved image
    Nchans=imhead(cube,mode='get',hdkey='shape')[2]

    # Define fraction of channels from the total assumed to be line-free
    Fchansfree=3

    # Determine range of line-free channels
    # ----> From Chrsitian (!)
    nrange=max([1,int(Nchans/Fchansfree)])
    Rchans="0~"+str(nrange-1)+","+str(Nchans-nrange)+"~"+str(Nchans-1)

    # Remove file if exists
    # ---> Raise a warning. Maybe Try and except? (!)
    os.system("rm -rf "+cube+".line")
    os.system("rm -rf "+cube+".cont")

    # Subtract continuum using 'imcontsub' task
    imcontsub(imagename=cube,
    linefile=cube+".line",
    contfile=cube+".cont",
    fitorder=0,
    chans=Rchans)

    return None


def create_moments(model,lineID):

    """
    Create moment maps from a convolved and continuum-subtracted
    data cube for the line emission.

    Parameters
    ----------
    model   : path to the ProDiMo model directory (str)
    lineID  : identificator of the emission line fits file produced by ProDiMo (str)

    Output
    ------
    A CASA image with the moment cero of the data set.

    """

    # Printing info
    print("\n Collapsing the cube and creating moments... \n")

    default('immoments')

    # Define variable 'cube'
    cube=model+"/LINE_3D_"+lineID+".fits.conv.line"

    # Remove file if exists
    # ---> Raise a warning. Maybe Try and except? (!)
    os.system("rm -rf "+cube+".mom0")

    # Computing moment cero map of the convolved line cube
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
    model   : path to the ProDiMo model directory [str]
    lineID  : identificator of the emission line fits file produced by ProDiMo [str]
    nu0     : rest frequency of the emission line [Hz] (float)
    bmaj    : the beam's major axis [arcsec] (float)
    bmin    : the beam's minor axis [arcsec] (float)

    Output
    ------
    A CASA image for the integrated line emission in units of K*km/s

    """

    # Printing info
    print("\n Converting units... \n")

    default('immath')

    # Define the variable 'cube'
    cube=model+"/LINE_3D_"+lineID+".fits.conv.line.mom0"

    # Rest frequency in GhZ
    nu0=nu0/1e9

    # Rest frequency squared
    nu0_pow_2=nu0**2

    # Remove file if exists
    # ---> Raise a warning. Maybe Try and except? (!)
    os.system("rm -rf "+cube+".Tb")

    # From intensity to brightness temperature. This method is used here:
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


def convert_to_fits(model,lineID):

    """

    Convert a .Tb CASA image to fits format in order to be analyzed with GoFish

    Parameters
    ----------
    model   : path to the ProDiMo model directory [str]
    lineID  : identificator of the emission line fits file produced by ProDiMo [str]

    """

    # Printing info
    print("\n Converting final product to fits format... \n")

    default('exportfits')

    # Define the 'cube' variable
    cube=model+"/LINE_3D_"+lineID+".fits.conv.line.mom0.Tb"

    # Convert CASA image into a fits image
    exportfits(imagename=cube,
    fitsimage=cube+'.fits',
    overwrite=True)

    return None


################################################################################
# Running the pipeline

model='/Users/bportilla/Documents/project2/ProDiMo_models/run07'
lineID='003'
bmaj=0.21
bmin=0.18
pa=73.1
nu0=3.56734276e11

convolve(model,lineID,bmaj,bmin,pa)
subtract_continuum(model,lineID)
create_moments(model,lineID)
convert_units(model,lineID,nu0,bmaj,bmin)
convert_to_fits(model,lineID)
