"""

Python script for processing ProDiMo models using CASA.
--------------------------------------------------------------------------------
This script assumes a full-installation of CASA and it has only been tested with
CASA version 6.4.

How to use this module?
--------------------------------------------------------------------------------
There are several methods that can be used independently. For a flawless
run however, it is adviced to run the whole set of routines in the following
order:

    - rotate()
    - convolve()
    - subtract_continuum()
    - create_moments()
    - convert_units()
    - convert_to_fits()

"""

import sys

def rotate(file,angle):

    """

    Rotates a cube by an angle measured east-from-north. Note that this
    operation first rotates the two direction axes of an cube and then, the
    input image is regridded to the rotated coordinate system. Therefore, an
    internal interpolation of the data takes place. If the image brightness
    units are Jy/pixel, the output is scaled to conserve flux.

    Parameters
    ----------
    model   : path to the ProDiMo model directory (str)
    lineID  : identificator of the emission line fits file produced by ProDiMo (str)
    angle   : the amount the two axes of the image will be rotated by. This angle
            is reckoned east-of-north [deg] (float).

    Output
    ------
    A CASA image with the rotated data cube.

    """

    # Printing info
    print("\n Rotating the file... \n")

    # Define the variables 'cube' and 'angle'
    cube = file
    angle = str(angle)+"deg"

    # Converting fits image to CASA format
    ia.fromfits(infile=cube)

    # Rotating image
    imrot=ia.rotate(outfile=cube+".rot",
            pa=angle,
            overwrite=True)

    # Detaching tool from image
    imrot.done()
    ia.done()

    return None


def convolve(file,bmaj,bmin,pa):

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
    cube=file

    # Convolution
    bmaj=str(bmaj)+"arcsec"
    bmin=str(bmin)+"arcsec"
    pa=str(pa)+"deg"

    # Attach the Image Analysis tool to the specified image.
    ia.open(cube)

    # Convolving image with beam
    imconv=ia.convolve2d(outfile=cube+".conv",
            major=bmaj,minor=bmin,pa=pa,
            overwrite=True)

    # Detaching tool from image
    imconv.done()
    ia.done()

    return None


def subtract_continuum(file):

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
    cube=file

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


def create_moments(file):

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
    cube = file

    # Remove file if exists
    # ---> Raise a warning. Maybe Try and except? (!)
    os.system("rm -rf "+cube+".mom0")

    # Computing moment cero map of the convolved line cube
    immoments(imagename=cube,
    moments=[0],
    outfile=cube+".mom0")

    return None


def convert_units(file,nu0,bmaj,bmin):

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
    cube=file

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


def convert_to_fits(file):

    """

    Convert a CASA image with extension .prefix to fits format.

    Parameters
    ----------
    model   : path to the ProDiMo model directory [str]
    lineID  : identificator of the emission line fits file produced by ProDiMo [str]
    prefix  : The prefix of the CASA file to be converted into fits format.
            The prefix are those characters placed to the right of the rightmost
            period in the measurement set's name.

    """

    # Printing info
    print("\n Converting %s file into fits format... \n"%(file))

    default('exportfits')

    # Define the 'cube' variable
    cube=file

    # Convert CASA image into a fits image
    exportfits(imagename=cube,
    fitsimage=cube+'.fits',
    overwrite=True)

    return None


################################################################################
# Running the pipeline
#file='~/Documents/Models/prodimo/PDS-70/Portilla-Revelo-2023-reduced.v2/image-cont-855.fits'
file='/Users/bportilla/Documents/Models/prodimo/PDS-70/Portilla-Revelo-2023-reduced.v2/LINE_3D_001.fits'
#file='/Users/bportilla/Documents/Models/prodimo/PDS-70/Portilla-Revelo-2023-reduced.v2/image-cont-855.fits'
angle=70.4 # Angle for cube rotation ---> Explain this in detail (!)

"""
bmaj    = 0.067
bmin    = 0.050
pa      = 61.5
"""

bmaj    = 0.13
bmin    = 0.10
pa      = -83.75
nu0     = 2.30537939e11

rotate(file,angle)
convolve(file+'.rot',bmaj,bmin,pa)
convert_to_fits(file+'.rot.conv')
subtract_continuum(file+'.rot.conv')
convert_to_fits(file+'.rot.conv.line')
create_moments(file+'.rot.conv.line')
convert_to_fits(file+'.rot.conv.line.mom0')
convert_units(file+'.rot.conv.line.mom0',nu0,bmaj,bmin)
convert_to_fits(file+'.rot.conv.line.mom0.Tb')

