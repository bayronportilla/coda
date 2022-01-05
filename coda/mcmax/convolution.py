from coda.mcmax.header import *

def convolve_model(data,fov,npix,beam_x,beam_y,PA_beam):

    ############################################################
    #
    # data: MCMax3D output matrix.
    # fov: the value declared in Image.out. In arcsec.
    # npix: the value declared in Image.out.
    # beam_x: beam width (FWHM) along the major axis in arcsec.
    # beam_y: beam width (FWHM) along the minor axis in arcsec.
    # PA: beam orientation respect to north (direction east)
    # in deg.
    #
    ############################################################

    pxsize=(fov/npix) # Pixel scale: arcsec/px.

    # Convert beam sizes to pixels
    theta_maj=beam_x/pxsize # px
    theta_min=beam_y/pxsize # px

    # Beam orientation
    angle=((90.0+PA_beam)*units.deg).to(units.rad).value

    # Standard deviations
    sigma_maj=theta_maj/2.3548 # px
    sigma_min=theta_min/2.3548 # px

    # Set up kernel
    kernel=Gaussian2DKernel(x_stddev=sigma_maj,y_stddev=sigma_min,theta=angle)

    print("Hi! I'm convolving your model...")
    convolved_data=convolve(data,kernel)

    return convolved_data # mJy/arcsec^2


def convolve_observation(data,pxsize,beam_x,beam_y,PA_beam):

    ############################################################
    #
    # data: matrix to convolve.
    # pxsize: pixel scale in arcsec/px.
    # beam_x: beam width (FWHM) along the major axis in arcsec.
    # beam_y: beam width (FWHM) along the minor axis in arcsec.
    # PA: beam orientation respect to north (direction east)
    # in deg.
    #
    ############################################################

    # Convert beam sizes to pixels
    theta_maj=beam_x/pxsize # px
    theta_min=beam_y/pxsize # px

    # Beam orientation
    angle=((90.0+PA_beam)*units.deg).to(units.rad).value

    # Standard deviations
    sigma_maj=theta_maj/2.3548 # px
    sigma_min=theta_min/2.3548 # px

    # Set up kernel
    kernel=Gaussian2DKernel(x_stddev=sigma_maj,y_stddev=sigma_min,theta=angle)

    print("Hi, I'm convolving your observation...")
    convolved_data=convolve(data,kernel)

    return convolved_data
