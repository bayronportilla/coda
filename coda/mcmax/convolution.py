from coda.mcmax.header import *

def convolve_image(data):

    # Start convolution
    pxsize=4.0 # mas/px
    beam_x=74.0 # mas
    beam_y=57.0 # mas

    beam_x=beam_x/pxsize # px
    beam_y=beam_y/pxsize # px

    PA=63.0
    angle=0.0#((PA+90.0)*units.deg).to(units.rad).value

    kernel=Gaussian2DKernel(x_stddev=beam_x,y_stddev=beam_y,theta=angle)
    convolved_data=convolve(data,kernel)

    return convolved_data

def convolve_image_jband(data):

    # Start convolution
    pxsize=12.26 # mas/px
    beam_x=24.5 # (mas) 50% of the image resolution in px
    beam_y=beam_x # mas

    beam_x=beam_x/pxsize
    beam_y=beam_y/pxsize

    #beam_x=beam_x/pxsize # px
    #beam_y=beam_y/pxsize # px

    PA=63.0
    angle=1.94

    kernel=Gaussian2DKernel(x_stddev=beam_x,y_stddev=beam_y,theta=angle)
    convolved_data=convolve(data,kernel)

    return convolved_data

