from coda.mcmax.header import *

'''

This module is largely dependant on the Gofish package by Richard Teague.

'''

def get_profile(image,dist,pa,inc,aperture,
                visual=None,write=None,fov=None,bunit=None,
                pxscale=None,ebar=None,mstar=None,unit=None):

    """

    Extract the signal using user-defined properties to mask the disk.

    Parameters
    ----------
    image       : path to the data image or data cube (str)
    dist        : distance to the disk [pc] (float)
    pa          : the position angle of the disk [deg] (float)
    inc         : inclination of the disk [deg] (float)
    aperture    : defines orientation, width, and radial stepsize of the
                aperture [deg,deg,arcsec] (list). This quantities are defined
                on the disk frame. For example, orientation=+10 deg means that
                the axis of symetry of the aperture is 10 deg apart from the
                disk's major axis (defined by the position angle argument)

    Output
    ------
    A three-column plain text file containing the distance from the star in
    arcsec, and the average intensity and standard deviation in the units
    specified by the input data cube.

    Example
    -------
    If you want to create an azimuthally averaged profile just use:
    >>> aperture=[0,360,stepsize]

    """

    # Load the data cube
    cube=imagecube(image,FOV=fov,bunit=bunit,pixel_scale=pxscale)

    if cube.data.ndim==2:

        print("\n Working with a two-dimensional fits image \n")

        # Determine limits of the conical aperture
        padir=aperture[0]
        widir=aperture[1]
        dr=aperture[2]
        PAmin=padir-0.5*widir
        PAmax=padir+0.5*widir


        # Extracting radial profile with GoFish
        xm, ym, dym = cube.radial_profile(inc=inc,PA=pa,dist=dist,
                                          x0=0.0,y0=0.0,assume_correlated=False,
                                          PA_min=PAmin,PA_max=PAmax,dr=dr)

        # Plot the aperture?
        if visual:
            vmax=np.percentile(cube.data,99)
            vmin=np.percentile(cube.data,1)
            fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5))
            ax1.imshow(cube.data,origin='lower',extent=cube.extent,vmin=vmin,vmax=vmax)
            cube.plot_mask(ax=ax1,PA_min=PAmin,PA_max=PAmax,
                            inc=inc,PA=pa,mask_frame='disk',r_max=1.5)
            if ebar is True:
                ax2.errorbar(xm,ym,dym)
            if ebar is None:
                ax2.plot(xm,ym)
            ax2.set(xlabel="r (arcsec)",ylabel="Intensity/flux density")
            plt.show()

        # Write to file?
        if write:
            file=open(image+".radial","w")
            for i,j,k in zip(xm,ym,dym):
                file.write("%.15e %.15e %.15e\n"%(i,j,k))
            file.close()

    else:

        print("\n Working with a data cube \n")

        # Determine limits of the conical aperture
        padir=aperture[0]
        widir=aperture[1]
        dr=aperture[2]
        PAmin=padir-0.5*widir
        PAmax=padir+0.5*widir

        # Extracting radial profile with GoFish

        xm, ym, dym = cube.radial_profile(inc=inc,PA=pa,dist=dist,
                                          x0=0.0,y0=0.0,assume_correlated=False,
                                          PA_min=PAmin,PA_max=PAmax,dr=dr,
                                          unit=unit)
        """
        xm, ym, dym = cube.radial_profile(inc=inc,PA=pa,dist=dist,
                                          unit=unit)
        """

        file=open("13CO.radial.cube.WithAperture","w")
        for i,j,k in zip(xm,ym,dym):
            file.write("%.15e %.15e %.15e\n"%(i,j,k))
        file.close()

        #plt.plot(xm,ym,'.')
        #plt.show()






    return None
