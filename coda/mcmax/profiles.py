from coda.mcmax.header import *

'''

This module depends on the GoFish package by Richard Teague
https://joss.theoj.org/papers/10.21105/joss.01632

'''

def get_profile(image,dist,pa,inc,aperture,
                visual=None,write=None,fov=None,bunit=None,
                pxscale=None,ebar=None,mstar=None,unit=None,esurf=None):

    """

    Extract the signal using user-defined properties to mask the disk.

    Parameters
    ----------
    image       : path to the data image or data cube (str)
    dist        : distance to the disk [pc] (float)
    pa          : the position angle of the disk [deg] (float)
    inc         : inclination of the disk [deg] (float)
    aperture    : defines orientation, angular width, and radial stepsize of the
                aperture [deg,deg,arcsec] (list). This quantities are defined
                on the disk frame. For example, orientation=+10 deg means that
                the axis of symetry of the aperture is 10 deg shifted from the
                disk's major axis (defined by the position angle argument)
    esurf       : emitting surface file. Single column file containing the parameters
                of the surface of emission as required by GoFish. Rows are: z0,
                psi, r_taper, q_taper, and r_cavity. Distances in arcsec, exponents
                are dimensionaless.

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

    # Loading parameters for surface of emission
    if esurf is not None:
        params=pd.read_csv(esurf,comment='#',sep=' ',header=None)
        z0=params.loc[0,0]          # arcsec
        psi=params.loc[1,0]
        r_taper=params.loc[2,0]     # arcsec
        q_taper=params.loc[3,0]
        r_cavity=params.loc[4,0]    # arcsec
    else:
        z0=None
        psi=None
        r_taper=None
        q_taper=None
        r_cavity=None

    # Parametrized surface of emission
    def Z(r):                       # r is in arcsec
        if esurf is not None:
            num = z0 * (r-r_cavity)**psi * np.exp(-((r-r_cavity)/r_taper)**q_taper)
            value = num/r
        else:
            value = np.zeros(len(r))
        return value


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
                                          PA_min=PAmin,PA_max=PAmax,dr=dr,
                                          z0=z0,psi=psi,r_cavity=r_cavity,
                                          r_taper=r_taper,q_taper=q_taper)

        # Plot the aperture?
        if visual:

            vmax=np.percentile(cube.data,99)
            vmin=np.percentile(cube.data,1)

            fig,((ax1,ax2,ax5),(ax3,ax4,ax6))=plt.subplots(2,3,figsize=(11,7))

            # Image with mask
            ax1.imshow(cube.data,origin='lower',extent=cube.extent,vmin=vmin,vmax=vmax)
            cube.plot_mask(ax=ax1,PA_min=PAmin,PA_max=PAmax,
                            inc=inc,PA=pa,mask_frame='disk',r_max=1.5)

            # Emitting surface
            cube.plot_surface(ax=ax2,inc=inc,PA=pa,z0=z0,psi=psi,r_cavity=r_cavity,
                                r_taper=r_taper,q_taper=q_taper)

            # Image + surface
            ax3.imshow(cube.data,origin='lower',extent=cube.extent,vmin=vmin,vmax=vmax)
            cube.plot_surface(ax=ax3,inc=inc,PA=pa,z0=z0,psi=psi,r_cavity=r_cavity,
                                r_taper=r_taper,q_taper=q_taper)

            # Emitting surface one dimensional profile
            r_array = np.linspace(xm[0],xm[-1],100)
            Z_array = Z(r_array)
            ax5.plot(r_array*dist,Z_array,".")
            ax5.set_xlabel("distance (au)")
            ax5.set_ylabel("z/r")



            # Radial profile
            if ebar is True:
                ax4.errorbar(xm,ym,dym)
            if ebar is None:
                ax4.plot(xm,ym)
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
        padir   = aperture[0]
        widir   = aperture[1]
        dr      = aperture[2]
        PAmin   = padir-0.5*widir
        PAmax   = padir+0.5*widir

        # Extracting radial profile with GoFish
        xm, ym, dym = cube.radial_profile(inc=inc,PA=pa,dist=dist,
                                          x0=0.0,y0=0.0,assume_correlated=False,
                                          PA_min=PAmin,PA_max=PAmax,dr=dr,
                                          z0=z0,psi=psi,r_cavity=r_cavity,
                                          r_taper=r_taper,q_taper=q_taper,
                                          unit=unit)

        # Write to file?
        if write:
            file=open(image+"."+unit+".radial","w")
            for i,j,k in zip(xm,ym,dym):
                file.write("%.15e %.15e %.15e\n"%(i,j,k))
            file.close()



        # Plotting teardrop plot
        #cube.plot_teardrop(inc=inc, PA=pa,mstar=0.76, dist=dist)
        #cube.plot_teardrop(inc=51.7, PA=160.4,mstar=1.0, dist=dist)

        #plt.show()

        #plt.plot(xm,ym,'.')
        #plt.show()

    return None
