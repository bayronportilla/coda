from coda.mcmax.header import *

@dataclass
class File:
    x:float
    y:float
    name:str

    def calculate_mass(self,rlim=None):

        """

        Integrate the density profile to find the mass.

        Parameters
        ----------
        rlim    : (rmin,rmax) to only integrate between rmin <= r <= rmax.
                If rlim is None, calculate the total mass. [au,au] (tuple).

        """

        x=(self.x*u.au).to(u.cm).value
        y=self.y

        if rlim is not None:
            # Convert limits to cm
            rmin,rmax=rlim
            rmin=(rmin*u.au).to(u.cm).value
            rmax=(rmax*u.au).to(u.cm).value

            # Start integration
            xint,yint=[],[]
            for i,j in zip(x,y):
                if i>=rmin and i<=rmax:
                    xint.append(i)
                    yint.append(j)
            xint,yint=np.array(xint),np.array(yint)
            dust_mass=(2*np.pi*simps(xint*yint,xint)*u.g).to(u.Msun)
        else:
            dust_mass=(2*np.pi*simps(x*y,x)*u.g).to(u.Msun)

        return dust_mass

    def plot(self):
        x,y=self.x,self.y
        fig,ax=plt.subplots(1,1)
        ax.plot(x,y,'.')
        ax.set(xscale="log",yscale="log")
        plt.show()

    def rescale_mass(self,k,rlim=None):

        """

        Rescale surface density profile by a factor of k.

        Parameters
        ----------

        k       : rescaling factor. [dimensionless] (float).
        rlim    : (rmin,rmax) to only rescale the points with rmin <= r <= rmax.
                If rlim is None, rescale the entire profile. [au,au] (tuple).

        """

        x=(self.x*u.au).to(u.cm).value
        y=self.y

        Mold=(2*np.pi*simps(x*y,x)*u.g).to(u.Msun)

        # Rescale the whole profile
        if rlim is None:
            ynew=k*y

        # Rescale between rmin and rmax
        if rlim is not None:
            rmin,rmax=rlim
            ynew=[]
            for i in range(len(self.x)):
                if self.x[i]>=rmin and self.x[i]<=rmax:
                    ynew+=[k*y[i]]
                else:
                    ynew+=[1*y[i]]

        Mnew=(2*np.pi*simps(x*ynew,x)*u.g).to(u.Msun)

        ''' Print info '''
        print("\nInitial mass:",Mold)
        print("Rescaled mass:",Mnew)

        ''' Create rescaled mass file '''
        f=open(self.name+".rescaled","w")
        for i,j in zip(self.x,ynew):
            f.write("%.15e %.15e\n"%(i,j))
        f.close()

        return None

    def plot_MCSpec(self,unit=None):
        if not unit:
            x,y=self.x,self.y
            fig,ax=plt.subplots(1,1)
            ax.plot(x,y)
            ax.set(xscale="log",yscale="log")
            plt.show()


def get_size_distribution(model):

    '''

    Returns an array with the particle's bin sizes of a particle size
    distribution.

    Caution!
    --------
    * This routine only works when a single type of dust grain is present,
    i.e. when only one computepartXX keyword exists in the input.dat file

    '''

    ''' Read file names '''
    filenames=[]
    folder=model+"/Particles/"
    for filename in os.listdir(folder):
        if fnmatch.fnmatch(filename,('*.fits.gz')):
            filenames.append(filename)

    psize=np.zeros((len(filenames),3))
    i=0
    for file in filenames:
        hdulist=fits.open(model+'/Particles/'+file)
        hdu=hdulist[0]
        hdr=hdu.header
        a=hdr['A1']
        amin=hdr["R_MIN"]
        amax=hdr["R_MAX"]
        psize[i]=[a,amin,amax]
        i+=1
    return psize


def upload_file(file):

    '''
    Create a File object from a two-column file
    '''

    data=np.loadtxt(file)
    x=np.reshape(data[:,0:1],data.shape[0])
    y=np.reshape(data[:,1:2],data.shape[0])
    name=file

    return File(x,y,name)


def density_profile(Rin,Rout,Mdust,
                    Rtap=None,epsilon=None,gamma=None):
    k=10
    N_points=2**k+1
    Rin,Rout,Mdust=Rin*u.au,Rout*u.au,Mdust*u.Msun
    rarray=np.linspace(Rin,Rout,N_points)
    Rtap=Rtap*u.au if Rtap is not None else 100*u.au
    if epsilon is None:epsilon=1.0
    if gamma is None:gamma=1.0
    def integrand(x):
        return x*(Rtap.value/x)**epsilon*np.exp(-(x/Rtap.value)**(2-gamma))
    Sigma0=(Mdust/(2*np.pi*quad(integrand,Rin.value,Rout.value)[0]*u.au**2)).to(u.g/u.cm**2)
    Sigma=[(Sigma0.value*(Rtap/i)**epsilon*np.exp(-(i/Rtap)**(2-gamma))).value for i in rarray]
    file=open("density.dat","w")
    for i,j in zip(rarray,Sigma):
        file.write("%.15e %.15e\n"%(i.value,j))
    file.close()
    return None


def convert_comp(fc,porosity,qtype) :

    '''

    Convert grain composition from the MCMax3D standard to the ProDiMo standard
    and vice-versa.

    Parameters
    ----------
    fc: fraction of carbon.

    porosity: the fraction of vaccum the grain is made of. If the amount of
        vacuum is say 10%, then you must enter porosity=0.1.

    qtype: The type of the input quantities. It can be either 'mcmax' or
        'prodimo'. If qtype='mcmax', then the routine will return the correspondent
        ProDiMo quantities. If qtype='prodimo', then the routine will return the
        correspondent MCMax3D quantities.

    '''

    #Bulk constants. Do not touch them unless the hard-coded quantities in
    #MCMax3D have changed.
    rho_c=1.8   # g cm^-2
    rho_si=3.01 # g cm^-2

    #From MCMax3D to ProDiMo
    if qtype=="mcmax":
        fcM=fc
        ''' Derived quantities '''
        fsiM=1-fcM
        fcV=(1-rho_c/rho_si*(1-fcM**(-1)))**(-1)
        fsiV=1-fcV

        ''' Converting to ProDiMo quantities '''
        FsiV=(1-porosity)*(fcV/fsiV+1)**(-1)
        FcV=(fcV/fsiV)*FsiV

        print("\nThe correspondent ProDiMo quantities are:")
        print("FsiV=%.15f"%(FsiV))
        print("FcV=%.15f"%(FcV))
        print("FsiV+FcV+porosity=",FsiV+FcV+porosity)

    #From ProDiMo to MCMax3D
    elif qtype=="prodimo":
        FcV=fc
        ''' Derived quantities '''
        FsiV=1-FcV-porosity

        ''' Converting to MCMax3D quantities '''
        fsiV=(1+FcV/FsiV)**(-1)
        fcV=1-fsiV
        fcM=(1-(rho_si/rho_c)*(1-fcV**-1))**(-1)
        fsiM=1-fcM

        print("\nThe correspondent MCMax3D quantities are:")
        print("fcM=%.15f"%(fcM))
        print("fsiM=%.15f"%(fsiM))
        print("fsiM+fcM=%.15f"%(fsiM+fcM))

    return None


def convert_density_file(model,g2d=None,visual=None,find_dust_mass=None):

    '''

    Convert a standard MCMax3D grid density output into a 1D input sdfile for
    ProDiMo.

                        r
                     ------>
                    |       ^
    M_mgc =  theta  |       | z
                    v       |
                     -------

                        r
                     ------>
                    |       ^
    M_pgc =      z  |       | theta
                    v       |
                     -------

    Caution!
    --------
    *   The length of the r_array will determine the resolution of the ProDiMo
        model. Choose it wisely!

    *   This routine only works when the MCMax3D model has only one type of
        particle i.e. when there is only one keyword of the type 'computepartXX'
        in the input.dat file.

    *   When using this routine, I assume you have already ran a 'stop-after-init'
        ProDiMo simulation. This way you have a grid where MCMax3D quntities will
        be interpolated to.

    *   The density profile inside the folder must be named as:
        'surface_density_PDS70_70_cropped.dat'

    *   This routine works only with 2D models.

    *   MCMax3D model must have exclussively two zones; no one, no three...

    *   Rin and Rout must be equal to 0.04 au and 130 au respectively.

    *   The number of radial points of the ProDiMo grid must be equal to that in
        the MCMax3D model.

    *   For production results, it is advisable to run the MCMax3D model with a
        refined radial grid, Nr>=60 per zone as long as this is computationally
        feasible.

    Parameters
    ----------
    model   : path to the MCMax3D model directory. (str).
    g2d     : The gas-to-dust ratio array. If g2d is None assume to be 100
            everywhere. If len(g2d)=1, use a constant g2d ratio of g2d[0]. If
            len(g2d)!=1 then g2d *must* be defined at every radial point in the
            grid.
    visual  : If true, shows the computed density profile for a selected dust
            sizes and also the reconstructed profile compared to the original one.
            (Bool).

    Example
    -------
    1. Run a stop_after_init ProDiMo model to create a cylindrical grid if needed.
    Remember to run with interface_1D=.false. Remember also to set the Nxx variable
    to match the total number of radial points of the MCMax3D model.
    >>> cd <path-to-ProDiMo-model>
    >>> nano Parameter.in
    >>> prodimo

    2. In the directory of the ProDiMo model, open an ipython session and import
    the modules
    >>> ipython
    >>> from coda.mcmax import input

    3. Call the module passing the appropriate parameters
    >>> input.convert_density_file("<path-to-MCMax3D-model>",visual=True,find_dust_mass=False)


    '''

    @dataclass
    class CellMcmax:
        r:float             # cm
        theta:float
        phi:float
        r_max:float
        r_min:float
        theta_max:float
        theta_min:float
        phi_max:float
        phi_min:float
        comp:float

        def dV(self):
            # Units: cm^3

            dr=self.r_max-self.r_min
            dtheta=self.theta_max-self.theta_min
            dphi=self.phi_max-self.phi_min

            value=self.r**2*np.sin(self.theta)*dr*dtheta*dphi

            return value

        def dz(self):
            # Units: cm
            value=self.r*self.dtheta
            return value

        def zsphere(self):

            '''

            Finds the projection of the radial spherical coordinate onto the
            vertical axes.

            Output
            ------
            Projected radius in cm.

            '''

            if self.theta==(90*u.deg).to(u.rad).value:
                value=0.0
            else:
                value=self.r*np.cos(self.theta)

            return value

        def rsphere(self):

            '''

            Finds the projection of the radial spherical coordinate onto the
            midplane

            Output
            ------
            Projected radius in cm.

            '''

            #if self.r==(0.04*u.au).to(u.cm).value
            value=self.r*np.sin(self.theta)

            return value

    @dataclass
    class CellProdimo:
        r:float
        z_min:float
        z_max:float
        comp:float

        def dz(self):

            '''

            Vertical size of the cylindric cell

            Output
            ------
            Vertical size in cm.

            '''

            value=((self.z_max-self.z_min)*u.au).to(u.cm).value

            return value


    # Computing psize
    print("\n Computing psize... \n")
    psize=get_size_distribution(model)
    ai_array=psize[:,0] # microns

    def func_fsize(zoneID):
        # Import data from zones
        hdu_1=fits.open(model+"/output/Zone000%d.fits.gz"%(zoneID))

        # Mass density matrix (g cm^-3)
        # Remember that C[3]=pmid(rad), C[4]=tmid(rad), C[5]=rmid(au)
        C=hdu_1[6].data

        # Coordinates at cell center
        rmidAU=hdu_1[0].data[0,0,0,:]
        rmid=(rmidAU*u.au).to(u.cm).value # cm
        tmid=hdu_1[0].data[1,0,:,0]       # rad
        pmid=hdu_1[0].data[2,:,0,0]       # rad

        # Coordinates at cell border
        rlim=hdu_1[1].data   # cm
        tlim=hdu_1[2].data   # rad
        plim=hdu_1[3].data   # rad

        # Correcting for radius
        if zoneID==1:
            rmidAU_min=0.04
            rmid_min=(rmidAU_min*u.au).to(u.cm).value
            rmidAU[0]=rmidAU_min
            rmid[0]=rmid_min

        if zoneID==2:
            rmidAU_max=130.0
            rmid_max=(rmidAU_max*u.au).to(u.cm).value
            rmidAU[-1]=rmidAU_max
            rmid[-1]=rmid_max

        # Correcting for theta
        theta_midplane=(90.0*u.deg).to(u.rad).value
        tmid[int(len(tmid)*0.5-1)]=theta_midplane

        # Declare matrix required for interpolation
        # Dim: len(theta)/2 x len(radius)
        M_mgc=np.empty((int(len(tmid)*0.5),len(rmid)),dtype='object')

        if find_dust_mass is True:
            # The big loop
            for i in range(C.shape[5]):             # Along r
                for j in range(C.shape[4]):         # Along theta
                    for k in range(C.shape[3]):     # Along phi
                        comp=C[0,:,0,k,j,i] # g cm^-3
                        GP=CellMcmax(rmid[i],tmid[j],pmid[k],
                                    rlim[i+1],rlim[i],
                                    tlim[j+1],tlim[j],
                                    plim[k+1],plim[k],
                                    comp)
                        if k==0.0 and j<len(tmid)*0.5:
                            M_mgc[j,i]=GP

        if find_dust_mass is False:
            # The not-too-big loop
            for i in range(C.shape[5]):             # Along r
                for j in range(C.shape[4]):         # Along theta
                    if j<int(len(tmid)*0.5):
                        comp=C[0,:,0,0,j,i] # g cm^-3
                        GP=CellMcmax(rmid[i],tmid[j],pmid[0],
                                    rlim[i+1],rlim[i],
                                    tlim[j+1],tlim[j],
                                    plim[1],plim[0],
                                    comp)
                        M_mgc[j,i]=GP
                    else:
                        continue

        return rmidAU,M_mgc

    # Concatenate fsizes and rmidAU
    rmidAU_1,M_mgc_1=func_fsize(1)[0],func_fsize(1)[1]
    rmidAU_2,M_mgc_2=func_fsize(2)[0],func_fsize(2)[1]
    M_mgc_full=np.concatenate((M_mgc_1,M_mgc_2),axis=1)
    r_array=np.concatenate((rmidAU_1,rmidAU_2))

    # Read ProDiMo grid. This steps requires prodimopy!
    model_prodimo=read_prodimo()
    r_array_p=np.reshape(model_prodimo.x[:,0:1],model_prodimo.x.shape[0])   # au
    z_matrix=model_prodimo.z                             # au, dim:len(r)*len(z)

    # Declaring and filling in M_pgc_full matrix for ProDiMo. Note that
    # composition is initialized to an array of ones.
    M_pgc_full=np.empty((z_matrix.shape[1],z_matrix.shape[0]),dtype='object')
    for i in range(M_pgc_full.shape[0]-1):       # Over z
        for j in range(M_pgc_full.shape[1]):     # Over radius
            GP_prodimo=CellProdimo(r_array_p[j],
                                    z_matrix[j,i],z_matrix[j,i+1],
                                    np.ones(psize.shape[0]))
            M_pgc_full[i,j]=GP_prodimo

    # Validate this step ---> (!)
    for i in range(M_pgc_full.shape[0]-1):
        for j in range(M_pgc_full.shape[1]):
            M_pgc_full[i,j].r=(M_mgc_full[-1,j].r*u.cm).to(u.au).value

    # Printing grid's limits
    # Radial limits MCMax
    print("\n Following are the coordinates for the cell's center (for MCMax)!!!")
    print("----------------------------------------------------------")
    rmin_grid_mcmax=(M_mgc_full[-1,0].r*u.cm).to(u.au).value
    rmax_grid_mcmax=(M_mgc_full[-1,-1].r*u.cm).to(u.au).value

    # Vertical limits MCMax
    zmin_grid_mcmax=(M_mgc_full[-1,0].zsphere()*u.cm).to(u.au).value
    zmax_grid_mcmax=(M_mgc_full[0,-1].zsphere()*u.cm).to(u.au).value

    # Radial limits ProDiMo
    rmin_grid_prodimo=M_pgc_full[0,0].r
    rmax_grid_prodimo=M_pgc_full[0,-1].r

    # Vertical limits ProDiMo
    zmin_grid_prodimo=M_pgc_full[0,0].z_min
    zmax_grid_prodimo=M_pgc_full[-2,-1].z_max

    if rmin_grid_prodimo<rmin_grid_mcmax:
        print("\nWarning!")
        print("r_min_prodimo < r_min_mcmax %.15f and %.15f"%(rmin_grid_prodimo,rmin_grid_mcmax))
        #print("Correcting...")
        #M_mgc_full[-1,0].r=(rmin_grid_prodimo*au).to(u.cm).value

    if rmax_grid_prodimo>rmax_grid_mcmax:
        print("\nWarning!")
        print("r_max_prodimo > r_max_mcmax %.15f and %.15f"%(rmax_grid_prodimo,rmax_grid_mcmax))
        #print("Correcting...")
        #M_mgc_full[-1,-1].r=(rmax_grd_prodimo*au).to(u.cm).value

    if zmin_grid_prodimo<zmin_grid_mcmax:
        print("\nWarning!")
        print("z_min_prodimo < z_min_mcmax %.15f and %.15f"%(zmin_grid_prodimo,zmin_grid_mcmax))
        #print("Correcting...")
        #M_mgc_full[-1,0].

    if zmax_grid_prodimo>zmax_grid_mcmax:
        print("\nWarning!")
        print("z_max_prodimo > z_max_mcmax %.15f and %.15f"%(zmax_grid_prodimo,zmax_grid_mcmax))

    # Do 2D interpolation. Sampling MCMax3D info into 1D arrays.
    #for k in range(psize.shape[0]):
    for k in range(psize.shape[0]):
        rsph_array=[]
        zsph_array=[]
        csph_array=[]
        for i in range(M_mgc_full.shape[0]):                           # Over theta
            for j in range(M_mgc_full.shape[1]):                       # Over radius
                rsph_array.append((M_mgc_full[i,j].rsphere()*u.cm).to(u.au).value)
                zsph_array.append((M_mgc_full[i,j].zsphere()*u.cm).to(u.au).value)
                csph_array.append(np.log10(M_mgc_full[i,j].comp[k]))

        # Uncomment the preferred interpolator. Rbf seems to work nicely.
        #f=interp2d(rsph_array,zsph_array,csph_array,kind='linear',bounds_error=False)
        f=Rbf(rsph_array,zsph_array,csph_array,function='linear')

        # Interpolate logarithm of densities into ProDiMo grid
        for i in range(M_pgc_full.shape[0]-1):                 # Over z
            for j in range(M_pgc_full.shape[1]):               # Over radius
                rprodimo=M_pgc_full[i,j].r                     # au
                zprodimo=M_pgc_full[i,j].z_min                 # au
                M_pgc_full[i,j].comp[k]=f(rprodimo,zprodimo)

    #Interpolate density profile
    fobj=upload_file(model+"/surface_density_PDS70_70_cropped.dat")
    x_array=fobj.x
    y_array=fobj.y
    cs=CubicSpline(x_array,y_array)
    S_array=np.array([cs(i) for i in r_array])

    # Declaring fsize matrix
    fsize=np.zeros((psize.shape[0],M_pgc_full.shape[1]))
    for k in range(psize.shape[0]):               # Over composition
        for j in range(M_pgc_full.shape[1]):      # Over radius
            S_val=0.0
            for i in range(M_pgc_full.shape[0]-1):  # Over height
                rho_val=10**float(M_pgc_full[i,j].comp[k])
                S_val+=rho_val*M_pgc_full[i,j].dz()
            # Multiply by 2 to account for cells below the midplane
            fsize[k,j]=2*S_val

    if visual:
        ''' Quick plot to check things out '''
        fig,(ax1,ax2)=plt.subplots(1,2,figsize=(15,5))
        for i in range(psize.shape[0]):
            if i==0:
                ax1.plot(r_array,fsize[i],'.',label="%.3f micron"%(psize[i,0]))
            elif i==psize.shape[0]-1:
                ax1.plot(r_array,fsize[i],'.',label="%.3f micron"%(psize[i,0]))
            else:
                if i%10==0:
                    ax1.plot(r_array,fsize[i],'.')

        density_reconstructed=np.zeros(fsize.shape[1])

        for j in range(fsize.shape[1]):
            density_reconstructed[j]=np.sum(np.reshape(fsize[:,j:j+1],fsize.shape[0]))

        ax2.plot(fobj.x,fobj.y,'+',label='original')
        ax2.plot(r_array,density_reconstructed,'.',label='reconstructed')

        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax1.legend()
        ax2.legend()
        plt.show()

    # The gas-to-dust ratio
    if g2d is None:
        g2d=100.0
    else:
        if type(g2d)==list:
            if len(g2d)==1:
                g2d=g2d[0]
            else:
                print("Option no yet available! Try again.")
        elif type(g2d)==str:
            print("Reading g2d from file: %s"%(g2d))
            data_from_file=np.loadtxt(g2d)
            r_from_file=np.reshape(data_from_file[:,0:1],data_from_file.shape[0])
            g2d_from_file=np.reshape(data_from_file[:,1:2],data_from_file.shape[0])
            print(r_from_file)
            print(g2d_from_file)

    sys.exit()

    g2d_array=g2d*np.ones(len(S_array))


    print("\n Working with g2d = ",g2d)


    # Converting to ProDiMo units
    ai_array=(ai_array*u.micron).to(u.cm)
    r_array=(r_array*u.au).to(u.cm)

    # Calling prodimopy
    write("sdprofile.in",r_array.value,S_array*g2d,g2d_array,ai_array.value,fsize)

    return None



def Lfuv(M,R,Mdot,Rin=None):

    '''
    This routine computes the FUV luminosity of the
    star.

    M: stellar mass in Msun
    R: stellar radius in Rsun
    Mdot: mass accretion rate in Msun/yr
    Rin (optional): inner radius of the accretion disk in AU.
                    The default value is Rin=5*R.

    The routine returns the FUV luminosity in Lsun.
    Recall that ProDiMo needs L_fuv/L_star, so
    do not forget to divide!
    '''

    M=M*u.Msun
    R=R*u.Rsun
    Mdot=Mdot*u.Msun/u.yr

    if Rin:
        Rin=(Rin*u.au).to(u.Rsun)
        ratio=R/Rin
    else:
        ratio=1/5.

    '''
    Derive accretion luminosity using eq.(8) of Gullbring et al. 1998
    '''
    L_acc=(c.G*M*Mdot/R*(1-ratio)).to(u.Lsun)


    ''' Computes FUV luminosity using Table 6 of Yang et al. 2012 '''
    pow=-1.670+0.836*np.log10(L_acc.value)
    L_fuv=(10**pow)*u.Lsun

    return L_fuv

def Rhill(ap,Mp,Ms):
    ap=ap*u.au
    Mp=Mp*u.Msun
    Ms=Ms*u.Msun

    Rh=ap*(Mp/(3*Ms))**(1./3)

    return Rh


def iposition(PA_disk,ri,PAi):
    '''
    This routine returns the x,y values for the position of a
    planet to be included in the MCMax3D input file such that the
    planet gets correclty positionated in the output image.

    IMPORTANT: the phi value in the Image.out file must be zero.

    PA_disk: the position angle of the disk measured from
    north to east in deg.
    ri: radial separation (projected) of the planet in AU.
    PAi: position angle (projected) of the planet in deg.
    '''

    # Poision vector of the object
    thetai=((PAi+90)*u.deg).to(u.rad).value
    posi=np.array([ri*np.cos(thetai),ri*np.sin(thetai)])

    # Rotation matrix
    theta=((PA_disk-90)*u.deg).to(u.rad).value
    M=np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])

    # Clockwise rotation by PAi
    pos_rotated=np.dot(M,posi)

    # Converting to MCMax3D coordinates
    x_mcmax=-1.0*pos_rotated[1]
    y_mcmax=pos_rotated[0]

    return (x_mcmax,y_mcmax)


def convert_flux(path_to_file,output_file_name):

    ############################################################
    #
    # This module converts the file ('path_to_file') generated
    # by MCMax3D into a two column file ('output_file_name')
    # with col1: wavelenght in microns and col2: flux in
    # W/m2/micron. The input file is assumed to contain two
    # columns with col1: wavelength in microns and col2:
    # flux in Jy. This is the standar MCMax3D units for the
    # output observables.
    #
    ############################################################

    # Loading MCMax3D output file
    data=np.loadtxt(path_to_file)
    x_array=[]
    y_array=[]
    for i in range(0,data.shape[0]):
        factor=c.c.value/((data[i][0]*1e-6)**2) # (m/s/m^2)
        flux_min=factor*(data[i][1]*1e-26) # from Jy to W/m^3
        x_array.append(data[i][0]) # microns
        y_array.append(flux_min/1e6) # W/m^2/microns
    file=open('%s'%output_file_name,'w')
    for i in range(0,len(x_array)):
        file.write('%.15e %.15e\n'%(x_array[i],y_array[i]))
    return file
