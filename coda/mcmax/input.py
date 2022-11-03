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
                If rlim is None, then it calculates the total mass. [au,au] (tuple).

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


    def plot(self,xlog=None,ylog=None,style=None,
        xlabel=None,ylabel=None):
        x,y=self.x,self.y
        fig,ax=plt.subplots(1,1)

        # Line style
        if style is None:
            ax.plot(x,y,'.')
        else:
            ax.plot(x,y,style)

        # Modifiers
        ax.set_xscale("log")
        ax.set_yscale("log")
        if xlog is False:
            ax.set_xscale("linear")
        if ylog is False:
            ax.set_yscale("linear")
        if xlabel is not None:
            ax.set_xlabel(r"%s"%xlabel)
        if ylabel is not None:
            ax.set_ylabel(r"%s"%ylabel)

        plt.show()


    def rescale_mass(self,k=None,rlim=None,files=None,reeplace=None,type=None,
                    ylim=None):


        """

        Rescale surface density profile by a factor of k.

        Parameters
        ----------

        k           : rescaling factor. [dimensionless] (float).
        rlim        : (rmin,rmax) to only rescale the points with rmin <= r <= rmax.
                    If rlim is None, rescale the entire profile. [au,au] (tuple).
        files       : This performs a point-by-point rescaling of the density profile
                    according to the ratio of R=file[0]/file[1]. It uses a linear
                    interpolator to infer the values of the files at each self.x point.
                    If out of bounds interpolation is required, it uses the extreme
                    values of the files as boundary conditions. (list).
        reeplace    : It has to be used alongside k and rlim. If True, it reeplaces
                    all the values between rlim by k.
        type        :

        """

        x=(self.x*u.au).to(u.cm).value
        y=self.y


        Mold=(2*np.pi*simps(x*y,x)*u.g).to(u.Msun)


        # Rescale the whole profile by a constant factor
        if k is not None and rlim is None:
            ynew=k*y


        # Rescale between rmin and rmax by a constant factor
        if k is not None and rlim is not None:
            """
            try:
                if reeplace is False:
                    rmin,rmax=rlim
                    ynew=[]
                    for i in range(len(self.x)):
                        if self.x[i]>=rmin and self.x[i]<=rmax:
                            ynew+=[k*y[i]]
                        else:
                            ynew+=[1*y[i]]

                elif reeplace is True:
                    rmin,rmax=rlim
                    ynew=[]
                    for i in range(len(self.x)):
                        if self.x[i]>=rmin and self.x[i]<=rmax:
                            ynew+=[k]
                        else:
                            ynew+=[1*y[i]]
            except:
                print("No 'reeplace' found. Assuming reeplace = False")
            """
            if reeplace is True:
                rmin,rmax=rlim
                ynew=[]
                for i in range(len(self.x)):
                    if self.x[i]>=rmin and self.x[i]<=rmax:
                        ynew+=[k]
                    else:
                        ynew+=[1*y[i]]

            else:
                rmin,rmax=rlim
                ynew=[]
                for i in range(len(self.x)):
                    if self.x[i]>=rmin and self.x[i]<=rmax:
                        ynew+=[k*y[i]]
                    else:
                        ynew+=[1*y[i]]

        # Point-by-point rescaling using two input files
        if files is not None:

            # Load data
            data_file_1 = np.loadtxt(files[0])
            data_file_2 = np.loadtxt(files[1])

            x_file_1 = np.reshape(data_file_1[:,0],data_file_1.shape[0])
            y_file_1 = np.reshape(data_file_1[:,1],data_file_1.shape[0])
            x_file_2 = np.reshape(data_file_2[:,0],data_file_2.shape[0])
            y_file_2 = np.reshape(data_file_2[:,1],data_file_2.shape[0])

            # Interpolate files (linearly) onto model's grid
            y_file_1_interp = np.interp(self.x, x_file_1, y_file_1)
            y_file_2_interp = np.interp(self.x, x_file_2, y_file_2)

            # Do a check?
            ms=1
            plt.plot(x_file_1,y_file_1,'+',markersize=5,label='$%s$ input'%files[0],color='black')
            plt.plot(self.x,y_file_1_interp,'-',markersize=ms,label='$%s$ interpolated'%files[0],color='blue')

            plt.plot(x_file_2,y_file_2,'+',markersize=5,label='$%s$ input'%files[1],color='red')
            plt.plot(self.x,y_file_2_interp,'-',markersize=ms,label='$%s$ interpolated'%files[1],color='green')

            #plt.plot(self.x,self.y,label='model',color='red')
            plt.legend(frameon=False)
            plt.show()

            # Built R coefficients
            R_array = abs(y_file_1_interp/y_file_2_interp)

            # Built new density profile
            ynew=[i*j for i,j in zip(self.y,R_array)]

            # Plot final result
            plt.plot(self.x,self.y,label='initial density')
            plt.plot(self.x,ynew,label='corrected density')
            plt.legend()

            plt.yscale('log')
            plt.show()

        # Reeplace points by a straight line
        if type=='linear':
            print(rlim)
            rmin,rmax   = rlim
            ymin,ymax   = ylim
            m           = (ymax-ymin)/(rmax-rmin)
            b           = ymin-m*rmin

            ynew=[]
            for i in range(len(self.x)):
                if self.x[i]>=rmin and self.x[i]<=rmax:
                    ynew+=[m*self.x[i]+b]
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

        return None


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


def make_density_profile(Rin,Rout,Mdust,
                        Rtap=None,epsilon=None,gamma=None,write=None):

    """

    Creates a surface density profile for a given value of the mass assuming an
    exponential and tappered profile (see Portilla-Revelo et al. 2022, Eq.(1))

    Parameters
    ----------
    Rin     : inner disk radius (au)
    Rout    : outer disk radius (au)
    Mdust   : Mass of the disk (Msun)
    R_tap   : tappering radius of the disk (au)
    epsilon : exponent of the linear-decay part
    gamma   : exponent of the exponentia-decay part

    Output
    ------
    The array for the radial distance in au and the array for the surface density
    in g/cm2. If write=True, it saves to a file 'density.dat'. If plot=True,
    it plots the profile.

    """

    # Number of sampling points
    k = 10
    N_points = 2**k+1

    # Dealing with input arguments
    Rin = Rin*u.au
    Rout = Rout*u.au
    Mdust = Mdust*u.Msun
    Rtap=Rtap*u.au if Rtap is not None else 100*u.au
    if epsilon is None:epsilon=1.0
    if gamma is None:gamma=1.0

    # Array of radial points
    rarray=np.linspace(Rin,Rout,N_points)

    # Finding Sigma0
    def integrand(x):
        return x * (Rtap.value/x)**epsilon * np.exp(-(x/Rtap.value)**(2-gamma))
    Sigma0=(Mdust / (2*np.pi * quad(integrand,Rin.value,Rout.value)[0]*u.au**2))
    Sigma0=Sigma0.to(u.g/u.cm**2)

    # Sigma array
    Sigma=[]
    for i in range(len(rarray)):
        r=rarray[i].value
        Sigma.append( Sigma0.value * (Rtap.value/r)**epsilon * np.exp(-(r/Rtap.value)**(2-gamma)) )
    Sigma=np.array(Sigma)
    Sigma=Sigma*(u.g/u.cm**2)

    # Write to a file
    if write:
        file=open("density.dat","w")
        for i,j in zip(rarray,Sigma):
            file.write("%.15e %.15e\n"%(i.value,j.value))
        file.close()

    return rarray,Sigma


def convert_comp(fc,porosity,qtype) :

    '''

    Convert grain composition from the MCMax3D standard to the ProDiMo standard
    and vice-versa.

    Parameters
    ----------
    fc: fraction of carbon (either per volume or per mass, depending on your needs).

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

        print("\nThe correspondent ProDiMo quantities (per volume) are:")
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

        print("\nThe correspondent MCMax3D quantities (per mass) are:")
        print("fcM=%.15f"%(fcM))
        print("fsiM=%.15f"%(fsiM))
        print("fsiM+fcM=%.15f"%(fsiM+fcM))

    return None


def convert_density_file(model,Nzones,g2d=None,visual=None,find_dust_mass=None,leftg2d=None):

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
    g2d     : The gas-to-dust ratio. If g2d is None assume to be 100
            everywhere. If g2d is a list and len(g2d)=1, use a constant g2d
            ratio of g2d[0]. If g2d is a two-column file containing the g2d ratio
            (second column) in function of the distance in AU (first column), then
            the g2d ratio will be interpolated at each radial point in the grid.
            A file containing the interpolated g2d ratio profile will be created.
    visual  : If true, shows the computed density profile for a selected dust
            sizes and also the reconstructed profile compared to the original one.
            (Bool).
    leftg2d : This sets the value of the g2d at those radial grid points to the
            left of the minimum distance included in the g2d file. (float)

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
    print("\n Retrieving dust size distribution from MCMax3D model... \n")
    psize=get_size_distribution(model)
    ai_array=psize[:,0] # microns


    # Reading input file to get zone limits
    Rin_array=[]
    Rout_array=[]
    infile=open(model+"/input.dat").readlines()

    for i in range(1,Nzones+1):

        for line in infile:

            if line.split("=")[0]==("zone%d:Rin"%(i)):
                Rin=float(line.split("=")[1])
                Rin_array.append(Rin)

            if line.split("=")[0]==("zone%d:Rout"%(i)):
                Rout=float(line.split("=")[1])
                Rout_array.append(Rout)


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

        """
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
        """
        # Correcting for radius at the inner and outer edge of the PPD
        if zoneID==1:
            rmidAU_min=Rin_array[zoneID-1]
            rmid_min=(rmidAU_min*u.au).to(u.cm).value
            rmidAU[0]=rmidAU_min
            rmid[0]=rmid_min

        if zoneID==Nzones:
            rmidAU_max=Rout_array[-1]
            rmid_max=(rmidAU_max*u.au).to(u.cm).value
            rmidAU[-1]=rmidAU_max
            rmid[-1]=rmid_max


        # Correcting for theta
        theta_midplane=(90.0*u.deg).to(u.rad).value
        tmid[int(len(tmid)*0.5-1)]=theta_midplane


        # Declare matrix required for interpolation
        # Dim: len(theta)/2 x len(radius). Each entry is an instance of the
        # CellMcmax class.
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


    # Read and concatenate fsizes and rmidAU
    for i in range(1,Nzones+1):
        if i==1:
            r_array=func_fsize(i)[0]
            M_mgc_full=func_fsize(i)[1]
        else:
            r_array=np.concatenate((r_array,func_fsize(i)[0]))
            M_mgc_full=np.concatenate((M_mgc_full,func_fsize(i)[1]),axis=1)

    """
    rmidAU_1, M_mgc_1 = func_fsize(1)[0], func_fsize(1)[1]
    rmidAU_2, M_mgc_2 = func_fsize(2)[0], func_fsize(2)[1]

    r_array=np.concatenate((rmidAU_1,rmidAU_2))
    M_mgc_full=np.concatenate((M_mgc_1,M_mgc_2),axis=1)
    """


    # Read ProDiMo grid. This steps requires prodimopy!
    model_prodimo=read_prodimo()
    r_array_p=np.reshape(model_prodimo.x[:,0:1],model_prodimo.x.shape[0])   # au
    z_matrix=model_prodimo.z                                                # au, dim:len(r)*len(z)


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
    print("\n Following are the coordinates for the cell's center (for MCMax3D)!")
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
        # Quick plot to check things out
        fig,(ax1,ax2)=plt.subplots(1,2,figsize=(15,5))
        fig.suptitle('Dust column density reconstruction')
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
        ax2.plot(r_array,density_reconstructed,'-',label='reconstructed')

        ax1.set_yscale("log")
        ax2.set_yscale("log")

        ax1.legend(frameon=False)
        ax2.legend(frameon=False)

        ax1.set_xlabel("distance (au)")
        ax2.set_xlabel("distance (au)")
        ax1.set_ylabel("dust surface density (g/cm2)")
        ax2.set_ylabel("dust surface density (g/cm2)")

        plt.show()

    # The gas-to-dust ratio
    if g2d is None:
        g2d=100.0
        g2d_array=g2d*np.ones(len(S_array))
    else:
        if type(g2d)==list:
            if len(g2d)==1:
                g2d=g2d[0]
                g2d_array=g2d*np.ones(len(S_array))
            else:
                print("Option no yet available! Try again.")
        elif type(g2d)==str:

            print("Reading g2d from file: %s"%(g2d))
            data_from_file=np.loadtxt(g2d)
            r_from_file=np.reshape(data_from_file[:,0:1],data_from_file.shape[0])
            g2d_from_file=np.reshape(data_from_file[:,1:2],data_from_file.shape[0])

            if leftg2d is not None:
                g2d_array=np.interp(r_array,r_from_file,g2d_from_file,left=leftg2d) # Linear interpolation
            else:
                g2d_array=np.interp(r_array,r_from_file,g2d_from_file) # Linear interpolation

            # Creating output file for interpolated g2d profile
            fg2d=open(g2d+".interpolated","w")
            for i,j in zip(r_array,g2d_array):
                fg2d.write("%.5f    %.5f\n"%(i,j))
            fg2d.close()

    # Converting to ProDiMo units
    ai_array=(ai_array*u.micron).to(u.cm)
    r_array=(r_array*u.au).to(u.cm)

    # Calling prodimopy
    write("sdprofile.in",r_array.value,S_array*g2d_array,g2d_array,ai_array.value,fsize)

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

    """

    Computes the Hill's radius around a planet.

    Parameters
    ----------
    ap  : The planet's semimajor axis (au).
    Mp  : The planet's mass (Msun)
    Ms  : The star's mass (Msun)

    """

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

    Parameters
    ----------
    PA_disk : the position angle of the disk measured from north to east. (deg)
    ri      : radial separation (projected) of the planet. (au)
    PAi     : position angle (projected) of the planet. (deg)

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


def convert_flux(infile):

    """

    Converts units of Monte Carlo spectrum computed with MCMax 3D.

    Parameters
    ----------
    infile  : path to the MCSpec file. It contains two columns: wavelength (micron)
            and Flux (Jy). Those are the standart output units from an MCMax3D
            spectrum.

    Output
    ------
    Two column file named infile.converted with two columnds: wavelength (micron)
    and Flux (W/m2/micron)

    """

    # Loading MCMax3D output file
    data=np.loadtxt(infile)
    x_array=[]
    y_array=[]
    for i in range(0,data.shape[0]):
        factor=c.c.value/((data[i][0]*1e-6)**2) # (m/s/m^2)
        flux_min=factor*(data[i][1]*1e-26) # from Jy to W/m^3
        x_array.append(data[i][0]) # microns
        y_array.append(flux_min/1e6) # W/m^2/microns
    file=open('%s.converted'%infile,'w')
    for i in range(0,len(x_array)):
        file.write('%.15e %.15e\n'%(x_array[i],y_array[i]))
    return file


def Fnu_disk(r,wl,Temp,Sigma,kabs,d,inc):

    """

    Computes the monochromatic flux from a vertically isothermal disk in LTE
    assuming no background intensity.

    Parameters
    ----------
    d       : distance to the source (pc) [float]
    inc     : inclination of the disk respect to the plane of the sky (deg) [float]
    kabs    : absorption mass opacity at that wavelength (cm2/g)
    Temp    : temperature at each radial point (K) [array]
    Sigma   : surface density at each radial point (g/cm2) [array]
    wl      : wavelength at which the flux is computed (micron) [float]
    r       : distance from the star (au) [array]

    Output
    ------
    Monochromatic flux density in Jy

    """

    # Dealing with input arguments
    wl = wl*u.micron
    Temp = Temp*u.K
    Sigma = Sigma*(u.g/u.cm**2)
    kabs = kabs*(u.cm**2/u.g)
    d = (d*u.pc).to(u.cm)
    inc = (inc*u.deg).to(u.rad)
    r = (r*u.au).to(u.cm)

    # Computing the source function
    bb=models.BlackBody(temperature=Temp)
    bb_lambda=bb(wl)                                       # erg / (cm2 Hz s sr)

    # Compute optical depth
    tau=Sigma*kabs

    # Computing flux
    Rin = r[0].to(u.au)
    Rout = r[-1].to(u.au)
    integrand = r * bb_lambda * ( 1 - np.exp(-tau/np.cos(inc)) )
    print("\n Integrating flux between %.3f au and %.3f au"%(Rin.value,Rout.value))
    I = simps(integrand,r)
    mono_flux = 2*np.pi * np.cos(inc)/(d.value)**2 * I     # erg / (cm2 Hz s)

    # Manipulation of flux units
    mono_flux = ( mono_flux*( u.erg/(u.cm**2 * u.Hz * u.s) ) ).to(u.Jy)

    return mono_flux


def brightness_temperature_equiv(I,nu0,bmaj,bmin):

    """

    Converts intensity from mJy/beam to K.

    Parameters
    ----------
    I       : flux. (mJy/beam) [float].
    nu0     : rest frequency. (Hz) [float].
    bmaj    : beam's major axis. (arcsec) [float].
    bmin    : beam's minor axis. (arcsec) [float].

    """

    I=I/1e3                                     # Flux in Jy/beam
    nu0=nu0/1e9                                 # nu0 in GHz
    value=1.222e6 * I/(nu0**2 * bmaj * bmin)    # Tb in K

    return value


def thin_scale(d,Teq,wl,kabs,inc,Mdust=None,flux=None):

    """

    Parameters
    ----------
    d       : Distance to source. (pc) [float].
    inc     : Inclination of the source. (deg) [float].
    Teq     : Dust equilibrium temperature. (K) [float].
    wl      : Wavelength. (micron) [float].
    kabs    : Absportion opacity. (cm^2/g) [float].
    Mdust   : Disk's mass. (Msun) [float].
    flux    : Sub-mm flux. (Jy) [float].

    """

    # Handling units
    d       = d*u.pc
    inc     = (inc*u.deg).to(u.rad)
    Teq     = Teq*u.K
    wl      = wl*u.micron
    kabs    = kabs*(u.cm**2/u.g)

    # Black body model
    bb=models.BlackBody(temperature=Teq)

    if Mdust is not None:

        print("\n Computing flux \n")
        Mdust   = Mdust*u.Msun
        value   = Mdust*kabs*bb(wl)/(d**2*np.cos(inc))


    if flux is not None:

        print("\n Computing mass \n")
        flux    = flux*u.Jy
        value   = flux*d**2*np.cos(inc)/(kabs* (bb(wl)*1*u.sr) )
        value   = value.to(u.Msun)

    return value
