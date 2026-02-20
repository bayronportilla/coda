from __future__ import annotations
from coda.mcmax.header import *
from dataclasses import dataclass
import numpy as np
from astropy import units as u
from scipy.integrate import romb,simps,quad,simpson

@dataclass
class File:

    x:float
    y:float
    name:str


    def calculate_mass(self,rlim=None):

        """

        Integrate the density profile to find the mass (for gas or dust).

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
            mass=(2*np.pi*simps(xint*yint,xint)*u.g).to(u.Msun)
        else:
            mass=(2*np.pi*simps(x*y,x)*u.g).to(u.Msun)

        return mass


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
                    ylim=None,mass=None,epsilon=None,gamma=None,sigma=None,pivot=None):


        """

        Rescale surface density profile by a factor of k.

        Parameters
        ----------

        k           : rescaling factor. [dimensionless] (float).
        rlim        : (rmin,rmax) to only rescale the points with rmin <= r <= rmax.
                    If rlim is None, rescale the entire profile. [au,au] (tuple).
        files       : This performs a point-by-point rescaling of the density profile
                    according to the ratio of R=file[0]/file[1]. This is useful, for 
                    example, when file[0] and file[1] are intensity profiles from an 
                    optically thin tracer. The routine uses a linear interpolator to 
                    infer the values of the files at each self.x point. If out of bounds, 
                    interpolation will be required and it will use the extreme values 
                    of the files as boundary conditions. User has to make sure that the 
                    units in both columns of both files are the same. (list).
        reeplace    : It must be used alongside k and rlim. If True, it reeplaces
                    all the values within rlim by k.
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


        # Point-by-point rescale using input files
        if files is not None:
            # Load data
            file1 = np.loadtxt(files[0])
            file2 = np.loadtxt(files[1])

            x1 = np.reshape(file1[:,0],file1.shape[0])
            y1 = np.reshape(file1[:,1],file1.shape[0])
            x2 = np.reshape(file2[:,0],file2.shape[0])
            y2 = np.reshape(file2[:,1],file2.shape[0])

            # Interpolate files (linearly) onto model's grid
            y1_interp = np.interp(self.x, x1, y1)
            y2_interp = np.interp(self.x, x2, y2)

            # Do a check?
            ms=1
            plt.plot(x1,y1,'+',markersize=5,label='$%s$ input'%files[0],color='black')
            plt.plot(self.x,y1_interp,'-',markersize=ms,label='$%s$ interpolated'%files[0],color='blue')

            plt.plot(x2,y2,'+',markersize=5,label='$%s$ input'%files[1],color='red')
            plt.plot(self.x,y2_interp,'-',markersize=ms,label='$%s$ interpolated'%files[1],color='green')

            #plt.plot(self.x,self.y,label='model',color='red')
            plt.legend(frameon=False)
            plt.show()

            # Built R coefficients
            R_array = abs(y1_interp/y2_interp)

            # Built new density profile
            ynew=[i*j for i,j in zip(self.y,R_array)]

            print(len(self.x))
            print(len(ynew))

            # Plot final result
            plt.plot(self.x,self.y,label='initial density')
            plt.plot(self.x,ynew,label='corrected density')
            plt.legend()

            plt.yscale('log')
            plt.show()
        

        # Reeplace points by a straight line
        if type=='linear':

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
        '''
        # Reeplace points by a straight line
        if type=='linear':

            rmin,rmax   = rlim
            ymin,ymax   = ylim

            if pivot=='min':
                m  = (k-1.0)/(rmax-rmin)
                b  = 1.0-m*rmin

            elif pivot=='max':
                m  = (1.0-k)/(rmax-rmin)
                b  = 1.0-m*rmax

            """
            rlinear = np.linspace(rmin,rmax,100)
            ylinear = m*rlinear+b

            plt.plot(rlinear,ylinear)
            plt.show()
            """

            ynew=[]
            for i in range(len(self.x)):
                if self.x[i]>=rmin and self.x[i]<=rmax:
                    ynew+=[y[i]*(m*self.x[i]+b)]
                else:
                    ynew+=[1*y[i]]
        '''

        if type=='exponential':

            try:
                rmin,rmax   = rlim
                profile     = make_density_profile(rmin,rmax,mass,epsilon=epsilon,gamma=gamma)
                xprof       = profile[0].value
                yprof       = profile[1].value

                ynew=[]
                for i in range(len(self.x)):
                    if self.x[i]>=rmin and self.x[i]<=rmax:
                        ynew+=[np.interp(self.x[i],xprof,yprof)]
                    else:
                        ynew+=[1*y[i]]

            except:
                print('\n mass argument must be passed')

        # Gaussian: usually helpful when one of the pivots has y~0
        if type=='gaussian':

            rmin,rmax       = rlim
            ymin,ymax       = ylim

            r_gauss     = np.linspace(rmin,rmax,100)

            if pivot=='min':
                y_gauss = norm(loc=rmin,scale=sigma).pdf(r_gauss)
            elif pivot=='max':
                y_gauss = norm(loc=rmax,scale=sigma).pdf(r_gauss)
            y_gauss     = y_gauss/y_gauss.max()

            # Normalize to ymax and then multiply by gaussian PDF.
            ynew=[]
            for i in range(len(self.x)):
                if self.x[i]>=rmin and self.x[i]<=rmax:
                    if pivot=='min':
                        ynew+=[ymin*np.interp(self.x[i],r_gauss,y_gauss)]
                    elif pivot=='max':
                        ynew+=[ymax*np.interp(self.x[i],r_gauss,y_gauss)]
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
    
    def shift_profile(self,Dx):

        """
        Shifts a surface density file along the x-axis by a fixed amount 
        conserving the mass

        Parameters
        ----------

        Dx  : amount along the x-axis the profile will be shifted by. [float] (au)

        """

        # Shifting along x-axis
        xnew = self.x+Dx # au
        xnew_cm = ( xnew*u.au).to(u.cm).value

        # Computing masses
        Mold = self.calculate_mass()
        Mnew = (2*np.pi*simps(xnew_cm*self.y,xnew_cm)*u.g).to(u.Msun)
        
        # Rescaling density to conserve mass
        ynew = self.y*( Mold/Mnew )

        # Save to file
        f=open(self.name+".shifted","w")
        for i,j in zip(xnew,ynew):
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
    filenames=np.sort(filenames).tolist()

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


def make_density_profile(Rin,Rout,M,
                         Rtap=None,
                         epsilon=None,
                         gamma=None,
                         write=None,
                         Sigmain=None):

    """

    Creates a surface density profile for a given value of the mass assuming an
    exponential and tappered profile (see e.g., Eq. 1 in Portilla-Revelo et al. 2022)

    Parameters
    ----------
    Rin     : inner disk radius (au)
    Rout    : outer disk radius (au)
    M       : Mass of the disk (Msun)
    R_tap   : tappering radius of the disk (au)
    epsilon : exponent of the linear-decay part
    gamma   : exponent of the exponentia-decay part
    Sigmain : Boundary condition for the density at Rin (g/cm2)
 
    Output
    ------
    The array for the radial distance in au and the array for the surface density
    in g/cm2. If write=True, it saves to a file 'density.dat'. If plot=True,
    it plots the profile.

    """

    # Number of sampling points
    k           = 10
    N_points    = 2**k+1

    # Dealing with input arguments
    Rin     = Rin*u.au
    Rout    = Rout*u.au
    M       = M*u.Msun
    Rtap    = Rtap*u.au if Rtap is not None else 100*u.au

    if epsilon is None:epsilon=1.0
    if gamma is None:gamma=1.0

    # Array of radial points
    rarray=np.linspace(Rin,Rout,N_points)

    # Finding Sigma0
    def integrand(x):
        return x * (Rtap.value/x)**epsilon * np.exp(-(x/Rtap.value)**(2-gamma))
    Sigma0=(M / (2*np.pi * quad(integrand,Rin.value,Rout.value)[0]*u.au**2))
    Sigma0=Sigma0.to(u.g/u.cm**2)

    # Sigma array
    Sigma = []
    for i in range(len(rarray)):
        r = rarray[i].value
        Sigma.append( Sigma0.value * (Rtap.value/r)**epsilon * np.exp(-(r/Rtap.value)**(2-gamma)) )
    Sigma = np.array(Sigma)
    Sigma = Sigma*(u.g/u.cm**2)

    # Write to a file
    if write:
        file=open("density.dat","w")
        for i,j in zip(rarray,Sigma):
            file.write("%.15e %.15e\n"%(i.value,j.value))
        file.close()

    # If Sigmain is given
    if Sigmain is not None:

        Sigmain = Sigmain*(u.g/u.cm**2)
        Sigma0  = Sigmain/((Rtap/Rin)**epsilon * np.exp(-(Rin/Rtap)**(2-gamma)))
        Mass    = (Sigma0 * 2*np.pi * quad(integrand,Rin.value,Rout.value)[0]*u.au**2).to(u.Msun)

        Sigma = []
        for i in range(len(rarray)):
            r = rarray[i].value
            Sigma.append( Sigma0.value * (Rtap.value/r)**epsilon * np.exp(-(r/Rtap.value)**(2-gamma)) )
        
        Sigma = np.array(Sigma)
        Sigma = Sigma*(u.g/u.cm**2)
        
        print(Mass)
        
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


def convert_density_file(model,
                         Nzones,
                         g2d=None,
                         visual=None,
                         find_dust_mass=None,
                         leftg2d=None,
                         save=None):

    '''
    Convert an MCMax3D density field into a 1D sdfile.in for ProDiMo.

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
    *   The length of the r_array determines the resolution of the ProDiMo
        model. Choose it wisely!

    *   This routine only works when the MCMax3D model has *one* type of
        particle; this is, there is only one 'computepartXX' keyword
        in the input.dat file.

    *   Before using this routine, you must run a 'stop-after-init'
        ProDiMo simulation. This is because there must exist a grid where to 
        interpolate the MCMax3D quntities. Make sure this simulation has the
        interface1D keyword set to False and the sdprofile.in is commented.

    *   The density profile inside the folder must be named as:
        'surface_density_PDS70_70_cropped.dat'

    *   This routine works only with 2D models.

    *   MCMax3D model must have two zones---not one, not three...

    *   Rin and Rout must be equal to 0.04 au and 130 au respectively.

    *   The number of radial grid points (NXX) of the ProDiMo model must be 
        equal to its MCMax3D counterpart.

    *   For production results, it is advisable to run the MCMax3D model with a
        refined radial grid, Nr>=60 per zone, as long as this is computationally
        feasible.

    Parameters
    ----------
    model   : path to the MCMax3D model directory. (str).
    g2d     : The gas-to-dust ratio. If g2d is None, assume to be 100
            everywhere. If g2d is a list and len(g2d)=1, use a constant g2d
            ratio of g2d[0]. If g2d is a two-column file containing the g2d ratio
            (second column) as a function of the distance in AU (first column), then
            the g2d ratio will be interpolated at each radial point on the grid.
            A file containing the interpolated g2d ratio profile will be created.
    visual  : If true, shows the computed density profile for a selected dust
            sizes and also the reconstructed profile compared to the original one.
            (Bool).
    leftg2d : This sets the value of the g2d at those radial grid points to the
            left of the minimum distance included in the g2d file. (float)
    save    : Save to a fits file everything that is needed to create the sdprofile.in
            file. This way, it won't be needed to access the MCMax3D model 
            if one wants to try, for example, a different gas distribution while keeping
            the dust distribution unchanged. The file has 3 HDUs: hdu1=fsize,
            hdu2=r_array (cm), and hdu3=ai_array (cm).  

    Example
    -------
    1. Run a stop_after_init ProDiMo model to create a cylindrical grid if needed.
    Remember to run with interface_1D=.false. Remember also to set the Nxx variable
    to match the total number of radial points of the MCMax3D model.
    >>> cd <path-to-ProDiMo-model>
    >>> nano Parameter.in
    >>> prodimo

    2. In the ProDiMo model directory, open an ipython session and import
    the coda modules
    >>> ipython
    >>> from coda.mcmax import input

    3. Call the module passing the appropriate parameters
    >>> input.convert_density_file("<path-to-MCMax3D-model>",visual=True,find_dust_mass=False,Nzones=2)

    '''

    @dataclass
    class CellMcmax:
        r           : float             # cm
        theta       : float
        phi         : float
        r_max       : float
        r_min       : float
        theta_max   : float
        theta_min   : float
        phi_max     : float
        phi_min     : float
        comp        : float

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


    # Saving to a fits file
    if save is True:
        hdu1=fits.PrimaryHDU(data=fsize) # Primary HDU contains the fsize matrix
        hdu2=fits.ImageHDU(data=r_array.value) # cm
        hdu3=fits.ImageHDU(data=ai_array.value) # cm
        
        hdul=fits.HDUList([hdu1,hdu2,hdu3])
        hdul.writeto('sdprofile.fits',overwrite=True)
    

    # Calling prodimopy
    write("sdprofile.in",r_array.value,
          S_array*g2d_array,
          g2d_array,
          ai_array.value,
          fsize)

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

    This rountine computes mass from flux and vice-versa assuming optically
    thin emission.

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


def Tdust(Tsource,Rsource,d,nu,kabs):

    """

    Computes the equilibrium temperature of a distribution of dust grains assuming
    radiation heating and black-body cooling. The illumination sources is a black
    body with effective temperature Tsource and size Rsource. It is assumed that
    the location where the equilibrium temperature is being calculated is r >> Rsource.

    Paramters:
    ----------
    Tsource         : Effective temperature of the source of illumination. (K) [float].
    Rsource         : Radius of the illumination source. (Rsun) [float].
    d               : distance from the illumination source. (au) [float].
    nu              : frequency array over which the opacity is defined. (micron) [array].
    kabs            : Absoprtion opacity array. (cm2/g) [array].

    Ouput:
    ------
    Equilibriun temperature in K

    """


    # Absorption opacity averaged over black body with effective temperature T
    def kabs_ave(T):

        # Numerical integration
        bb          = models.BlackBody(temperature=T*u.K)                               # erg / cm2 Hz s sr
        y           = bb(nu*u.Hz)*(kabs*(u.cm**2/u.g))                                  # erg / Hz s sr g
        I           = simpson(y,nu)                                                     # Dimensionless
        num         = I*(u.erg/(u.s*u.sr*u.g))                                          # erg / s sr g
        den         = (c.sigma_sb.cgs.value/np.pi * T**4)*(u.erg/(u.cm**2*u.s*u.sr ))   # erg / cm2 s sr
        kabs_ave    = num/den                                                           # cm2 / g

        return kabs_ave


    # The equation to solve
    def f(x):
        frad = (0.5*(Rsource*u.Rsun)/(d*u.au)).decompose()
        return x**4 - frad**2 * kabs_ave(Tsource)/kabs_ave(x) * Tsource**4


    # Using Newton-Raphson method to find equilibrum temperature
    Td_sol = optimize.newton(f,100)

    return Td_sol




def Td_razor(Tsource,Rsource,d):

    """
    Finds the dust temperature of a razor disk irradiated by a central star.
    See Eq.(2.38) Armitage

    Parameters:
    -----------
    Tsource         : Effective temperature of the source of illumination. (K) [float].
    Rsource         : Radius of the illumination source. (Rsun) [float].
    d               : distance from the illumination source. (au) [float].

    Output:
    -------

    """

    Rsource     = Rsource*u.Rsun
    d           = d*u.au
    R2d         = (Rsource/d).decompose().value

    Tdisk = ((1/np.pi)*(np.arcsin(R2d) - (R2d)*(1-(R2d)**2)**0.5))**0.25 * Tsource

    return Tdisk


def Td_grain(Tsource,Rsource,d):

    """
    Finds the dust temperature of a single grain irradiated by a central star.

    Parameters:
    -----------
    Tsource         : Effective temperature of the source of illumination. (K) [float].
    Rsource         : Radius of the illumination source. (Rsun) [float].
    d               : distance from the illumination source. (au) [float].

    Output:
    -------

    """


    Lsource     = (4*np.pi*c.sigma_sb.cgs*(Rsource*u.Rsun)**2*(Tsource*u.K)**4)
    Tgrain      = ((Lsource/(4*np.pi*c.sigma_sb.cgs*(d*u.au)**2))**0.25).to(u.K)

    return Tgrain


def h_gas(r,
    hprm_list=None):

    """

    Gas scale height

    Parameters:
    -----------
    r           : radial distance. (au)
    prm_list    : parameter list. [r0(au),h0(au),beta]

    Output:
    -------

    """

    if hprm_list is None:
        hprm_list = [1.0, 0.1, 1.2]

    r0      = hprm_list[0]
    h0      = hprm_list[1]
    beta    = hprm_list[2]

    value = h0*(r/r0)**beta

    return value


def hd_with_settling(agrain,r,Sigmag,
    alpha=None,
    hprm_list=None,
    rho_grain=None):

    """

    Dust scale height per bin size

    Parameters:
    -----------
    agrain      : grain size. (micron)
    r           : radial distance. (au)
    Sigmag      : Gas surface density at r. (g/cm2)
    alpha       : (Optional) Turbulence strength.
    hprm_list   : (Optional) Parameter list for h_gas(). [r0(au),h0(au),beta]
    rho_grain   : (Optional) Gran mass density. (g/cm3)

    Output:
    -------
    hd          : Scale height for i-th grain. (au)

    """

    # Optional variables
    if hprm_list is None:
        hprm_list = [1.0, 0.1, 1.2]

    if rho_grain is None:
        rho_grain = 2

    if alpha is None:
        alpha = 1e-3

    # Asign units
    s           = agrain*(u.micron)
    rho_grain   = rho_grain*(u.g/u.cm**3)
    Sigmag      = Sigmag*(u.g/u.cm**2)

    # Gas scale height
    r0      = hprm_list[0]
    h0      = hprm_list[1]
    beta    = hprm_list[2]
    hg      = h0*(r/r0)**beta

    # Calculate f factor at point r for bin size s
    f_rs        = (alpha/((6*np.pi)**0.5*s*rho_grain) * Sigmag).decompose()

    # Calculate scale height for size s
    hd_rs = hg * (f_rs/(1+f_rs))**0.5

    return hd_rs


def compare(pdm: prodimopy.read.read_prodimo, 
            rdm: radmc3dPy.analyze.readData):

    xpdm = pdm.x[:,0]
    xrdm = (rdm.grid.x*u.cm).to(u.au).value


    """ Check dust mass """
    def check_dust_mass():
    
        mdust_pdm = (np.sum(pdm.rhod*pdm.vol)*u.g).to(u.Msun).value
        mdust_rdm = (rdm.getDustMass()*u.g).to(u.Msun).value
        
        err = 100 * ( abs(mdust_pdm - mdust_rdm)/mdust_pdm )

        print("Dust mass prodimo (Msun): ", mdust_pdm)
        print("Dust mass radmc3d (Msun): ", mdust_rdm)
        print(f"Relative error |mdust_pro-mdust_rad|/mdust_pro x 100 = {err:.2f}%")

        return None

    #check_dust_mass()


    """ Check surface density profile """
    def check_sigmad():

        rdm.getSigmaDust()

        sigmad_pdm = 2*pdm.sdd[:,0]
        sigmad_rdm = rdm.sigmadust[:,0]

        # Find shared spatial grid
        xmin = max(xpdm.min(),xrdm.min())
        xmax = min(xpdm.max(),xrdm.max())
        xcommon = np.linspace(xmin,xmax,1000)

        # Interpolate to common axis
        sigmad_pdm_int = np.interp(xcommon,xpdm,sigmad_pdm)  
        sigmad_rdm_int = np.interp(xcommon,xrdm,sigmad_rdm) 

        # Find relative error
        err = (sigmad_pdm_int-sigmad_rdm_int)/sigmad_pdm_int

        # Plot
        fig,(ax1,ax2) = plt.subplots(nrows=2)
        ax1.loglog(xpdm,sigmad_pdm,'.',label='prodimo')
        ax1.loglog(xrdm,sigmad_rdm,'+',label='radmc3d')
        ax1.set_ylabel("sigmad (g/cm2)")
        ax1.legend()
        
        ax2.plot(xcommon,err)
        ax2.set_xlabel("r (au)")
        ax2.axhspan(-0.1,0.1,color='grey',alpha=0.25)
        ax2.set_ylim(-1,1)
        ax2.set_ylabel("err")

        fig.suptitle("Check surface density")
        
        plt.show()
        
        return None
    
    #check_sigmad()

    
    def interpolate_rhod_rdm_onto_pdm():
        '''
        Docstring for spherical_to_cylindrical
        '''

        rs_rdm = rdm.grid.x # array of radii (cm)
        ts_rdm = rdm.grid.y # array of thetas (rad)
        ps_rdm = rdm.grid.z # array of phis (rad)

        if len(ps_rdm) == 1:
            # Create 2D interpolator to map rdm.rhodust onto cylindrical meshgrid            
            interp = RegularGridInterpolator((rs_rdm,ts_rdm),rdm.rhodust[:,:,0,0],bounds_error=False)

            # Get prodimo grid
            xs_pdm = (pdm.x*u.au).to(u.cm).value
            zs_pdm = (pdm.z*u.au).to(u.cm).value 

            # Convert ProDiMo cylyndrical grid onto radmc3d spherical grid            
            rs_pdm = (xs_pdm**2 + zs_pdm**2)**0.5
            ts_pdm = np.arctan(xs_pdm/zs_pdm)

            # Define array for radmc3d rhodust values interpolated on ProDiMo grid
            rhodust_rdm_onto_pdm = np.zeros(pdm.rhod.shape)            
            
            # Interpolate dust density onto ProDiMo grid
            count=0
            for i in range(rhodust_rdm_onto_pdm.shape[0]):
                for j in range(rhodust_rdm_onto_pdm.shape[1]):
                    rhodust_rdm_onto_pdm[i,j] = interp([rs_pdm[i,j],ts_pdm[i,j]])
                    count+=1
                    print("Points interpolated: %d/%d"%(count,rhodust_rdm_onto_pdm.size))
        else:
            print("Sorry, I cannot handle 3D geometries yet")
        
        return rhodust_rdm_onto_pdm


    """ Check density maps """
    def check_rhod():
        
        # Interpolate radmc3d rhod onto ProDiMo grid
        rhod_rdm_onto_pdm = interpolate_rhod_rdm_onto_pdm()
        residual_map = (pdm.rhod-rhod_rdm_onto_pdm)
        residual_map_flat = residual_map.flatten()
        residual_map_flat_nonans = residual_map_flat[~np.isnan(residual_map_flat)] 

        # Find stats 
        mu = np.mean(residual_map_flat_nonans)
        std = np.std(residual_map_flat_nonans)
        standard_residual_map = (residual_map-mu)/std
        standard_residual_map_masked = np.ma.masked_invalid(standard_residual_map)
        #print(standard_residual_map_masked.min())
        #print(standard_residual_map_masked.max())
        
        # Choose a point within the radmc3d-interpolated and prodimo domains
        xp = 6.6
        zp = 6.0
        idxp = abs(pdm.x[:,0]-xp).argmin()
        idzp = abs(pdm.z[idxp,:]-zp).argmin()
        rhod_p_pdm = pdm.rhod[idxp,idzp]
        rhod_p_rdm = rhod_rdm_onto_pdm[idxp,idzp]    
        diff = ( abs(rhod_p_rdm-rhod_p_pdm)/rhod_p_pdm )*100

        #print(residual_map_flat_nonans)
        #sys.exit()
        #residual_map_flat_positive = residual_map_flat[residual_map_flat>0]
        #residual_map_flat_negative = residual_map_flat[residual_map_flat<0]

        #print(len(residual_map_flat))
        #print(len(residual_map_flat_positive))
        #print(len(residual_map_flat_negative))

        #sys.exit()
        
        #plt.hist(standard_residual_map.flatten(),density=True,bins=25) # There are positive values, why?!!!!!
        #plt.show()

        #print(np.nanmin(residual_map))   
        #print(np.nanmax(residual_map))
    
        # Plot
        fig,(ax1,ax2,ax3) = plt.subplots(ncols=3,figsize=(12,4))
        xmax = 10 
        zmax = 10

        # ProDiMo
        cs1 = ax1.contourf(pdm.x,pdm.z,np.log10(pdm.rhod))
        cbar1 = fig.colorbar(cs1,ax=ax1)
        ax1.scatter(xp,zp,marker='x',color='white',linewidth=2.5)
        cbar1.set_label("log10 rhod")
        ax1.text(0.05,0.95,'xp = %d'%xp,transform=ax1.transAxes)
        ax1.text(0.05,0.90,'zp = %d'%zp,transform=ax1.transAxes)
        ax1.text(0.05,0.85,'rhod(xp,zp) = %.2e'%rhod_p_pdm,transform=ax1.transAxes)
        ax1.set_xlim(None,xmax)
        ax1.set_ylim(None,zmax)

        ax1.set_title("ProDiMo")

        # radmc3d
        cs2 = ax2.contourf(pdm.x,pdm.z,np.log10(rhod_rdm_onto_pdm))
        cbar2 = fig.colorbar(cs2,ax=ax2)
        ax2.scatter(xp,zp,marker='x',color='white',linewidth=2.5)
        cbar2.set_label("log10 rhod")
        ax2.text(0.05,0.95,'xp = %d'%xp,transform=ax2.transAxes)
        ax2.text(0.05,0.90,'zp = %d'%zp,transform=ax2.transAxes)
        ax2.text(0.05,0.85,'rhod(xp,zp) = %.2e'%rhod_p_rdm,transform=ax2.transAxes)
        ax2.set_xlim(None,xmax)
        ax2.set_ylim(None,zmax)
        ax2.set_title("radmc3d interpolated")

        # Residuals
        vmin = -0.5
        vmax = +0.5
        levels = np.linspace(vmin,vmax,21) 
        cs3 = ax3.contourf(pdm.x,pdm.z,standard_residual_map_masked,
                           vmin=vmin,vmax=vmax,levels=levels)
        ax3.scatter(xp,zp,marker='x',color='white',linewidth=2.5)
        cbar3 = fig.colorbar(cs3,ax=ax3)
        ax3.text(0.05,0.95,'xp = %d'%xp,transform=ax3.transAxes)
        ax3.text(0.05,0.90,'zp = %d'%zp,transform=ax3.transAxes)
        ax3.text(0.05,0.85,'relative diff. = %.1f%%'%diff,transform=ax3.transAxes)
        ax3.set_xlim(None,xmax)
        ax3.set_ylim(None,zmax)

        ax3.set_title("ProDiMo - radmc3d ")

        plt.tight_layout()
        plt.show()
        
        
        return None
        
        
    check_rhod()

    


    return None