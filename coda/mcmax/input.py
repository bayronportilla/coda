from coda.mcmax.header import *

@dataclass
class File:
    x:float
    y:float
    name:str

    def calculate_mass(self,rlim=None):

        x=(self.x*u.au).to(u.cm).value
        y=self.y

        if rlim is not None:
            rlim=(rlim*u.au).to(u.cm).value
            xlim,ylim=[],[]
            for i,j in zip(x,y):
                if i<rlim:
                    xlim.append(i)
                    ylim.append(j)
            xlim,ylim=np.array(xlim),np.array(ylim)
            dust_mass=(2*np.pi*simps(xlim*ylim,xlim)*u.g).to(u.Msun)
        else:
            dust_mass=(2*np.pi*simps(x*y,x)*u.g).to(u.Msun)

        return dust_mass

    def plot(self):
        x,y=self.x,self.y
        fig,ax=plt.subplots(1,1)
        ax.plot(x,y,'.')
        ax.set(xscale="log",yscale="log")
        plt.show()

    def rescale_mass(self,k,rmin=None):

        '''
        Rescale density profile by a factor of k.

        Parameters:
            k (float): rescaling factor.
            rmin (float): to only rescale the points r>=rmin or r<=rmin.
            If rmin is None, rescale the entire profile.
        '''

        x,y=(self.x*u.au).to(u.cm).value,self.y
        Mold=(2*np.pi*simps(x*y,x)*u.g).to(u.Msun)
        if rmin is None: ynew=k*y
        if rmin is not None:
            ynew=[]
            ynew+=[1*y[i] for i in range (len(self.x)) if self.x[i]<=rmin]
            ynew+=[k*y[i] for i in range (len(self.x)) if self.x[i]>rmin]

        Mnew=(2*np.pi*simps(x*ynew,x)*u.g).to(u.Msun)

        ''' Print info '''
        print("\nInitial mass:",Mold)
        print("Rescaled mass:",Mnew)
        #print("Mnew/Mold: %.3f"%(Mnew/Mold))

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
    ''' '''

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

    ''' Bulk constants. Do not touch them unless the hard-coded quantities in
    MCMax3D have changed. '''
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


def convert_density_file(fobj,psize=None,visual=None,mcmax_like=None):
    '''
    Convert a standard MCMax3D surface density profile
    into a ProDiMo 1D input sdfile. As an argument, it receives
    an object of the class File with the following
    characteristics:

    - fobj.x in au
    - fobj.y in g cm^-2

    The length of the r_array will determine the
    resolution of the ProDiMo model. Choose it wisely!
    '''

    @dataclass
    class GridPointMcmax:
        r:float
        theta:float
        phi:float
        dr:float
        dtheta:float
        dphi:float
        comp:float

        def dV(self):
            # Units: cm^3
            value=self.r**2*np.sin(self.theta)*self.dr*self.dtheta*self.dphi
            return value

        def dz(self):
            # Units: cm
            value=self.r*self.dtheta
            return value

        def zsphere(self):
            # Units: cm
            value=self.r*np.cos(self.theta)
            return value

        def rsphere(self):
            ### Units: cm
            value=self.r*np.sin(self.theta)
            return value

    @dataclass
    class GridPointProdimo:
        r:float
        z:float
        dz:float
        comp:float


    print("Retrieving the exact MCMax3D profile")
    model='/data/users/bportilla/runs_P1/final_runs/recalibration_ppd/run130/'
    #model='/data/users/bportilla/runs_P1/final_runs/recalibration_ppd/run128/'

    ''' Reading particle size array '''
    ai_array=psize[:,0] # microns

    def func_fsize(zoneID):
        ''' Working with zones'''
        ### Import data
        hdu_1=fits.open(model+"output/Zone000%d.fits.gz"%(zoneID))

        ### Composition matrix
        ### C[3]:pmid(rad), C[4]:tmid(rad), C[5]:rmid(au)
        C=hdu_1[6].data                   # g cm^-3

        ### Center coordinates
        rmidAU=hdu_1[0].data[0,0,0,:]
        rmid=(rmidAU*u.au).to(u.cm).value # cm
        tmid=hdu_1[0].data[1,0,:,0]       # rad
        pmid=hdu_1[0].data[2,:,0,0]       # rad

        ### Borders
        rlim=hdu_1[1].data   # cm
        tlim=hdu_1[2].data   # rad
        plim=hdu_1[3].data   # rad

        ### Declare matrix needed for interpolation
        MGP=np.empty((int(len(tmid)*0.5),len(rmid)),dtype='object')

        ### The big loop
        for i in range(C.shape[5]):             # Along r
            dr=rlim[i+1]-rlim[i]
            for j in range(C.shape[4]):         # Along theta
                dt=tlim[j+1]-tlim[j]
                for k in range(C.shape[3]):     # Along phi
                    dp=plim[k+1]-plim[k]
                    comp=C[0,:,0,k,j,i] # g cm^-3
                    GP=GridPointMcmax(rmid[i],tmid[j],pmid[k],dr,dt,dp,comp)
                    if k==0.0 and j<len(tmid)*0.5:
                        MGP[j,i]=GP

        return rmidAU,MGP

    ### Concatenate fsizes and rmidAU
    rmidAU_1,MGP_1=func_fsize(1)[0],func_fsize(1)[1]
    rmidAU_2,MGP_2=func_fsize(2)[0],func_fsize(2)[1]
    MGP=np.concatenate((MGP_1,MGP_2),axis=1)
    r_array=np.concatenate((rmidAU_1,rmidAU_2))

    ''' Interpolate density profile '''
    x_array=fobj.x
    y_array=fobj.y
    cs=CubicSpline(x_array,y_array)
    S_array=np.array([cs(i) for i in r_array])

    ''' Declare gas-to-dust ratio '''
    g2d=100.0
    g2d_array=g2d*np.ones(len(S_array))

    ### Declaring and filling in C matrix for MCMax3D. It only stores
    ### the density.
    C_mcmax=np.zeros((psize.shape[0],MGP.shape[0],MGP.shape[1]))
    for k in range(C_mcmax.shape[0]):           # Over grain size
        for i in range(C_mcmax.shape[1]):       # Over theta
            for j in range(C_mcmax.shape[2]):   # Over radius
                C_mcmax[k,i,j]=MGP[i,j].comp[k]

    ### Read ProDiMo grid
    model=read_prodimo()
    r_array_p=np.reshape(model.x[:,0:1],model.x.shape[0])   # au
    z_matrix=model.z                                        # au, dim:len(r)*len(z)

    ### Declaring and filling in C matrix for ProDiMo. Composition is
    ### initialised to None.
    C_prodimo=np.empty((psize.shape[0],z_matrix.shape[1],z_matrix.shape[0]),
                        dtype='object')
    for k in range(C_prodimo.shape[0]):             # Over grain size
        for i in range(C_prodimo.shape[1]):         # Over z
            for j in range(C_prodimo.shape[2]):     # Over radius
                if i<C_prodimo.shape[1]-1:
                    dz=((z_matrix[j,i+1]-z_matrix[j,i])*u.au).to(u.cm).value
                else:
                    dz=((z_matrix[j,i]-z_matrix[j,i-1])*u.au).to(u.cm).value
                GP_prodimo=GridPointProdimo(r_array_p[j],z_matrix[j,i],dz,None)
                C_prodimo[k,i,j]=GP_prodimo

    ''' Start 2D interpolation '''
    ### Sampling MCMax3D info into 1D arrays.
    for k in range(psize.shape[0]):
        rsph_array=[]
        zsph_array=[]
        csph_array=[]
        for j in range(MGP.shape[1]):                           # Over radius
            for i in range(MGP.shape[0]):                       # Over theta
                rsph_array.append((MGP[i,j].rsphere()*u.cm).to(u.au).value)
                zsph_array.append((MGP[i,j].zsphere()*u.cm).to(u.au).value)
                csph_array.append(np.log10(C_mcmax[k,i,j]))     # Note we store the
                                                                # log10 of the density

        ### Uncomment the preferred interpolator. Rbf seems to work nicely.
        #f=interp2d(rsph_array,zsph_array,csph_array,kind='linear',bounds_error=False)
        f=Rbf(rsph_array,zsph_array,csph_array,function='linear')

        ### Interpolate densities into ProDiMo grid
        for i in range(C_prodimo.shape[1]):                 # Over z
            for j in range(C_prodimo.shape[2]):             # Over radius
                rprodimo=C_prodimo[k,i,j].r                 # au
                zprodimo=C_prodimo[k,i,j].z                 # au
                C_prodimo[k,i,j].comp=f(rprodimo,zprodimo)  # g cm^-3

    ''' Declaring fsize matrix '''
    """
    fsize=np.zeros((C_prodimo.shape[0],C_prodimo.shape[2]))
    for k in range(C_prodimo.shape[0]):          # Over composition
        for j in range(C_prodimo.shape[2]):      # Over radius
            rho_vals,z_vals=[],[]
            for i in range(C_prodimo.shape[1]):  # Over height
                rho_vals.append(10**float(C_prodimo[k,i,j].comp))
                z_vals.append((C_prodimo[k,i,j].z*u.au).to(u.cm).value)
            fsize[k,j]=2*simps(rho_vals,z_vals)
    """
    fsize=np.zeros((C_prodimo.shape[0],C_prodimo.shape[2]))
    for k in range(C_prodimo.shape[0]):          # Over composition
        for j in range(C_prodimo.shape[2]):      # Over radius
            S_val=0.0
            for i in range(C_prodimo.shape[1]):  # Over height
                rho_val=10**float(C_prodimo[k,i,j].comp)
                S_val+=rho_val*C_prodimo[k,i,j].dz
            fsize[k,j]=2*S_val

    if visual:
        ''' Quick plot to check things out '''
        fig,(ax1,ax2)=plt.subplots(1,2,figsize=(15,5))
        for i in range(psize.shape[0]):
            if i==0:
                ax1.plot(r_array,fsize[i],'--',label="%.3f micron"%(psize[i,0]))
            elif i==psize.shape[0]-1:
                ax1.plot(r_array,fsize[i],'--',label="%.3f micron"%(psize[i,0]))
            else:
                if i%10==0:
                    ax1.plot(r_array,fsize[i])

        density_reconstructed=np.zeros(fsize.shape[1])
        for j in range(fsize.shape[1]):
            density_reconstructed[j]=np.sum(np.reshape(fsize[:,j:j+1],fsize.shape[0]))

        ax2.plot(fobj.x,fobj.y,label='original')
        ax2.plot(r_array,density_reconstructed,'.',label='reconstructed')

        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax1.legend()
        ax2.legend()
        plt.show()

    ''' Convert to prodimopy units '''
    ai_array=(ai_array*u.micron).to(u.cm)
    r_array=(r_array*u.au).to(u.cm)


    ''' Calling prodimopy '''
    write("sdprofile.in",r_array.value,S_array*g2d,g2d_array,ai_array.value,fsize)

    return None


def find_mass_temp(fobj,psize,cprodimo):
    ### Computing fsize matrix
    fsize=np.zeros((cprodimo.shape[0],cprodimo.shape[2]))
    r_array=[]
    for k in range(cprodimo.shape[0]): # Over composition
        for j in range(cprodimo.shape[2]): # Over radius
            if k==0:
                r_array.append(cprodimo[0,0,j].r)
            rho_vals,z_vals=[],[]
            for i in range(cprodimo.shape[1]): # Over height
                rho_vals.append(10**float(cprodimo[k,i,j].comp))
                z_vals.append((cprodimo[k,i,j].z*u.au).to(u.cm).value)
            fsize[k,j]=2*simps(rho_vals,z_vals)

    density_reconstructed=np.zeros(fsize.shape[1])
    for j in range(fsize.shape[1]):
        density_reconstructed[j]=np.sum(np.reshape(fsize[:,j:j+1],fsize.shape[0]))

    ''' Quick plot to check things out '''
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(15,5))
    for i in range(psize.shape[0]):
        if i==0:
            ax1.plot(r_array,fsize[i],'--',label="%.3f micron"%(psize[i,0]))
        elif i==psize.shape[0]-1:
            ax1.plot(r_array,fsize[i],'--',label="%.3f micron"%(psize[i,0]))
        else:
            if i%10==0:
                ax1.plot(r_array,fsize[i])
    ax2.plot(fobj.x,fobj.y,label='original')
    ax2.plot(r_array,density_reconstructed,'.',label='reconstructed')
    #plt.plot(fobj.x,fobj.y,label='original')
    #plt.plot(r_array,density_reconstructed,'.',label='reconstructed')

    #ax1.set_yscale("log")
    #plt.yscale("log")
    #ax1.legend()
    ax1.set_yscale("log")
    ax2.set_yscale("log")
    ax1.legend()
    ax2.legend()
    plt.show()
    #plt.show()




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