from coda.mcmax.header import *
import matplotlib.ticker

def read_opacities(model,Nzones,fvSI,fvC):

    ############################################################
    # Initializing stuff
    Nbins=np.zeros(Nzones)
    apows=np.zeros(Nzones)
    psizes_min=np.zeros(Nzones)
    psizes_max=np.zeros(Nzones)
    psizes=[]

    ############################################################
    # Working dir
    folder=model+'/Particles'
    path_to_input=model+'/input.dat'
    infile=open(path_to_input).readlines()
    for i in range(1,Nzones+1):
        for line in infile:
            if line.split("=")[0]==("computepart0%d:ngrains"%(i)):
                Nbins[i-1]=int(line.split("=")[1])
            if line.split("=")[0]==("computepart0%d:amin"%(i)):
                psizes_min[i-1]=float(line.split("=")[1])
            if line.split("=")[0]==("computepart0%d:amax"%(i)):
                psizes_max[i-1]=float(line.split("=")[1])
            if line.split("=")[0]==("computepart0%d:apow"%(i)):
                apows[i-1]=float(line.split("=")[1])

    for i in range(0,Nzones):
        psizes.append((psizes_min[i],psizes_max[i]))


    ############################################################
    # Creating lists for both cases
    case2=[]
    for filename in os.listdir(folder):
        if fnmatch.fnmatch(filename,('*_f%.2f_f%.2f.fits.gz'%(fvSI,fvC))):
            case2.append(filename)


    ############################################################
    # Define class 'Archivo'
    class Archivo:
        def __init__(self,filename):
            self.filename=filename

        def zone(self):
            return self.filename[8:12]

        def binnum(self):
            return self.filename[13:17]

        def f(self):
            hdulist=fits.open(model+"/Particles/%s"%(self.filename))
            hdu=hdulist[0]
            hdr=hdu.header
            amin_bin=hdr["R_MIN"]
            amax_bin=hdr["R_MAX"]
            apow=hdr["R_POW"]
            z=int(self.zone())
            bins=Nbins[z-1]
            amin=psizes[z-1][0]
            amax=psizes[z-1][1]
            f_num=amax_bin**(-apow+4)-amin_bin**(-apow+4)
            f_den=amax**(-apow+4)-amin**(-apow+4)
            value=f_num/f_den
            return value

        def getEntry(self,i,j):
            hdulist=fits.open(model+"/Particles/%s"%(self.filename))
            data=np.transpose(hdulist[0].data)
            value=data[i][j]
            return value

        def getWavelength(self):
            hdulist=fits.open(model+"/Particles/%s"%(self.filename))
            data=np.transpose(hdulist[0].data)
            value=np.reshape(data[:,0:1],data.shape[0])
            return value


    ############################################################
    # Create wavelength array
    p2=Archivo(case2[0])
    wl2=p2.getWavelength()


    ############################################################
    # Create matrices to store values of Ai,Bi and Ci for case 1
    # and case 2

    ext_matrix=np.zeros((Nzones,len(wl2),int(max(Nbins))))
    abso_matrix=np.zeros((Nzones,len(wl2),int(max(Nbins))))
    sca_matrix=np.zeros((Nzones,len(wl2),int(max(Nbins))))


    ############################################################
    # Filling matrices
    Nbin_index=np.zeros(Nzones)
    for particle in case2:
        p=Archivo(particle)
        zone_index=int(p.zone())
        j=int(Nbin_index[zone_index-1])
        for i in range(0,len(wl2)):
            ext_matrix[zone_index-1][i][j]=p.f()*p.getEntry(i,1)
            abso_matrix[zone_index-1][i][j]=p.f()*p.getEntry(i,2)
            sca_matrix[zone_index-1][i][j]=p.f()*p.getEntry(i,3)
        Nbin_index[zone_index-1]+=1


    ############################################################
    # Adding entries
    sum_ext=np.zeros((Nzones,len(wl2)))
    sum_abso=np.zeros((Nzones,len(wl2)))
    sum_sca=np.zeros((Nzones,len(wl2)))

    for j in range(0,Nzones):
        for i in range(0,len(wl2)):
            sum_ext[j][i]=np.sum(ext_matrix[j][i])
            sum_abso[j][i]=np.sum(abso_matrix[j][i])
            sum_sca[j][i]=np.sum(sca_matrix[j][i])


    ############################################################
    # Creating HDU's for case 2
    hdu_ext=np.zeros((Nzones+1,len(wl2)))
    hdu_abso=np.zeros((Nzones+1,len(wl2)))
    hdu_sca=np.zeros((Nzones+1,len(wl2)))

    hdu_ext[0]=wl2
    hdu_abso[0]=wl2
    hdu_sca[0]=wl2

    for j in range(0,Nzones):
        hdu_ext[j+1]=sum_ext[j]
        hdu_abso[j+1]=sum_abso[j]
        hdu_sca[j+1]=sum_sca[j]

    print(np.transpose(hdu_ext).shape)
    print(np.transpose(hdu_ext))

    hdu_ext=fits.PrimaryHDU(np.transpose(hdu_ext))
    hdu_abso=fits.PrimaryHDU(np.transpose(hdu_abso))
    hdu_sca=fits.PrimaryHDU(np.transpose(hdu_sca))

    hdu_ext.writeto(model+'/ext.fits',overwrite=True)
    hdu_abso.writeto(model+'/abs.fits',overwrite=True)
    hdu_sca.writeto(model+'/sca.fits',overwrite=True)

    return None

def plot_opacities(model,label=None,save=None):
    hdulist_ext=fits.open(model+'/ext.fits') # row: wl, col: z1,z2,z3
    hdulist_abso=fits.open(model+'/abs.fits')
    hdulist_sca=fits.open(model+'/sca.fits')

    ext=hdulist_ext[0].data
    abso=hdulist_abso[0].data
    sca=hdulist_sca[0].data

    Nzones=ext.shape[1]-1
    print("Nzones=%d"%Nzones)

    #fig=plt.figure(figsize=(13,4))
    fig=plt.figure(figsize=(8,3))
    gs=gridspec.GridSpec(1,Nzones,hspace=0.0)
    lw=2.0
    fsize=12
    for i in range(0,Nzones):
        ax=plt.subplot(gs[0,i])
        ax.plot(ext[:,0:1],ext[:,i+1:i+2],label=r"extinction",color="cornflowerblue",linewidth=lw)
        ax.plot(abso[:,0:1],abso[:,i+1:i+2],label=r"absorption",color="salmon",linewidth=lw)
        ax.plot(sca[:,0:1],sca[:,i+1:i+2],label=r"scattering",color="seagreen",linewidth=lw)
        ax.set(xscale='log',yscale='log')
        ax.set_xlabel(r"$\lambda \, (\mu \mathrm{m})$",fontsize=fsize)
        ax.margins(0.02)
        ax.tick_params(labelsize=fsize)
        if i==0:
            ax.set_ylabel(r"Dust opacities (cm$^2$/g(dust))",fontsize=fsize)
            ax.legend(fontsize=12,loc="lower left",frameon=False)
            ''' Following lines were used to create the figure in Portilla-Revelo et al. 2021 '''
            ymajor=matplotlib.ticker.LogLocator(base=1000,subs=[1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900])
            ax.yaxis.set_major_locator(ymajor)
        if label is not None:
            ax.annotate("%s"%label[i],(0.85,0.9),xycoords='axes fraction',ha='center',va='center',color="dimgray")
    gs.tight_layout(fig,w_pad=2)
    plt.show()
    ''' Save figure? '''
    if save:
        fig.savefig("opacities.ps",format='ps')
    return None



def psize_bin(model,fvSI,fvC):
    ''' Read names of particle files '''
    list_names=[]
    for filename in os.listdir(model+"/Particles"):
        if fnmatch.fnmatch(filename,('*_f%.2f_f%.2f.fits.gz'%(fvSI,fvC))):
            list_names.append(filename)


    ''' Reading header keyword '''
    R_MIN_list=[]
    R_MAX_list=[]
    A_list=[]
    for name in list_names:
        hdulist=fits.open(model+"/Particles/"+name)
        hdu=hdulist[0]
        R_MIN_list.append(hdu.header['R_MIN'])
        R_MAX_list.append(hdu.header['R_MAX'])
        A_list.append(hdu.header['A1'])


    ''' Sorting values '''
    A_list=np.sort(A_list)

    return A_list
