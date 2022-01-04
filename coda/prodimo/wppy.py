from coda.mcmax.header import *
from prodimopy.read import calc_surfd

''' Wppy (Wrapper for ProDiMopy) is largely based on the
reading routines availabe in ProDiMoPy. '''

def sdd(save=None):
    ''' Read the model '''
    model=read_prodimo()

    ''' Get the data. sdd containes the surface density
    of the upper hemisphere of the disk. tsdd containes
    the total surface density '''

    sdd=np.reshape(model.sdd[:,0:1],model.sdd.shape[0])
    tsdd=2*sdd
    x=np.reshape(model.x[:,0],model.x.shape[0])

    ''' Save to a file? '''
    if save:
        f=open("wppy_sdd.dat","w")
        for i,j in zip(x,tsdd):
            print(i,j)
            f.write("%.15e %.15e\n"%(i,j))
        f.close()

    ''' Plot '''
    fig,ax=plt.subplots(1,1)
    ax.plot(x,tsdd,'+')
    ax.set(yscale="log")
    plt.show()

    return None

def sigmaa(file=None):
    ''' Read the model '''
    model=read_prodimo()

    ''' Get the data. sdd containes the surface density
    of the upper hemisphere of the disk. tsdd containes
    the total surface density '''

    sigma_full=model.dust.sigmaa
    x=np.reshape(model.x[:,0],model.x.shape[0])
    a=model.dust.asize


    sdd=np.reshape(model.sdd[:,0:1],model.sdd.shape[0])
    tsdd=2*sdd
    x=np.reshape(model.x[:,0],model.x.shape[0])


    ''' Reconstruct density '''
    density_reconstructed=np.zeros(sigma_full.shape[1])
    for j in range(sigma_full.shape[1]):
        density_reconstructed[j]=np.sum(np.reshape(sigma_full[:,j:j+1],sigma_full.shape[0]))

    ''' Plot '''
    fig,(ax1,ax2,ax3)=plt.subplots(1,3,figsize=(25,5))
    ax1.plot(x,tsdd,'+')
    for i in range(sigma_full.shape[0]):
        if i==0:
            ax2.plot(x,sigma_full[i],'--',label="%.3f micron"%(a[i]))
        elif i==sigma_full.shape[0]-1:
            ax2.plot(x,sigma_full[i],'--',label="%.3f micron"%(a[i]))
        else:
            if i%10==0:
                ax2.plot(x,sigma_full[i])
    ax3.plot(x,2*density_reconstructed,label="reconstructed")
    if file:
        data_file=np.loadtxt(file)
        xfile=np.reshape(data_file[:,0:1],data_file.shape[0])
        yfile=np.reshape(data_file[:,1:2],data_file.shape[0])

        ax3.plot(xfile,yfile,label="input file")

    ax1.set(yscale="log")
    ax2.set(yscale="log")
    ax3.set(yscale="log")

    ax2.legend()
    ax3.legend()

    plt.show()

    return None

def td():
    model=read_prodimo()
    x=np.reshape(model.x[:,0],model.x.shape[0])
    td_z0=np.reshape(model.td[:,0],model.td.shape[0])


    ''' Plot '''
    fig,ax=plt.subplots(1,1)
    ax.plot(x,td_z0)
    ax.set(yscale='log')
    plt.show()
    return None

def tg():
    model=read_prodimo()
    x=np.reshape(model.x[:,0],model.x.shape[0])
    tg_z0=np.reshape(model.tg[:,0],model.tg.shape[0])


    ''' Plot '''
    fig,ax=plt.subplots(1,1)
    ax.plot(x,tg_z0)
    ax.set(yscale='log')
    plt.show()
    return None

def A1():
    model=read_prodimo()
    A_list=model.dust.asize
    plt.plot(A_list,'+')
    plt.yscale("log")
    plt.show()
    #print(A_list)
    return A_list

def plot_dust_opac():
    model=read_prodimo()

    print("PLOT: dust opacities ...")
    fig,ax=plt.subplots()
    x=model.dust.lam
    dust=model.dust

    ax.plot(x,dust.kabs,label="absorption")
    ax.plot(x,dust.ksca,label="scattering")
    ax.plot(x,dust.kext,label="extinction")

    ax.set_xlabel(r"wavelength $\mathrm{[\mu m]}$")
    ax.set_ylabel(r"opacity $\mathrm{[cm^2 g(dust)^{-1}]}$")

    ax.set_xlim(np.min(x),np.max(x))
    # ax.set_ylim(1.e-2,None)

    ax.semilogx()
    ax.semilogy()

    #self._dokwargs(ax,**kwargs)
    ax.legend()
    plt.show()
    return None
