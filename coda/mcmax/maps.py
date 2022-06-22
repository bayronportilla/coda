from coda.mcmax.header import *
#plt.style.use('fancy')

def axes1D(x,y,xlog=None,ylog=None,xlabel=None,
           ylabel=None):

    ''' Default plotting '''
    fsize=14
    fig,ax=plt.subplots()
    ax.plot(x,y)

    ax.set_xlabel(r"x",fontsize=fsize)
    ax.set_ylabel(r"y",fontsize=fsize)

    ax.set_xscale("log")
    ax.set_yscale("log")

    ''' Modifiers '''
    if xlog is False:
        ax.set_xscale("linear")

    if ylog is False:
        ax.set_yscale("linear")

    if xlabel is not None:
        ax.set_xlabel(r"%s"%xlabel,fontsize=fsize)

    if ylabel is not None:
        ax.set_ylabel(r"%s"%ylabel,fontsize=fsize)


    """
    ''' Cosmetics '''
    xmin,xmax=ax.get_xlim()
    ymin,ymax=ax.get_ylim()
    xy=(xmin,ymin)
    width=xmax-xmin
    height=ymax-ymin
    p=patches.Rectangle(xy,width,height,hatch='///',fill=False, zorder=-10)
    ax.add_patch(p)
    plt.tight_layout()

    ''' Save figure? '''
    if save:
        fig.savefig("map.ps",format='ps')
    """
    return ax


@dataclass
class Map:

    '''
    Base class that contains the reading and plotting routines of MCMax3D output
    fields.

    Accepted fieldnames are: 'temp',...

    Examples
    --------
    * To read and p lot the equilibrium temperature in the midplane do:

        1. Create the map. This will call mcmax3dpy
        >>> td=maps.create_maps("<path-to-model>",fieldname='temp')

        2. Compute the azimuthally averaged grain temperature in the midplane,
        visualize and save to a file 'field_ave.dat'. Note that second column
        is log10(td).
        >>> td.plot_midplane(ave=True)

    '''

    name:str
    x:float
    y:float
    z:float
    r:float
    phi:float
    f:float

    def plot_vertical(self,clines=None,fname=None,log=True,
                      ylim=None,save=None,ave=None,cmap=None):
        ''' Read coordinates and fields '''
        x=self.x[0,0:int(self.x.shape[1]*0.5),:]
        y=self.y[0,0:int(self.y.shape[1]*0.5),:]
        z=self.z[0,0:int(self.z.shape[1]*0.5),:]
        f=self.f[0,0:int(self.f.shape[1]*0.5),:]
        f=np.ma.masked_invalid(f)
        if ave:
            f_full=self.f[:,0:int(self.f.shape[1]*0.5),:]
            for i in range(1,f_full.shape[0]):
                f+=f_full[i,:,:]
            f=f/f_full.shape[0]
            f=np.ma.masked_invalid(f)


        ''' Derived properties '''
        rho=(x**2+y**2)**0.5
        rho0=np.ones((rho.shape[0],rho.shape[1]))*rho[-1]
        z0=np.flip(z[:,-1])
        z00=np.ones((rho0.shape[0],rho0.shape[1]))
        for i in range(0,z00.shape[0]):
            z00[i,:]=z0[i]

        ''' Default plotting '''
        fsize=14
        fig,ax=plt.subplots()
        lmin,lmax=f.min(),f.max()
        #lmin,lmax=-30,f.max()
        levels=np.linspace(lmin,lmax,100)
        if not cmap:
            CS=ax.contourf(rho,z/rho,f,levels=levels,extend='neither')
        if cmap:
            CS=ax.contourf(rho,z/rho,f,levels=levels,extend='neither',cmap=cmap)
        CL=ax.contour(rho,z/rho,f,levels=clines,colors="white",linestyles="dashed",
                      linewidths=1.0)
        CB=fig.colorbar(CS,format="%.1f")
        CB.set_label(fname,fontsize=fsize)
        CB.set_ticks(np.linspace(lmin,lmax,5))
        ax.set(xscale="log",xlim=(rho[-1].min(),rho[-1].max()),ylim=(0.0,0.3))
        ax.set_xlabel(r"$r$ (AU)",fontsize=fsize)
        ax.set_ylabel(r"$z/r$",fontsize=fsize)

        ''' Modifiers '''
        if not log:
            ax.set_xscale("linear")
        if ylim:
            ax.set_ylim(ylim)
        if clines:
            ax.clabel(CL,clines,inline=True,fontsize=14,colors="white",manual=True,
            fmt='%1.1f')

        ''' Cosmetics '''
        xmin,xmax=ax.get_xlim()
        ymin,ymax=ax.get_ylim()
        xy=(xmin,ymin)
        width=xmax-xmin
        height=ymax-ymin
        p=patches.Rectangle(xy,width,height,hatch='///',fill=False, zorder=-10)
        ax.add_patch(p)
        plt.tight_layout()

        ''' Save figure? '''
        if save:
            fig.savefig("map.ps",format='ps')

        plt.show()

        return None

    ''' Plot the values on the midplane '''
    def plot_midplane(self,clines=None,fname=None,log=True,ave=None,fullave=None,
                      xlog=None,ylog=None,xlabel=None,ylabel=None):

        """
        At the moment, this routine only supports plotting averaged quantities
        """

        def plot_statistics(matrix,array,F,visual=None):
            median,mad=np.median(array),median_absolute_deviation(array)
            if visual is True:
                fig,((ax1,ax2))=plt.subplots(1,2,figsize=(15,5))
                ax1.imshow(matrix)
                ax1.set(xlabel="radial",ylabel="azimuthal")
                ax2.plot(array,".")
                ax2.axhline(median)
                ax2.axhspan(median-F*mad,median+F*mad,alpha=0.2)
                ax2.text(0,median,"Median=%.6f"%median)
                ax2.text(0,median-F*mad,"$-%d\mathrm{MAD}$"%F)
                ax2.text(0,median+F*mad,"$+%d\mathrm{MAD}$"%F)
                ax2.set(xlabel="radial",ylabel="standard derivation")
                return fig
            else:
                return median,mad

        x=self.x[:,int(self.x.shape[1]*0.5),:]
        y=self.y[:,int(self.y.shape[1]*0.5),:]
        z=self.z[:,int(self.z.shape[1]*0.5),:]
        r=self.r[:,int(self.r.shape[1]*0.5),:]
        phi=self.phi[:,int(self.phi.shape[1]*0.5),:]
        f=(self.f[:,int(self.f.shape[1]*0.5),:]+self.f[:,int(self.f.shape[1]*0.5+1),:])*0.5
        r,phi=np.reshape(r[0],r.shape[1]),np.reshape(phi[:,0:1],phi.shape[0])

        if ave:

            std_array=[np.std(np.reshape(f[:,i:i+1],f.shape[0])) for i in range(0,f.shape[1])]

            while True:
                F=int(input("Enter the F-value: "))
                fig=plot_statistics(f,std_array,F,visual=True)
                fig.show()
                useF=bool(int(input("Use this F? (1 means yes, 0 means No): ")))
                if useF==True:
                    break
            median=plot_statistics(f,std_array,F,visual=False)[0]
            mad=plot_statistics(f,std_array,F,visual=False)[1]

            f_ave,r_ave=[],[]
            for i in range(0,f.shape[1]):
                if abs(std_array[i]-median)<F*mad and std_array[i]!=0:
                    f_ave.append(np.sum(np.reshape(f[:,i:i+1],f.shape[0]))/f.shape[0])
                    r_ave.append(r[i])
            f_ave=np.array(f_ave)
            r_ave=np.array(r_ave)

            file=open("field_ave.dat","w")
            for i,j in zip(r_ave,f_ave):
                file.write("%.15e %.15e\n"%(i,j))
            file.close()

        if fullave:
            std_array=[np.std(np.reshape(f[:,i:i+1],f.shape[0])) for i in range(0,f.shape[1])]
            fig=plot_statistics(f,std_array,0.0,visual=True)
            fig.show()
            f_ave,r_ave=[],[]
            for i in range(0,f.shape[1]):
                f_ave.append(np.sum(np.reshape(f[:,i:i+1],f.shape[0]))/f.shape[0])
                r_ave.append(r[i])
            f_ave,r_ave=np.array(f_ave),np.array(r_ave)
            file=open("field_ave.dat","w")
            for i,j in zip(r_ave,f_ave):
                file.write("%.15e %.15e\n"%(i,j))
            file.close()

        ''' Plotting '''
        ax=axes1D(r_ave,10**f_ave,xlog,ylog,xlabel,ylabel)
        plt.show()

        return None


def create_maps(model,fieldname,cpd=None):

    """
    Stores the three-dimensional map of every supported MCMax3D quantity
    to be plotted or analyzed with other routines.

    Let op!
    * The grid resolution in each zone of the PPD must be the same.
    * Only one CPD zone is supported.
    """

    if fieldname!="amean":
        zones=mread.read_zones(model+"/output/")
        if cpd is None:
            x=np.empty((int(zones[0].np),int(zones[0].nt),int(zones[0].nr)),dtype=float)
            y,z,r,phi,f=x,x,x,x,x
            for zone in zones:
                if zone.x0==0.0 and zone.y0==0.0:
                    x=np.append(x,zone.x,axis=2)
                    y=np.append(y,zone.y,axis=2)
                    z=np.append(z,zone.z,axis=2)
                    r=np.append(r,zone.r,axis=2)
                    phi=np.append(phi,zone.phi,axis=2)
                    field=np.log10(getattr(zone,fieldname))
                    f=np.append(f,field,axis=2)
            x,y,z=x[:,:,int(zones[0].nr):],y[:,:,int(zones[0].nr):],z[:,:,int(zones[0].nr):]
            r,phi=r[:,:,int(zones[0].nr):],phi[:,:,int(zones[0].nr):]
            f=f[:,:,int(zones[0].nr):]
        if cpd is True:
            for zone in zones:
                if zone.x0!=0.0 or zone.y0!=0.0:
                    field=np.log10(getattr(zone,fieldname))
                    x=zone.x
                    y=zone.y
                    z=zone.z
                    r=zone.r
                    phi=zone.phi
                    f=field

    if fieldname is "amean":
        zoneID=int(input("Introduce the zoneID: "))
        zones=mread.read_zones(model+"/output/")
        zone=zones[zoneID-1]
        x=np.empty((int(zones[0].np),int(zones[0].nt),int(zones[0].nr)),dtype=float)
        y,z,r,phi,f=x,x,x,x,x
        psize_mcmax=mread.read_particle_sizes_modified("./Particles",zoneID)
        print("min psize in zone %d = %.2f micron"%(zoneID,psize_mcmax.min()))
        print("max psize in zone %d = %.2f micron"%(zoneID,psize_mcmax.max()))
        zone.calc_amean_modified(psize_mcmax,zoneID)
        x=np.append(x,zone.x,axis=2)
        y=np.append(y,zone.y,axis=2)
        z=np.append(z,zone.z,axis=2)
        r=np.append(r,zone.r,axis=2)
        phi=np.append(phi,zone.phi,axis=2)
        f=np.append(f,zone.amean,axis=2)
        x,y,z=x[:,:,int(zone.nr):],y[:,:,int(zone.nr):],z[:,:,int(zone.nr):]
        r,phi=r[:,:,int(zone.nr):],phi[:,:,int(zone.nr):]
        f=f[:,:,int(zone.nr):]
        print("fmin %.3f, fmax %.3f"%(f.min(),f.max()))




    return Map(fieldname,x,y,z,r,phi,f)



def tau(model):
    '''
    Extract information about the optical depth
    for a given zone
    '''

    ''' Loading data '''
    data=np.loadtxt(model+"/output/heightR0003.dat")

    rgrid=np.loadtxt(model+"/output/radgrid0003.dat")
    tgrid=np.loadtxt(model+"/output/thetagrid0003.dat")
    pgrid=np.loadtxt(model+"/output/phigrid0003.dat")

    rgridMP=[(rgrid[i]+rgrid[i+1])*0.5 for i in range(len(rgrid)-1)]
    tgridMP=[(tgrid[i]+tgrid[i+1])*0.5 for i in range(len(tgrid)-1)]
    pgridMP=[(pgrid[i]+pgrid[i+1])*0.5 for i in range(len(pgrid)-1)]

    '''
    zones=mread.read_zones(model+"/output")
    r=zones[2].r[0,0:][0]
    t=zones[2].theta[0,:,0]
    p=np.reshape(zones[2].phi[:,:,0][:,0:1],zones[2].phi.shape[0])
    '''

    tau_r=np.reshape(data[:,0:1],data.shape[0])
    tau_h=np.reshape(data[:,1:2],data.shape[0])


    tau_x=[tau_r[i]*np.sin(tgridMP[i]*np.cos(pgridMP[0])) for i in range(len(tgridMP))]


    print(tau_x)
    plt.plot(tau_x,tau_h,'.')
    plt.xlabel("x (AU)")
    plt.ylabel("h (AU)")
    plt.show()

    return None
