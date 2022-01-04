from coda.mcmax.header import *
from matplotlib.colors import LogNorm
''' 
This library is largely based on the Gofish package by 
Richard Teague 
'''

def get_profile(image,padir,
                visual=None,write=None,compara=None,fov=None,
                bunit=None,pxscale=None):

    cube=imagecube(image,FOV=fov,bunit=bunit,pixel_scale=pxscale)
    widir=20.0
    PAmin=padir-0.5*widir
    PAmax=padir+0.5*widir
    xm, ym, dym = cube.radial_profile(inc=51.7,PA=158.6,dist=113.43,
                                      x0=0.0,y0=0.0,assume_correlated=False,
                                      PA_min=PAmin,PA_max=PAmax,dr=0.025,
                                      unit='Jy/beam')
    xm=xm*113.43

    if visual:
        vmax=np.percentile(cube.data,99)
        vmin=np.percentile(cube.data,1)
        fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5))
        ax1.imshow(cube.data,origin='lower',extent=cube.extent,vmin=vmin,vmax=vmax)
        cube.plot_mask(ax=ax1,r_min=0.02,r_max=3.0,PA_min=PAmin,PA_max=PAmax,inc=51.7,PA=158.6,mask_frame='disk')
        ax2.errorbar(xm,ym,dym)
        ax2.set(xlabel="r (AU)",ylabel="flux density")
        plt.show()
    if write:
        file=open("extracted_signal.dat","w")
        for i,j,k in zip(xm,ym,dym):
            file.write("%.15e %.15e %.15e\n"%(i,j,k))
        file.close()
    if compara:
        data_keppler=np.loadtxt("../cut_along_c_obs.dat")
        x=np.reshape(data_keppler[:,0:1],data_keppler.shape[0])
        y=np.reshape(data_keppler[:,1:2],data_keppler.shape[0])
        dy=np.reshape(data_keppler[:,2:3],data_keppler.shape[0])
        
        fig,ax=plt.subplots()
        ax.errorbar(x,y,dy,label="keppler")
        ax.errorbar(x,y*0.75,dy*0.75,label="keppler x 0.75")
        ax.errorbar(xm,ym*1000,dym*1000,label="isella")
        ax.set(xlabel="r (AU)",ylabel="flux density (mJy/beam)")
        ax.legend()
        plt.show()
        
    return None

#get_profile(visual=True,write=False,compara=False)