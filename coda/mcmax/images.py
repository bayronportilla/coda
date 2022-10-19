from coda.mcmax.header import *
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization import astropy_mpl_style
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle
plt.style.use(astropy_mpl_style)

"""

Assumes a fits image (ProDiMo output) is an instance of a class.

Works only on 2D images with 1 HDU

Example
-------
>>> importlib.reload(coda.mcmax.images)
>>> im=images.Image('../run064/LINE_3D_005.fits.rot.conv.line.mom0.Tb.fits')
>>> im.visualize()
"""

@dataclass
class Cube:

    name:str
    #data:float

    def visualize(self):

        # Print structure of the file
        fits.info(self.name)

        hdu     = fits.open(self.name)
        data    = hdu[0].data
        wcs     = WCS(hdu[0].header)


        c = SkyCoord('04h32m00.00s', '+24d00m00s', frame='icrs')

        d           = 113.47*u.pc   # Must be in pc
        semiaxis    = 94           # Must be in au
        inc         = 51.7*u.deg    # Inclination

        ecc         = np.sin( inc.to(u.rad) )
        a           = (semiaxis/d.value)*u.arcsec # Semimajor axis in arcsec
        b           = a*(1-ecc**2)**0.5

        adeg = a.to(u.deg).value
        bdeg = b.to(u.deg).value

        print(a,b)

        print(c)
        print(c.ra.deg)
        fig=plt.figure()
        ax=plt.subplot(projection=wcs,slices=('x','y',0))
        ax.imshow(data[0,:,:])
        ax.scatter(c.ra.deg,c.dec.deg,transform=ax.get_transform('icrs'))

        e = Ellipse((c.ra.deg,c.dec.deg),2*adeg,2*bdeg,-70.4,
            edgecolor='black', facecolor='none', linestyle=':',
            linewidth=5,alpha=0.2,
            transform=ax.get_transform('icrs'))



        ax.add_patch(e)

        #ax.plot(c)
        #ax.plot(c)
        plt.show()

        return data
