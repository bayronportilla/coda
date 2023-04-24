from coda.mcmax.header import *

def wavelength_grid(lam_min,lam_max):

    '''

    Creates a one column file with the value of discrete wavelength points where
    continuum radiative transfer will be performed.

    '''

    Nlam        = 150
    wl_array    = np.logspace(np.log10(lam_min),np.log10(lam_max),Nlam,)

    file = open("wavelength_micron.inp","w")
    file.write("%d\n"%Nlam)
    for i in wl_array:
        print(i)
        file.write(" %.5e \n"%i)

    return None

def stars()
