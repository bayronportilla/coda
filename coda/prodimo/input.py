from coda.mcmax.header import *

def calc_vmax(r,Mstar,
              theta=None,inc=None):
    r=r*u.au
    Mstar=Mstar*u.Msun

    if theta is not None: theta=(theta*u.deg).to(u.rad)
    else: theta=(90*u.deg).to(u.rad)
        
    if inc is not None: inc=(inc*u.deg).to(u.rad)
    else: inc=(45*u.deg).to(u.rad)
   
    v_r=(c.G*Mstar/r)**0.5*np.sin(theta)*np.sin(inc)

    print(v_r.to(u.km/u.s))
    return None
