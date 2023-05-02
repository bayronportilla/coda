from coda.mcmax.header import *
import radmc3dPy
from radmc3dPy import analyze

def temperature_midplane():

    """
    It supports only averaged quantities

    """

    d           = analyze.readData()
    d.readDustTemp()
    r   = (d.grid.x*u.cm).to(u.au).value



    # Average over azimuth
    Td_small  = np.sum(d.dusttemp[:,-1,:,0],axis=1)/d.dusttemp[:,-1,:,0].shape[1]
    Td_large  = np.sum(d.dusttemp[:,-1,:,1],axis=1)/d.dusttemp[:,-1,:,1].shape[1]

    #print(Td_small)
    #print(Td_large)

    #sys.exit()
    #plt.plot(d.grid.x,d.dusttemp[:,-1,0,0])
    #plt.show()
    #print(d.dusttemp.shape)
    #print(d.grid.y)

    print(r.shape)
    print(Td_small.shape)
    print(Td_large.shape)


    file=open("temperature_profile_averaged.dat","w")
    file.write("# r(cm)   Td_small (K)    Td_large (K)\n")
    for i in range(len(r)):
        file.write("%.15e   %.15e   %.15e\n"
            %(r[i],Td_small[i],Td_large[i]))
    file.close()
    return None


def Sigma_dust():

    """


    """

    d   = analyze.readData()
    d.getSigmaDust()
    r   = (d.grid.x*u.cm).to(u.au).value

    #print(d.sigmadust)
    #sys.exit()
    plt.plot(r,d.sigmadust)
    plt.show()
    #print(d.dusttemp.shape)
    #print(d.grid.y)

    """
    print(d.dusttemp[0,-1,:,0])
    sys.exit()
    file=open("temperature_profile.dat","w")
    file.write("r(cm)   T")

    for i in range(len(d.grid.x)):
        file.write("%.15e   %.15e   %.15e\n"
            %(d.grid.x[i],d.dusttemp[i,-1,0,0],d.dusttemp[i,-1,0,1]))
    file.close()
    """

    return None
