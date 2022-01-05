from coda.mcmax.header import *
import pandas as pd

def RTconv(file):
    ''' Read file and create dataframe '''
    data=pd.read_csv(file,sep='\s+',usecols=[1,3,5],names=['RT-iteration','J-change','T-change'])

    ''' Ignoring ok=T lines '''
    data=data[data['RT-iteration']!='ok=T']

    ''' Extracting column information '''
    nit=data['RT-iteration'].values
    J=data['J-change'].values
    T=data['T-change'].values
    
    ''' Removing non-desired characters in RT-iteration column''' 
    x=[]
    for i in nit:
        x.append(int(i.replace(":","")))
    x=np.array(x)

    ''' Plotting '''
    fig,((ax1,ax2))=plt.subplots(1,2,figsize=(15,5))
    ax1.plot(x,T,'.')
    ax2.plot(x,J,'.')
    ax1.set(xlabel='Iteration number',ylabel='T-change',yscale='log')
    ax2.set(xlabel='Iteration number',ylabel='J-change',yscale='log')
    plt.show()

    return None

def integrate_line(file):
    ''' A simple solution to find the integrated flux from 
    a two column file '''
    data=np.loadtxt(file)
    v=np.reshape(data[1:,0:1],data.shape[0]-1) # Remove v~1e99 km/s
    f=np.reshape(data[1:,1:2],data.shape[0]-1)

    ''' Integrate using Simpson's rule '''
    print(simps(f,v))
    
    

    return None
