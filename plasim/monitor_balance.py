import numpy as np
import os
import netCDF4 as nc
import glob


#This version lets the model relax.

def spatialmath(lt,ln,variable,mean=True,radius=6.371e6):
    lt1 = np.zeros(len(lt)+1)
    lt1[0] = 90
    for n in range(0,len(lt)-1):
        lt1[n+1] = 0.5*(lt[n]+lt[n+1])
    lt1[-1] = -90
    ln1 = np.zeros(len(ln)+1)
    ln1[0] = -2.8125
    for n in range(0,len(ln)-1):
        ln1[n+1] = 0.5*(ln[n]+ln[n+1])
    ln1[-1] = 360.0-2.8125
    
    lt1*=np.pi/180.0
    ln1*=np.pi/180.0
    
    darea = np.zeros((len(lt),len(ln)))
    for jlat in range(0,len(lt)):
        for jlon in range(0,len(ln)):
            dln = ln1[jlon+1]-ln1[jlon]
            darea[jlat,jlon] = (np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*dln
    
    svar = variable*darea
    if mean:
        outvar = np.sum(svar)/np.sum(darea)
    else:
        outvar = np.sum(svar) * radius**2
    
    return outvar



def energybalance():
    files = sorted(glob.glob(".nc"))
    sbalance = np.zeros(len(files))
    toabalance = np.zeros(len(files))
    for n in range(0,len(files)):
        ncd = nc.Dataset(files[n],"r")
        ntr = ncd.variables['ntr'][:]
        hfns = ncd.variables['hfns'][:]
        ncd.close()
        topt = np.zeros(12)
        bott = np.zeros(12)
        for m in range(0,12):
            topt[m] = spatialmath(ntr[m,:,:])
            bott[m] = spatialmath(hfns[m,:,:])
        sbalance[n] = np.mean(bott)
        toabalance[n] = np.mean(topt)
    return (sbalance,toabalance)
        
    
        

if __name__=="__main__":
  os.system("echo 'SURFACE      TOA'>balance.log")
  surf,toa = energybalance()
  for n in range(0,len(surf)):
      os.system("echo '%02.6f  %02.6f'>>balance.log"%(surf[n],toa[n]))
      
    
