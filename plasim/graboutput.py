import numpy as np
import glob
import netCDF4 as nc
import sys

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

if __name__=="__main__":
    keys = []
    xy = []
    name = "output"
    for arg in sys.argv[1:]:
        if arg[0:4]=="key=":
            keys.append(arg.split('=')[1])
        if arg[0:3]=="yx=":
            val = arg.split('=')[1].split(',')
            x = int(val[0])
            y = int(val[1])
            xy.append((x,y))
        if arg[0:5]=="name=":
            name = arg.split('=')[1]
        
    ncoords = len(xy)+1
    nseries = len(keys)
    dfiles = sorted(glob.glob("*.nc"))
    ncd = nc.Dataset(dfiles[0],"r")
    lon = ncd.variables['lon'][:]
    lat = ncd.variables['lat'][:]
    ncd.close()
    npoints = 12*len(dfiles)
    nfiles = len(dfiles)
    output = np.zeros((nseries,ncoords,npoints))
    for n in range(0,nfiles):
        ncd = nc.Dataset(dfiles[n],"r")
        for k in range(0,nseries):
            dvar = ncd.variables[keys[k]][:]
            for jk in range(0,ncoords-1):
                output[k,jk,n*12:(n+1)*12] = dvar[:,xy[jk][0],xy[jk][1]]
            for m in range(0,12):
                output[k,-1,n*12+m] = spatialmath(lat,lon,dvar[m,:,:])
        ncd.close()
        print n
        
    np.save(name+".npy",output)
    
            