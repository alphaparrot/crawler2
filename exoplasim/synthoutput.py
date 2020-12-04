import numpy as np
import netCDF4 as nc
import sys
import glob
from exoplasim import gcmt

if __name__=="__main__":
    
    prefix = sys.argv[1]
    
    start = int(sys.argv[2])
    try:
        end = int(sys.argv[3])    
    except:
        end = len(glob.glob(prefix+"*.nc"))
    
    keys = ["clt" ,
            "hfls",
            "hfns",
            "hfss",
            "nbr" ,
            "prw" ,
            "rls" ,
            "rlut",
            "rss" ,
            "rst" ,
            "rsut",
            "ssru",
            "stru",
            "sic",
            "sit",
            "ts"  ]
    
    nrecs = end - start + 1
    
    history = {}
    
    for k in keys:
        history[k] = np.zeros(nrecs)
    
    for n in range(start,end+1):
        ncd = nc.Dataset(prefix+".%04d.nc"%n,"r")
        lat = ncd.variables['lat'][:]
        lon = ncd.variables['lon'][:]
        for k in keys:
            try:
                varx = ncd.variables[k][:]
                nt = varx.shape[0]
                x = gcmt.spatialmath(varx,lon=lon,lat=lat)
                history[k][n-start] = x
            except:
                pass
        ncd.close()
        print n
        
    np.save(prefix+"_history.npy",history)
    
    
    