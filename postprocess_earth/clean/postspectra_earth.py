import postprocess as pp
import netCDF4 as nc
import sys

if __name__=="__main__":
    color=False
    makemap=False
    if "color" in sys.argv[:]:
        color=True
    if "map" in sys.argv[:]:
        makemap=True
    specs = nc.Dataset("spectra.nc","r")
    pp.makephase(specs,color=color,makemap=makemap)

