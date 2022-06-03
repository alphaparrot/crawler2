import glob
import sys


if __name__=="__main__":
    outputdir = sys.argv[1]
    lat1 = int(sys.argv[2])
    lat2 = int(sys.argv[3])
    lon1 = int(sys.argv[4])
    lon2 = int(sys.argv[5])

    spectra = glob.glob(outputdir+"/sbout.*")
    
    complete = True
    
    for jlat in range(lat1,lat2):
        for jlon in range(lon1,lon2):
            sfile = outputdir+"/sbout.%02d_%02d"%(jlat,jlon)
            if sfile not in spectra:
                print("Missing sbout.%02d_%02d"%(jlat,jlon))
                complete = False
                

