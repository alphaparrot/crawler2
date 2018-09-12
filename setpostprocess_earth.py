import numpy as np
import os
import netCDF4 as nc
from batch_system import SUB, BATCHSCRIPT

#This is really a minimal working example. The user has no choice in which fields to provide.
#Any missed fields, and the script will crash. It is however unusual in that the 'workdir' is
#not a folder called jobXX within the postprocess tree, and is probably instead in some *other*
#code's tree--in a way, a little parasitic.

def prep(job):
    workdir = job.parameters["workdir"]
    gtype = job.parameters["type"]
    gcm = job.parameters["gcm"]
    
    cwd = os.getcwd()
    
    notify = 'ae'
    
    os.system("cp postprocess_earth/clean/* "+workdir+'/')
    os.system("cp release.py "+workdir+'/')
    os.system("cp crawldefs.py "+workdir+"/")
    os.system("cp identity.py "+workdir+"/")
    
    if gtype=="plasim":
        data = nc.Dataset("hopper/"+gcm,"r")
        lats = data.variables['lat'][:]
        lons = data.variables['lon'][:]
    if gtype=="lmdz":
        data = np.load("hopper/"+gcm).item()
        #Or do we need to use rlatv,rlonu? latitude and longitude give us edges I think. But the
        #grid seems to be discretized that way too?
        lats = data['latitude'][1:-1] #Get rid of poles
        lons = data['longitude'][:-1] #Get rid of redundant longitude
    np.save(workdir+"/latitudes.npy",lats)
    np.save(workdir+"/longitudes.npy",lons)
    
    if "lon0" in job.parameters:
        lon0 = job.parameters["lon0"]
    else:
        lon0 = "0"
    
    if "NCORES" in job.parameters:
        ncpu = job.parameters["NCORES"]
    else:
        ncpu = '1'
        
    if "QUEUE" in job.parameters:
        queue = job.parameters["QUEUE"]
    else:
        queue = "workq"
        
    if "notify" in job.parameters:
        notify = job.parameters['notify']
        
    part2 = False
    if "part2" in job.parameters:
        if job.parameters["part2"]=="0":
            part2=False
        else:
            part2=True
    
    color=False
    if "color" in job.parameters:
        if job.parameters["color"]=="True":
            color=True
        else:
            color=False
    
    makemap=False
    if "map" in job.parameters:
        if job.parameters["map"]=="True":
            makemap=True
        else:
            makemap=False
    
    tag = ''
    if color:
        tag+="color "
    if makemap:
        tag+="map "
            
    if not part2:
        jobscript =(BATCHSCRIPT(job,notify)+
                    "module load gcc/4.9.1                                            \n"+
                    "module load python/2.7.9                                         \n"+
                    "mv "+cwd+"/postprocess/job"+str(job.home)+"/job.npy ./           \n"+
                    "mkdir /mnt/node_scratch/paradise/postprocess_earth/                    \n"+
                    "mkdir /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/  \n"+
                    "cp -a * /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/  \n"+
                    "cd /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/  \n"+
                    "python postprocess.py "+lon0+" "+tag+"                           \n"+
                    "cp -a /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/* $PBS_O_WORKDIR/ \n"+
                    "rm -rf *                                                         \n"+
                    "cd $PBS_O_WORKDIR                                                \n"+
                    "cp spectra.nc "+cwd+"/postprocess_earth/output/"+job.name+"_spectra.nc \n"+
                    "cp phasecurve.nc "+cwd+"/postprocess_earth/output/"+job.name+"_phasecurve.nc \n"+
                    "python release.py \n")
    else:
        jobscript =(BATCHSCRIPT(job,notify)+
                    "module load gcc/4.9.1                                            \n"+
                    "module load python/2.7.9                                         \n"+
                    "mv "+cwd+"/postprocess_earth/job"+str(job.home)+"/job.npy ./           \n"+
                    "mkdir /mnt/node_scratch/paradise/postprocess_earth/                    \n"+
                    "mkdir /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/  \n"+
                    "cp -a * /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/  \n"+
                    "cd /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/  \n"+
                    "python postspectra.py "+tag+"                                    \n"+
                    "cp -a /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/* $PBS_O_WORKDIR/ \n"+
                    "rm -rf *                                                         \n"+
                    "cd $PBS_O_WORKDIR                                                \n"+
                    "cp phasecurve.nc "+cwd+"/postprocess_earth/output/"+job.name+"_phasecurve.nc \n"+
                    "python release.py \n")
        
    
     
    rs = open(workdir+"/runpostprocess_earth","w")
    rs.write(jobscript)
    rs.close()

def submit(job):
    workdir = job.parameters["workdir"]
  
    os.system("cd "+workdir+" && "+SUB+" runpostprocess_earth && cd ../../")