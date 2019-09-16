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
    
    homedir = cwd+"/postprocess_earth/job"+str(job.home)
    
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
    
    if "times" in job.parameters:
        times = job.parameters["times"]
    else:
        times="0"
        
    if "angles" in job.parameters:
        angles = job.parameters["angles"]
    else:
        angles="Z"
    
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
                    "cd "+workdir+"                      \n"+
                    "cp "+cwd+"/postprocess_earth/job"+str(job.home)+"/job.npy ./           \n"+
                    "mkdir /mnt/node_scratch/paradise/postprocess_earth/                    \n"+
                    "mkdir /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/  \n"+
                    "tar cvzf stuff.tar.gz *                                         \n"+
                    "rsync -avzhur stuff.tar.gz /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/  \n"+
                    "rm -rf stuff.tar.gz            \n"+
                    "cd /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/  \n"+
                    "tar xvzf stuff.tar.gz                \n"+
                    "rm -rf stuff.tar.gz                   \n"+
                    "python postprocess_earth.py "+times+" "+angles+" "+tag+"                           \n"+
                    "tar cvzf stuff.tar.gz *                 \n"+
                    "rsync -avzhur stuff.tar.gz "+workdir+"/ \n"+
                    "rm -rf *                                                         \n"+
                    "cd "+workdir+"                                                \n"+
                    "tar xvzf stuff.tar.gz                                \n"+
                    "rm stuff.tar.gz          \n"+
                    "cp spectra.nc "+cwd+"/postprocess_earth/output/"+job.name+"_spectra.nc \n"+
                    "cp phases.nc "+cwd+"/postprocess_earth/output/"+job.name+"_phases.nc \n"+
                    "python release.py \n")
    else:
        tag+="phases "
        jobscript =(BATCHSCRIPT(job,notify)+
                    "module load gcc/4.9.1                                            \n"+
                    "module load python/2.7.9                                         \n"+
                    "cd "+workdir+"                      \n"+
                    "cp "+cwd+"/postprocess_earth/job"+str(job.home)+"/job.npy ./           \n"+
                    "mkdir /mnt/node_scratch/paradise/postprocess_earth/                    \n"+
                    "mkdir /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/  \n"+
                    "tar cvzf stuff.tar.gz *                                         \n"+
                    "rsync -avzhur stuff.tar.gz /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/  \n"+
                    "rm -rf stuff.tar.gz            \n"+
                    "cd /mnt/node_scratch/paradise/postprocess_earth/job"+str(job.home)+"/  \n"+
                    "tar xvzf stuff.tar.gz                \n"+
                    "rm -rf stuff.tar.gz                   \n"+
                    "python postprocess_earth.py "+times+" "+angles+" "+tag+"                                    \n"+
                    "tar cvzf stuff.tar.gz *                 \n"+
                    "rsync -avzhur stuff.tar.gz "+workdir+"/ \n"+
                    "rm -rf *                                                         \n"+
                    "cd "+workdir+"                                                \n"+
                    "tar xvzf stuff.tar.gz                                \n"+
                    "rm stuff.tar.gz          \n"+
                    "cp phases.nc "+cwd+"/postprocess_earth/output/"+job.name+"_phases.nc \n"+
                    "python release.py \n")
        
    
     
    rs = open(homedir+"/runpostprocess_earth","w")
    rs.write(jobscript)
    rs.close()

def submit(job):
    workdir = job.parameters["workdir"]
    cwd = os.getcwd()
    homedir = cwd+"/postprocess_earth/job"+str(job.home)
  
    os.system("cd "+homedir+" && "+SUB+" runpostprocess_earth && cd ../../")