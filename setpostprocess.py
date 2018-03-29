import numpy as np
import os
import netCDF4 as nc

#This is really a minimal working example. The user has no choice in which fields to provide.
#Any missed fields, and the script will crash. It is however unusual in that the 'workdir' is
#not a folder called jobXX within the postprocess tree, and is probably instead in some *other*
#code's tree--in a way, a little parasitic.

def prep(job):
    workdir = job.parameters["workdir"]
    gtype = job.parameters["type"]
    gcm = job.parameters["gcm"]
    
    cwd = os.getcwd()
    
    os.system("cp postprocess/clean/postprocess.py "+workdir+'/')
    os.system("cp release.py "+workdir+'/')
    os.system("cp crawldefs.py "+workdir+"/")
    
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
    
    jobscript =("#!/bin/bash -l                                                  \n"+
                "#PBS -l nodes=1:ppn=1                                            \n"+
                "#PBS -q workq                                                    \n"+
                "#PBS -r n                                                        \n"+
                "#PBS -l walltime=48:00:00                                        \n"+
                "#PBS -m abe                                                      \n"+
                "#PBS -N "+job.name+"                                             \n"
                "# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE"+
                " nodes,ppn,walltime and my_job_name VALUES                       \n"+
                "cd $PBS_O_WORKDIR                                                \n"+
                "module load gcc/4.9.1                                            \n"+
                "module load python/2.7.9                                         \n"+
                "python postprocess.py "+lon0+"                                   \n"+
                "cp spectra.nc "+cwd+"/postprocess/output/"+job.name+"_spectra.nc \n"+
                "cp phasecurve.nc "+cwd+"/postprocess/output/"+job.name+"_phasecurve.nc \n"+
                "mv "+cwd+"/postprocess/job"+str(job.home)+"/job.npy ./           \n"+
                "python release.py \n")
    
     
    rs = open(workdir+"/runpostprocess","w")
    rs.write(jobscript)
    rs.close()

def submit(job):
    workdir = job.parameters["workdir"]
  
    os.system("cd "+workdir+" && qsub runpostprocess && cd ../../")