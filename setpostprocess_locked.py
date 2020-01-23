import numpy as np
import os, sys
import netCDF4 as nc
from batch_system import SUB, BATCHSCRIPT
from identity import USER, SCRATCH
from crawldefs import Job

#This is really a minimal working example. The user has no choice in which fields to provide.
#Any missed fields, and the script will crash. It is however unusual in that the 'workdir' is
#not a folder called jobXX within the postprocess tree, and is probably instead in some *other*
#code's tree--in a way, a little parasitic.

def prep(job):
    workdir = job.parameters["workdir"]

    gtype = job.parameters["type"]
    gcm = job.parameters["gcm"]
    
    cwd = os.getcwd()
    
    homedir = cwd+"/postprocess_locked/job"+str(job.home)
    
    notify = 'ae'
    
    os.system("cp postprocess_locked/clean/* "+workdir+'/')
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
                    "cp "+cwd+"/postprocess_locked/job"+str(job.home)+"/job.npy ./           \n"+
                    "mkdir "+SCRATCH+"/postprocess_locked/                    \n"+
                    "mkdir "+SCRATCH+"/postprocess_locked/job"+str(job.home)+"/  \n"+
                    "cp -a * "+SCRATCH+"/postprocess_locked/job"+str(job.home)+"/  \n"+
                    "cd "+SCRATCH+"/postprocess_locked/job"+str(job.home)+"/  \n"+
                    "python postprocess_locked.py "+times+" "+angles+" "+tag+"                           \n"+
                    "cp -a "+SCRATCH+"/postprocess_locked/job"+str(job.home)+"/* "+workdir+"/ \n"+
                    "rm -rf *                                                         \n"+
                    "cd "+workdir+"                                                \n"+
                    "cp spectra.nc "+cwd+"/postprocess_locked/output/"+job.name+"_spectra.nc \n"+
                    "cp phases.nc "+cwd+"/postprocess_locked/output/"+job.name+"_phases.nc \n"+
                    "python release.py \n")
    else:
        tag+="phases "
        jobscript =(BATCHSCRIPT(job,notify)+
                    "module load gcc/4.9.1                                            \n"+
                    "module load python/2.7.9                                         \n"+
                    "cd "+workdir+"                      \n"+
                    "cp "+cwd+"/postprocess_locked/job"+str(job.home)+"/job.npy ./           \n"+
                    "mkdir "+SCRATCH+"/postprocess_locked/                    \n"+
                    "mkdir "+SCRATCH+"/postprocess_locked/job"+str(job.home)+"/  \n"+
                    "cp -a * "+SCRATCH+"/postprocess_locked/job"+str(job.home)+"/  \n"+
                    "cd "+SCRATCH+"/postprocess_locked/job"+str(job.home)+"/  \n"+
                    "python postprocess_locked.py "+times+" "+angles+" "+tag+"                                    \n"+
                    "cp -a "+SCRATCH+"/postprocess_locked/job"+str(job.home)+"/* "+workdir+"/ \n"+
                    "rm -rf *                                                         \n"+
                    "cd "+workdir+"                                                \n"+
                    "cp phases.nc "+cwd+"/postprocess_locked/output/"+job.name+"_phases.nc \n"+
                    "python release.py \n")
        
    
     
    rs = open(homedir+"/runpostprocess_locked","w")
    rs.write(jobscript)
    rs.close()

def submit(job):
    workdir = job.parameters["workdir"]
    cwd = os.getcwd()
    homedir = cwd+"/postprocess_locked/job"+str(job.home)
  
    if "DEPENDENCIES" in job.parameters:
        dlist = job.parameters["DEPENDENCIES"].split(',')
        priorjobs = []
        for d in dlist:
            with open(d+".id","r") as f:
                priorjobs.append(f.read().split('\n')[0].split()[0])
            os.system("echo %s >> "%job.pid+d+".id") #indicate that we depend on this job
        os.system("cd "+homedir+" && "+HOLD(priorjobs)+" runpostprocess_locked > %s/%s.id && cd ../../"%(job.top,job.pid))
    else:
        os.system("cd "+homedir+" && "+SUB+" runpostprocess_locked > %s/%s.id && cd ../../"%(job.top,job.pid))

def _prep(job):
    workdir = job[0]
    cwd = os.getcwd()
    top = cwd
    
    gtype = "plasim" #job.parameters["type"]
    jobname = job[1]
    gcm = top + "/plasim/output/"+jobname+"_snapshot.nc"
    
    homedir = cwd+"/postprocess_locked/"+jobname
    
    notify = 'ae'
    
    os.system("cp postprocess_locked/clean/* "+workdir+'/')
    os.system("cp release.py "+workdir+'/')
    os.system("cp crawldefs.py "+workdir+"/")
    os.system("cp identity.py "+workdir+"/")
    
    data = nc.Dataset(gcm,"r")
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:]
    np.save(workdir+"/latitudes.npy",lats)
    np.save(workdir+"/longitudes.npy",lons)
    
    times = job[2]
    angles = job[3]
    
    ncpu = '1'
    
    queue = job[4]
    
    notify = 'ae'
    
    color = True
    makemap = True
    
    tag = ''
    if color:
        tag+="color "
    if makemap:
        tag+="map "
            
    dummyjob = Job("# PID MODEL JOBNAME STATE NCORES QUEUE","%d pipeline ppc_%s 0 1 %s"%(-9999,jobname,queue),-1)
    jobscript =(BATCHSCRIPT(dummyjob,notify)+
                    "module load gcc/4.9.1                                            \n"+
                    "module load python/2.7.9                                         \n"+
                    "cd "+workdir+"                      \n"+
                    #"cp "+cwd+"/postprocess_locked/job"+str(job.home)+"/job.npy ./           \n"+
                    "mkdir "+SCRATCH+"/postprocess_locked/                    \n"+
                    "mkdir "+SCRATCH+"/postprocess_locked/"+jobname+"/  \n"+
                    "tar cvzf stuff.tar.gz *                                         \n"+
                    "rsync -avzhur stuff.tar.gz "+SCRATCH+"/postprocess_locked/"+jobname+"/  \n"+
                    "rm -rf stuff.tar.gz            \n"+
                    "cd "+SCRATCH+"/postprocess_locked/"+jobname+"/  \n"+
                    "tar xvzf stuff.tar.gz                \n"+
                    "rm -rf stuff.tar.gz                   \n"+
                    "python postprocess_locked.py "+times+" "+angles+" "+tag+"                           \n"+
                    "tar cvzf stuff.tar.gz *                 \n"+
                    "rsync -avzhur stuff.tar.gz "+workdir+"/ \n"+
                    "rm -rf *                                                         \n"+
                    "cd "+workdir+"                                                \n"+
                    "tar xvzf stuff.tar.gz                                \n"+
                    "rm stuff.tar.gz          \n"+
                    "cp spectra.nc "+cwd+"/postprocess_locked/output/"+jobname+"_spectra.nc \n"+
                    "cp phases.nc "+cwd+"/postprocess_locked/output/"+jobname+"_phases.nc \n\n"+
                    "cd "+cwd+"/postprocess_locked/output/   \n"+
                    "python orthoprojection.py "+jobname+"_phases.nc 0 \n"+
                    "mv "+jobname+"*/*.png . \n"+
                    "for p in $(ls "+jobname+"*.png)    \n"+
                    "do   \n"
                    "     mogrify -trim +repage $p    \n"+
                    "     mogrify -shave 1x1 +repage $p   \n"+
                    "done    \n"+
                    "rm -rf "+jobname+"*/   \n"+
                    "cp "+top+"/plasim/output/"+jobname+".nc .  \n"+
                    "cd "+workdir+"       \n"+
                    "python release.py     \n")

    rs = open(workdir+"/runpostprocess_locked","w")
    rs.write(jobscript)
    rs.close()

    
def _submit(job):
    workdir = job[0]
    cwd = os.getcwd()
    homedir = cwd+"/postprocess_locked/"+job[1]
    os.system("cd "+workdir+" && "+SUB+" runpostprocess_locked && cd ../../")
    

if __name__=="__main__":
    job = sys.argv[1:]
    _prep(job)
    _submit(job)
    