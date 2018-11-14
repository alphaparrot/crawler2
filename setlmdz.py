import os
import numpy as np
import time
from batch_system import SUB, BATCHSCRIPT

def edit_def(jid,filename,arg,val):
    f=open("lmdz/job"+jid+"/"+filename,"r")
    fnl = f.read().split('\n')
    found = False
    for l in range(0,len(fnl)):
        line = fnl[l].split('=')
        if len(line)>1:
            keyarg = line[0].replace(" ",'')
            if keyarg==arg:
                found=True
                fnl[l] = arg+'       = '+val
                break
    if not found:
        fnl.append(arg+'         = '+val)
    nml = '\n'.join(fnl)
    f=open("lmdz/job"+jid+"/"+filename,"w")
    f.write(nml)
    f.close()
    return

def edit_gases(jid,gasnames,gasfracs):
    gastext = ("# gases \n"+
               "%d      \n"%len(gasnames))
    for g in gasnames:
        gastext += g+"\n"
    
    for f in gasfracs:
        gastext += str(f)+"\n"
        
    gastext += "\n"
    
    gastext += ("# First line is number of gases              \n"+
                "# Followed by gas names (always 3 characters)\n"+
                "# and then mixing ratios.                    \n"+
                "# mixing ratio -1 means the gas is variable. \n")
    
    f=open("lmdz/job"+jid+"/gases.def","w")
    f.write(gastext)
    f.close()

    
def prep(job):
    
    sig = job.name
    jid = str(job.home)
    pid = job.pid
    args = job.args
    fields = job.fields
    
    workdir = job.top+"/lmdz/job"+jid
    
    if "template" in job.parameters:
        template = job.parameters["template"]
    else:
        template = "proxima"
  
    print "Setting stuff for job "+sig+" in lmdz/job"+jid+" which is task number "+pid
    print "Arguments are:",fields    
  
    scriptfile = "lmdz_loop.sh"
  
    os.system("cp "+job.top+"/lmdz/"+template+"/* "+workdir+"/")
    os.system("cp "+job.top+"/lmdz/"+scriptfile+" "+workdir+"/")
    
    recompile = False
    
    paramdict = np.load(job.top+"/lmdz/"+template+"/compiledfields.npy").item()
    
    daylen = paramdict["daylen"]
    preff = paramdict["preff"]
    radius = paramdict["radius"]
    gravity = paramdict["gravity"]
    surface = paramdict["surface"]
    
    co2 = 360.0
    if "co2" in job.parameters:
        co2 = float(job.parameters["co2"])
    
    nruns = 20
    
    if "daylen" in job.parameters: #seconds
        daylen = float(job.parameters["daylen"])
        recompile = True
    
    if "preff" in job.parameters: #Pa
        preff = float(job.parameters["preff"])
        fco2 = co2*0.1/preff
        fn2 = 1.0-fco2
        edit_gases(jid,["N2_","CO2","H2O"],[fn2,fco2,-1])
        recompile = True
        
    if "radius" in job.parameters: #in R_earth
        radius = float(job.parameters["radius"])
        recompile = True
        
    if "gravity" in job.parameters: #in m/s^2
        gravity = float(job.parameters["gravity"])
        recompile = True
    
    if "surface" in job.parameters:
        surface = job.parameters["surface"]
        recompile = True
    
    if "startemp" in job.parameters:
        edit_def(jid,"callphys.def","stelTbb",job.parameters["startemp"])
        edit_def(jid,"callphys.def","stelbbody",".true.")    
    
    if "corrkdir" in job.parameters:
        edit_def(jid,"callphys.def","corrkdir",job.parameters["corrkdir"])
        
    if "tlocked" in job.parameters:
        edit_def(jid,"callphys.def","tlocked",job.parameters["tlocked"])
        if job.parameters['tlocked']==".true.":
            edit_def(jid,"callphys.def","diurnal",".false.")

    if "nres" in job.parameters:
        edit_def(jid,"callphys.def","nres",job.parameters["nres"])
        
    if "season" in job.parameters:
        edit_def(jid,"callphys.def","season",job.parameters["season"])
        
    if "insol" in job.parameters:
        edit_def(jid,"callphys.def","Fat1AU",job.parameters["insol"])
        
    if "albedosnow" in job.parameters:
        edit_def(jid,"callphys.def","albedosnow",job.parameters["albedosnow"])
            
    if "maxicethick" in job.parameters:
        edit_def(jid,"callphys.def","maxicethick",job.parameters["maxicethick"])
        
    if "co2cond" in job.parameters:
        edit_def(jid,"callphys.def","co2cond",job.parameters["co2cond"])
        
    if "nday" in job.parameters:
        edit_def(jid,"run.def","nday",job.parameters["nday"])
        
    if "day_step" in job.parameters:
        edit_def(jid,"run.def","day_step",job.parameters['day_step'])
        
    if recompile:
        os.system("cd "+workdir+"/ && ./autostart.e "+str(preff)+" "+str(radius)+" "+
                  str(daylen)+" "+str(gravity)+" "+surface+" && cp restart.nc start.nc &&"+
                  "cp restartfi.nc startfi.nc && cd ../../")
        
    notify = "a"
    if "notify" in job.parameters:
        notify = job.parameters["notify"]
        
    if "nruns" in job.parameters:
        nruns = int(job.parameters["nruns"])
        
    print "Arguments and boundary conditions set."
    
    jobscript = (BATCHSCRIPT(job,notify)+
              "module load gcc/4.9.1                                          \n"+
              "module load python/2.7.9                                       \n"+
              "mkdir /mnt/node_scratch/paradise/lmdz_"+job.name+"/             \n"+
              "cp -a * /mnt/node_scratch/paradise/lmdz_"+job.name+"/           \n"+
              "cd /mnt/node_scratch/paradise/lmdz_"+job.name+"/                \n"+
              "./"+scriptfile+" "+job.name+" "+str(nruns)+"                        \n")
    
    rs = open(workdir+"/runlmdz","w")
    rs.write(jobscript)
    rs.close()
    
def submit(job):
    os.system("cd "+job.top+"/lmdz/job"+str(job.home)+" && "+SUB+" runlmdz && cd "+job.top)
    time.sleep(1.0)
    tag = job.getID()
    job.write()
    
        