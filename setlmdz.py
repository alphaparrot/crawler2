import os
import numpy as np
import time

jobscript1 =("#!/bin/bash -l                                                  \n"+
            "#PBS -l nodes=1:ppn=8                                            \n"+
            "#PBS -q workq                                                    \n"+
            "#PBS -r n                                                        \n"+
            #"#PBS -m ae                                                       \n"+
            "#PBS -l walltime=48:00:00                                        \n")
jobscript2 =("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE"+
             " nodes,ppn,walltime and my_job_name VALUES                       \n"+
            "cd $PBS_O_WORKDIR                                                \n"+
            "module unload intel/intel-17                                     \n"+
            "module unload openmpi/2.0.1-intel-17                             \n"+
            "module load gcc/4.9.1                                            \n"+
            "module load python/2.7.9                                         \n")

def edit_def(jid,filename,arg,val):
    f=open("lmdz/job"+jid+"/"+filename,"r")
    fnl = f.read().split('\n')
    found = False
    for l in range(0,len(fnl)):
        line = fnl[l].split('=')
        if len(line)>1:
            keyarg = fnl[l][0].split()[0]
            if keyarg==arg:
                found=True
                fnl[l] = arg+'       = '+val
                break
    if not found:
        fnl.append(arg+'         = '+val)
    nml = '\n'.join(fnl)
    f=open("lmdz/job"+jid+"/"+filename,"r")
    f.write(nml)
    f.close()
    return

def prep(job):
    
    sig = job.name
    jid = job.home
    pid = job.pid
    args = job.args
    fields = job.fields
    
    if "template" in job.parameters:
        template = job.parameters["template"]
    else:
        template = "proxima"
  
    print "Setting stuff for job "+sig+" in lmdz/job"+jid+" which is task number "+pid
    print "Arguments are:",fields    
  
    emailtag = "#PBS -m ae \n"
    scriptfile = "run.sh"
  
    os.system("cp lmdz/"+template+"/* lmdz/job"+jid+"/")
    os.system("cp lmdz/"+scriptfile+" lmdz/job"+jid+"/")
    
    for name in job.fields:
        val = job.parameters[name]
        found = False
        
        