import numpy as np
import netCDF4 as nc
import os
import time
from batch_system import SUB, BATCHSCRIPT, HOLD

#This interface won't actually build the job--instead it'll spawn a sunnyvale worker that
#will pick up the saved job info, build the job, and then submit it to the queue.

def prep(job):
  workdir = job.top+"/sbdart_earth/job"+str(job.home)
  os.system("mkdir "+workdir)
  np.save(workdir+"/jobdat.npy",job)
  os.system("cp sbdart_earth/buildsbdart.py "+workdir+"/")
  os.system("cp crawldefs.py "+workdir+"/")
  os.system("cp crawldefs.py sbdart_earth/")
  os.system("cp identity.py "+workdir+"/")
  os.system("cp batch_system.py "+workdir+"/")
  os.system("cp batch_system.py sbdart_earth/")
  os.system("cp torque.py "+workdir+"/")
  os.system("cp torque.py sbdart_earth/")
  os.system("cp slurm.py "+workdir+"/")
  os.system("cp slurm.py sbdart_earth/")
  #job.name = "build_"+job.name
  jobscript =(BATCHSCRIPT(job,'a')+
              "module load gcc/4.9.1                                          \n"+
              "module load python/2.7.9                                       \n"+
              "python buildsbdart.py   \n")
  
  rs = open(workdir+"/runsbdart_earth","w")
  rs.write(jobscript)
  rs.close()
          
 
def submit(job):
  workdir = "sbdart_earth/job"+str(job.home)
  
  if "DEPENDENCIES" in job.parameters:
      dlist = job.parameters["DEPENDENCIES"].split(',')
      priorjobs = []
      for d in dlist:
          with open(d+".id","r") as f:
              priorjobs.append(f.read().split('\n')[0].split()[0])
          os.system("echo %d >> "%job.pid+d+".id") #indicate that we depend on this job
      os.system("cd "+workdir+" && "+HOLD(priorjobs)+" runsbdart_earth > %s/%d.id && cd "%(job.top,job.pid)+job.top)
  else:
      os.system("cd "+workdir+" && "+SUB+" runsbdart_earth > %s/%d.id && cd "%(job.top,job.pid)+job.top)