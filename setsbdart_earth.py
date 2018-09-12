import numpy as np
import netCDF4 as nc
import os
import time
from batch_system import SUB, BATCHSCRIPT

#This interface won't actually build the job--instead it'll spawn a sunnyvale worker that
#will pick up the saved job info, build the job, and then submit it to the queue.

def prep(job):
  workdir = job.top+"/sbdart_earth/job"+str(job.home)
  os.system("mkdir "+workdir)
  np.save(workdir+"/jobdat.npy",job)
  os.system("cp sbdart_earth/buildsbdart.py "+workdir+"/")
  os.system("cp crawldefs.py "+workdir+"/")
  os.system("cp identity.py "+workdir+"/")
  os.system("cp batch_system.py "+workdir+"/")
  os.system("cp torque.py "+workdir+"/")
  os.system("cp slurm.py "+workdir+"/")
  job.name = "build_"+job.name
  jobscript =(BATCHSCRIPT(job,'a')+
              "module load gcc/4.9.1                                          \n"+
              "module load python/2.7.9                                       \n"+
              "python buildsbdart.py   \n")
  
  rs = open(workdir+"/runsbdart_earth","w")
  rs.write(jobscript)
  rs.close()
          
 
def submit(job):
  workdir = "sbdart_earth/job"+str(job.home)
  
  os.system("cd "+workdir+" && "+SUB+" runsbdart_earth && cd "+job.top)
