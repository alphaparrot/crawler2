import os
import time

def newtask(job):
  modeldir = job.model+"/"
  workdir  = modeldir+"job"+str(job.home)+"/"
  
  modelroutines = "import set"+job.model+" as setjob"
  exec(modelroutines)
  
  #Create the working directory if it doesn't already exist, and clean it if necessary
  try:
    os.system("mkdir "+workdir)
  except:
    pass
  try:
    os.system("rm -rf "+workdir+"*")
  except:
    pass
  
  os.system("cp -r "+modeldir+"clean/* "+workdir)
  os.system("cp packup.py "+workdir)
  
  setjob.prep(job)
  setjob.submit(job)