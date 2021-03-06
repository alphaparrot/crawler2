import os
import time

def newtask(job,dryrun=False):
  modeldir = job.model+"/"
  workdir  = modeldir+"job"+str(job.home)+"/"
  
  modelroutines = "import set"+job.model+" as setjob"
  exec(modelroutines)
  
  #Create the working directory if it doesn't already exist, and clean it if necessary
  try:
    os.system("mkdir "+workdir)
  except:
    pass
  #try:
  #  os.system("rm -rf "+workdir+"*")
  #except:
  #  pass
  
  #os.system("cp -r "+modeldir+"clean/* "+workdir)
  os.system("cp release.py "+workdir)
  os.system("cp crawldefs.py "+workdir)
  
  setjob.prep(job)
  if not dryrun:
    setjob.submit(job)