import os
import numpy as np

USER = "paradise"

MODELS = {"plasim":1,                #tasks per node (1 workq node on Sunnyvale has 8 threads)
          "sbdart":8,           #Here we use 'task' to mean a Sunnyvale job, as opposed to the
          "postprocess":8,      #HPC convention of a task being a thread or process. This way our
          "lmdz":8}             #code is MPI/OpenMP-agnostic.

def getjobs():
    os.system("qstat -u paradise > cjobs.tmp")
    cjf = open("cjobs.tmp","r")
    joblist = cjf.read().split('\n')[5:-1]
    cjf.close()
    os.system("rm cjobs.tmp")
    resources={}
    for m in MODELS.keys():
        resources[m] = np.zeros(256)
    tags = []
    for j in joblist:
        job = j.split()
        tags.append(job[0])
    for t in tags:
        os.system("qstat -f "+t+" > jinfo.tmp")
        jf = open("jinfo.tmp","r")
        jinfo = jf.read().split('\n')[1:-2]
        jf.close()
        os.system("rm jinfo.tmp")
        ncpus = 1
        for l in jinfo:
            if l.split()[0]=="init_work_dir":
                workdir = l.split()[2]
            if l.split()[0]=="Resource_List.ncpus":
                ncpus = int(l.split()[2])
        job = np.load(workdir+"/job.npy").item()
        jid = job.home
        if jid>=len(resources[job.model]):
            tmp = np.zeros(jid+100)
            tmp[:len(resources[job.model])] = resources[job.model][:]
            resources[job.model] = tmp
        resources[job.model][jid] = float(ncpus)/MODELS[job.model]
    
    return resources


class Job:
  def __init__(self,header,args,resource):
    self.args = args.split()
    self.pid   = self.args[0]
    self.model = self.args[1]
    self.name  = self.args[2]
    self.stat  = self.args[3]
    self.args = self.args[4:]
    self.ncores = int(self.args[0])
    self.queue = self.args[1]
    self.fields = header.split()[5:]
    self.top = os.getcwd()
    
    self.parameters = {}
    n=0
    for n in range(0,len(self.fields)):
      self.parameters[self.fields[n]] = self.args[n]
      
    self.home = resource
    self.jobname = self.name+".cl"
    
  def getID(self):
    os.system("qstat -u "+USER+" > cjobs.tmp")
    jf = open("cjobs.tmp","r")
    jlist = jf.read().split('\n')[5:-1]
    jf.close()
    os.system("rm cjobs.tmp")
    tag = None
    if len(jlist)>0:
      for j in jlist:
        job = j.split()
        name = job[3]
        if name==self.jobname:
          tag = job[0]
          break
    self.tag = tag
    return tag

  def write(self):
    jf = open(self.model+"/job"+str(self.home)+"/job.crwl","w")
    try:
        jt = self.pid+'\n'+' '.join(self.fields)+'\n'+' '.join(self.args)+'\n'+self.name+'\n'+self.tag
    except:
        jt = self.pid+'\n'+' '.join(self.fields)+'\n'+' '.join(self.args)+'\n'+self.name+'\nFAILED'
    jf.write(jt)
    jf.close()

  def kill(self):
    tag = self.getID()
    if tag:
      os.system("qdel "+tag)
      return 1
    else:
      return 0