import os
import glob
import crawlset
import sys

USER = "paradise"

MODELS = {"plasim":1,                #tasks per node (1 workq node on Sunnyvale has 8 threads)
          "sbdart":8,           #Here we use 'task' to mean a Sunnyvale job, as opposed to the
          "postprocess":8,      #HPC convention of a task being a thread or process. This way our
          "lmdz":8}             #code is MPI/OpenMP-agnostic.
class Job:
  def __init__(self,header,args,resource):
    self.args = args.split()
    self.model = self.args[0]
    self.pid   = self.args[1]
    self.name  = self.args[2]
    self.stat  = self.args[3]
    self.args = self.args[4:]
    self.fields = header.split()[5:]
    
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
    jf = open(self.model+"/job"+self.home+"/job.crwl","w")
    jt = self.pid+'\n'+' '.join(self.fields)+'\n'+' '.join(self.args)+'\n'+self.name+'\n'+self.tag
    jf.write(jt)
    jf.close()

  def kill(self):
    tag = self.getID()
    if tag:
      os.system("qdel "+tag)
      return 1
    else:
      return 0


if __name__=="__main__":    
  #Check which resources are in use and compare to our max allotment
  
  if "DRYRUN" in sys.argv[:]:
      dryrun=True
  else:
      dryrun=False
  
  f=open("nnodes.crwl","r")
  nnodes=int(f.read().split('\n')[0])
  f.close()
  resources={}
  for m in MODELS:
    rf=open("running_"+m+".crwl","r")
    resources[m] = rf.read().split('\n')[:-1]
    rf.close()
    
  running = 0
  for r in resources:
    for n in range(0,len(resources[r])):
      running+=int(resources[r][n])*1.0/MODELS[r]
      
  print str(running)+" of "+str(nnodes)+" nodes used"
  

  while running < nnodes: #We are using less than our full allocation
    
    #Get next task
    
    f=open("tasks.crwl","r")
    tasks=f.read().split('\n')
    f.close()
    queued=False
    ready=False
    while not ready:       #Search for a job to run 
      for i in range(0,len(tasks)-1):
        if tasks[i]!='':
          if tasks[i][0]!="#":
            task = tasks[i].split()
            if int(task[3])==0:
              queued=int(task[0])
              taskmodel = task[1]
              if running+(1.0/MODELS[taskmodel])<=nnodes:  #We might be at 1 cpu less than capacity,
                mark=i                                     #so a plasim job might put us over the limit.
                task[3] = '1'
                task = ' '.join(task)
                f=open("tasklog.crwl","a")
                f.write("\nFound job "+task)
                f.close()
                break
              else:
                queued = False
      ready=False
      if queued:
        f=open("tasklog.crwl","a")
        f.write("\nSearching for header...")
        f.close()
        header=False
        for i in range(mark,-1,-1): #Find header for the task
          if tasks[i][0]=="#":
            header = tasks[i]
            if (len(header.split())-1)==(len(task.split())): #first header with the right number of args
              break
            else:
              header=False
        if not header:                #None match!
          f=open("tasklog.crwl","a")
          f.write("\nTask "+str(queued)+" header mismatch; skipping")
          f.close()
        else:                         #We found the header. Onward!
          ready=True
          f=open("tasklog.crwl","a")
          f.write("\nTask "+str(queued)+" header found")
          f.close()          
      else:                           #We didn't find an available task
        f=open("tasklog.crwl","a")
        f.write("\nNo open jobs in task queue; all done!")
        f.close()
        running=nnodes
        break

    #Engage next task
    if ready:
      tasks[mark] = task
      tasks='\n'.join(tasks)
      f=open("tasks.crwl","w")
      f.write(tasks)
      f.close()
      
      goahead = False
      
      for i in range(0,len(resources[taskmodel])):
        if int(resources[taskmodel][i])==0:
          rid = i
          resources[taskmodel][i]="1"
          goahead = True
          break
      if not goahead:                 #No open slot found for this model
        if running+1.0/MODELS[taskmodel] <= nnodes: #We do have space for one more though
          resources[taskmodel].append("1")
          rid = len(resources[taskmodel])-1
          goahead = True
        else:                                       #Nope, pack up and go home
          f=open("running_"+MODELS[taskmodel]+".crwl","w")
          f.write('\n'.join(resources[taskmodel]))
          f.close()
      if goahead:                           #Found a job slot
        newjob = Job(header,task,rid)       #Collect and organize the job parameters
        crawlset.newtask(newjob,dryrun=dryrun)            #Set up the job and submit it
        np.save(taskmodel+'/job'+newjob.home+'/job.npy',newjob)
        if not dryrun:
            newjob.getID()
        else:
            newjob.tag='xxxxx.doug'
        newjob.write()
        running += 1.0/MODELS[taskmodel]    #Note that we are *that* much closer to the limit.
        
    folders = []
    for m in MODELS:
      fr=open("running_"+m+".crwl","w")
      fr.write('\n'.join(resources[m]))
      fr.close()
      folders+=sorted(glob.glob(m+"/job*/"))
    
    #Keep track of the current jobs
    
    htext = "Current job array:\n"
    for f in folders:
      jf = open(f+"/job.crwl")
      jc = jf.read().split('\n')
      jf.close()
      htext+=f+" "+jc[3]+"\n"
    cjb = open("currentjobs.crwl","w")
    cjb.write(htext)
    cjb.close()
    
    os.system("echo 0>inuse.crwl")         #Release ownership of crawler.py
     