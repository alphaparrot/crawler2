import os
import glob
import crawlset
import sys
import numpy as np
from crawldefs import *


if __name__=="__main__":    
  #Check which resources are in use and compare to our max allotment
  
  if "DRYRUN" in sys.argv[:]:
      dryrun=True
  else:
      dryrun=False
  
  f=open("nnodes.crwl","r")
  nnodes=float(f.read().split('\n')[0])
  f.close()
  #resources={}
  #for m in MODELS.keys():
    #rf=open("running_"+m+".crwl","r")
    #resources[m] = rf.read().split('\n')[0]
    #resources[m] = resources[m].split()
    ##print resources[m]
    #rf.close()
    
  resources = getjobs()  
    
  running = 0
  for r in resources.keys():
    #for n in range(0,len(resources[r])):
      #running+=int(resources[r][n])*1.0/MODELS[r]
    running += np.sum(resources[r])  
  print str(running)+" of "+str(nnodes)+" nodes used"
  
  priority = True
  
  while running < nnodes and priority: #First do jobs placed in priority.crwl.
    
    #Get next task
    #NOTE JOBS IN PRIORITY.CRWL SHOULD BE ASSIGNED NEGATIVE PIDs TO AVOID COLLISIONS
    
    f=open("priority.crwl","r")
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
              queued=True
              qpid=int(task[0])
              taskmodel = task[1]
              if running+(1.0/MODELS[taskmodel])<=nnodes:  #We might be at 1 cpu less than capacity,
                mark=i                                     #so a plasim job might put us over the limit.
                f=open("tasklog.crwl","a")
                f.write("\nFound job "+' '.join(task))
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
            if (len(header.split())-1)==(len(task)): #first header with the right number of args
              break
            else:
              header=False
        if not header:                #None match!
          f=open("tasklog.crwl","a")
          f.write("\nTask "+str(qpid)+" header mismatch; skipping")
          f.close()
        else:                         #We found the header. Onward!
          ready=True
          f=open("tasklog.crwl","a")
          f.write("\nTask "+str(qpid)+" header found")
          f.close()          
      else:                           #We didn't find an available task
        f=open("tasklog.crwl","a")
        f.write("\nNo open jobs in task queue; all done!")
        print "No jobs in priority queue"
        f.close()
        #running=nnodes
        priority = False
        break

    #Engage next task
    if ready:
      print "Engaging a priority task"
      
      goahead = False
      newjob = Job(header,' '.join(task),-1)       #Collect and organize the job parameters
      
      for i in range(0,len(resources[taskmodel])):
        if float(resources[taskmodel][i])==0.0:
          rid = i
          print resources[taskmodel]
          resources[taskmodel][i]=newjob.ncores/float(MODELS[taskmodel])
          print resources[taskmodel]
          newjob.home = rid
          goahead = True
          break
      if not goahead:                 #No open slot found for this model
        if running+float(newjob.ncores)/MODELS[taskmodel] <= nnodes: #We do have space for one more though
          print resources[taskmodel]
          tmp = np.zeros(len(resources[taskmodel])+100)
          tmp[:len(resources[taskmodel])] = resources[taskmodel][:]
          n0 = len(resources[taskmodel])
          resources[taskmodel] = tmp
          resources[taskmodel][n0] = newjob.ncores/float(MODELS[taskmodel])
          print resources[taskmodel]
          rid = n0
          newjob.home = rid
          goahead = True
        else:                                       #Nope, pack up and go home
          print "At capacity."
          priority = False
      if goahead:                           #Found a job slot
        task[3] = '1'
        task = ' '.join(task)
        tasks[mark] = task
        tasks='\n'.join(tasks)
        f=open("priority.crwl","w")
        f.write(tasks)
        f.close()
      
        crawlset.newtask(newjob,dryrun=dryrun)            #Set up the job and submit it
        np.save(taskmodel+'/job'+str(newjob.home)+'/job.npy',newjob)
        if not dryrun:
            newjob.getID()
        else:
            newjob.tag='xxxxx.doug'
        newjob.write()
        f=open("tasklog.crwl","a")
        f.write("\nQueued job "+newjob.name+" in "+newjob.queue+" for "+str(newjob.ncores)+" cores.")
        f.close()
        running += float(newjob.ncores)/MODELS[taskmodel]    #Note that we are *that* much closer to the limit.
  
  
  print "Moving on to normal tasks"
  print running,nnodes
  while running < nnodes: #We are using less than our full allocation, and the priority list is empty.
    
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
                f=open("tasklog.crwl","a")
                f.write("\nFound job "+' '.join(task))
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
            if (len(header.split())-1)==(len(task)): #first header with the right number of args
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
      
      goahead = False
      newjob = Job(header,' '.join(task),-1)       #Collect and organize the job parameters
      
      for i in range(0,len(resources[taskmodel])):
        if float(resources[taskmodel][i])==0.0:
          rid = i
          print resources[taskmodel]
          resources[taskmodel][i]=newjob.ncores/float(MODELS[taskmodel])
          print resources[taskmodel]
          newjob.home = rid
          goahead = True
          break
      if not goahead:                 #No open slot found for this model
        if running+float(newjob.ncores)/MODELS[taskmodel] <= nnodes: #We do have space for one more though
          print resources[taskmodel]
          tmp = np.zeros(len(resources[taskmodel])+100)
          tmp[:len(resources[taskmodel])] = resources[taskmodel][:]
          n0 = len(resources[taskmodel])
          resources[taskmodel] = tmp
          resources[taskmodel][n0] = newjob.ncores/float(MODELS[taskmodel])
          print resources[taskmodel]
          rid = n0
          newjob.home = rid
          goahead = True
        else:                                       #Nope, pack up and go home
          print "At capacity."
      if goahead:                           #Found a job slot
        task[3] = '1'
        task = ' '.join(task)
        tasks[mark] = task
        tasks='\n'.join(tasks)
        f=open("tasks.crwl","w")
        f.write(tasks)
        f.close()
        
        crawlset.newtask(newjob,dryrun=dryrun)            #Set up the job and submit it
        np.save(taskmodel+'/job'+str(newjob.home)+'/job.npy',newjob)
        if not dryrun:
            newjob.getID()
        else:
            newjob.tag='xxxxx.doug'
        newjob.write()
        running += float(newjob.ncores)/MODELS[taskmodel]    #Note that we are *that* much closer to the limit.
        
  #folders = []
  ##print MODELS.keys()
  #for m in MODELS.keys():
    #fr=open("running_"+m+".crwl","w")
    ##print resources[m]
    ##print '--'
    #fr.write(' '.join(resources[m])+'\n')
    #fr.close()
    #folders+=sorted(glob.glob(m+"/job*/"))
  
  ##Keep track of the current jobs
  
  #htext = "Current job array:\n"
  #for f in folders:
    #jf = open(f+"/job.crwl")
    #jc = jf.read().split('\n')
    #jf.close()
    #htext+=f+" "+jc[3]+"\n"
  #cjb = open("currentjobs.crwl","w")
  #cjb.write(htext)
  #cjb.close()
  
  os.system("echo '0'>inuse.crwl")         #Release ownership of crawler.py
     