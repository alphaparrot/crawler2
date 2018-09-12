import glob
import sys
import os
import time
import random
from crawldefs import Job

def addtopriority(top,outputdir,jlat,jlon):
    inuse = True
    actime = 0.0
    
    job = np.load("job.npy").item()
    
    while inuse:
        inf=open(top+"/inuse.crwl","r")
        if int(inf.read()[0])=="0":
            inuse=False
            os.system("echo '1'>"+top+"/inuse.crwl")
        else:
            nt = random.random()
            actime +=nt
            time.sleep(nt)
            if actime>600.0:
                break
    if not inuse:
      taskf = open(top+"/priority.crwl","r")
      tasks = taskf.read().split('\n')
      taskf.close()
      last = tasks[-1]
      n=-1
      while last=='':
          n-=1
          last = tasks[n]
          if n==-len(tasks):
              last = '0 DUMMY'
              break
      name = outputdir.split('/')[1]
      lastpid = int(last.split()[0])
      head = '# PID MODEL JOBNAME STATUS '+' '.join(job.fields)
      job.args[4] = str(jlat)+","+str(jlat+1)
      job.args[5] = str(jlon)+","+str(jlon+1)
    
      jobt = str(lastpid-1)+" sbdart "+name+"%02d 0 "%(lastpid-1)+" ".join(job.args)

      os.system("echo '"+head+"' >> "+top+"/priority.crwl")
      os.system("echo '"+jobt+"' >> "+top+"/priority.crwl")
      os.system("echo '0'>"+top+"/inuse.crwl")

if __name__=="__main__":
    outputdir = sys.argv[1]
    lat1 = int(sys.argv[2])
    lat2 = int(sys.argv[3])
    lon1 = int(sys.argv[4])
    lon2 = int(sys.argv[5])
    lon0 = sys.argv[6]
    top = sys.argv[7]
    model = sys.argv[8]
    gcm = sys.argv[9]
    
    spectra = glob.glob(outputdir+"/sbout.*")
    
    complete = True
    
    for jlat in range(lat1,lat2):
        for jlon in range(lon1,lon2):
            sfile = outputdir+"/sbout.%02d_%02d"%(jlat,jlon)
            if sfile not in spectra:
                print "Missing sbout.%02d_%02d"%(jlat,jlon)
                addtopriority(top,outputdir,jlat,jlon)
                complete = False
                
    
    #if complete:
        #taskf = open(top+"/priority.crwl","r")
        #tasks = taskf.read().split('\n')
        #taskf.close()
        
        #last = tasks[-1]
        #n = -1
        #while last=='':
            #n-=1
            #last = tasks[n]
            #if n==-len(tasks):
                #last = '0 DUMMY'
                #break
        
        #name = outputdir.split('/')[1]
        #lastpid = int(last.split()[0])
        #head = "# PID MODEL JOBNAME STATUS NCORES QUEUE type gcm workdir lon0"
        #jobt = str(lastpid+1)+" postprocess "+name+" 0 1 workq "+model+" "+gcm+" "+outputdir+' '+lon0
        
        #usef = open(top+"/inuse.crwl","r")
        #inuse = int(usef.read()[0])
        #usef.close()
        #if inuse==0:
            #os.system("echo '1'>"+top+"/inuse.crwl")
            #os.system('echo "'+head+'">>'+top+"/priority.crwl")
            #os.system('echo "'+jobt+'">>'+top+"/priority.crwl")
            #os.system("echo '0'>"+top+"/inuse.crwl")
            
    #Clean up
    os.system("python release.py")
