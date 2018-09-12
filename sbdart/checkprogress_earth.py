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
      job.args[2] = "dom"
      job.args[5] = str(jlat)+","+str(jlat+1)
      job.args[6] = str(jlon)+","+str(jlon+1)
    
      jobt = str(lastpid-1)+" sbdart "+name+"%02d 0 "%(lastpid-1)+" ".join(job.args)

      os.system("echo '"+head+"' >> "+top+"/priority.crwl")
      os.system("echo '"+jobt+"' >> "+top+"/priority.crwl")
      f=open(top+"/inuse.crwl","w")
      f.write("0")
      f.close()
      
def nearest(array,val,lower=True): #if lower, we'll take the closest index that has a lower value
    diff=1.0e10                    #else, we'll find the index with the closest value
    index = -1
    for n in range(0,len(array)): 
        dd = val-array[n]
        if not lower:
            dd = abs(dd)
        if dd<diff and dd>0:
            diff=dd
            index = n
    return index

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
    
    color=False
    makemap=False
    if "color" in sys.argv[:]:
        color=True
    if "map" in sys.argv[:]:
        makemap=True
    
    spectra = glob.glob(outputdir+"/sbout.*")
    
    complete = True
    
    dtokens = sorted(glob.glob(outputdir+"/finished/token*.crwl"))
    rtokens = sorted(glob.glob(outputdir+"/running/token*.crwl"))
    
    dpids = []
    rpids = []
    
    ln1d = []
    ln2d = []
    
    ln1r = []
    ln2r = []
    
    lt1d = []
    lt2d = []
    
    lt1r = []
    lt2r = []
    
    if len(dtokens)>0:
    
        for t in dtokens:
            name = t.split('/')[-1]
            name = name.split('_')
            coords = name[1].split('-')
            dpids.append(int(name[0].split("token")[-1]))
            ln1d.append(int(coords[0]))
            ln2d.append(int(coords[1]))
            lt1d.append(int(coords[2]))
            lt2d.append(int(coords[3]))
            
        if len(rtokens)>0:    
            for t in rtokens:
                name = t.split('/')[-1]
                name = name.split('_')
                coords = name[1].split('-')
                rpids.append(int(name[0].split("token")[-1]))
                ln1r.append(int(coords[0]))
                ln2r.append(int(coords[1]))
                lt1r.append(int(coords[2]))
                lt2r.append(int(coords[3]))
        
        for jlat in range(lat1,lat2):
            for jlon in range(lon1,lon2):
                for nang in range(0,4):
                    for vw in ['Z','N','E','S','W']:
                        sfile = outputdir+"/sbout.%02d_%02d_%1d_%s"%(jlat,jlon,nang,vw)
                        if sfile not in spectra:
                            n=nearest(lt1r,jlat)
                            if ln2d[n]<jlon:
                                n+=1
                            fixpid = rpids[n]
                            if fixpid in dpids:
                                addtopriority(top,outputdir,jlat,jlon)
                            complete = False
                
    
    if complete:
        taskf = open(top+"/priority.crwl","r")
        tasks = taskf.read().split('\n')
        taskf.close()
        
        last = tasks[-1]
        n = -1
        while last=='':
            n-=1
            last = tasks[n]
            if n==-len(tasks):
                last = '0 DUMMY'
                break
        
        name = outputdir.split('/')[-1]
        os.system("echo "+name+">>"+top+"/inuse.log")
        lastpid = int(last.split()[0])
        head = "# PID MODEL JOBNAME STATUS NCORES QUEUE type gcm workdir lon0"
        jobt = str(lastpid-1)+" postprocess "+name+" 0 1 workq "+model+" "+gcm+" "+outputdir+' '+lon0
        
        if color:
            head+=" color"
            jobt+=" True"
        if makemap:
            head+=" map"
            jobt+=" True"
        
        usef = open(top+"/inuse.crwl","r")
        inuse = int(usef.read()[0])
        usef.close()
        if inuse==0:
            os.system("echo '1'>"+top+"/inuse.crwl")
            os.system('echo "'+head+'">>'+top+"/priority.crwl")
            os.system('echo "'+jobt+'">>'+top+"/priority.crwl")
            f=open(top+"/inuse.crwl","w")
            f.write("0")
            f.close()
            
            
    #Clean up
    os.system("cat "+top+"/inuse.crwl>>"+top+"/inuse.log")
    os.system("python release.py")
