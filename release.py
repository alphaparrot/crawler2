import os
import time
import numpy as np
from crawldefs import Job

if __name__=="__main__":
    
    job = np.load("job.npy").item()
    sig = job.name
    jid = job.home
    pid = job.pid
    args = job.args
    fields = job.fields
    params = job.parameters
    model = job.model
    
    os.system("echo "+sig+"@"+pid+">>../../waitlist.crwl")
    
    waitf = open("../../waitlist.crwl","r")
    waitl = waitf.read().split('\n')
    waitf.close()
    
    usef = open("../../inuse.crwl","r")
    inuse = int(usef.read()[0])
    usef.close()
    
    while waitl[1].split()[0]!=(sig+"@"+pid) and inuse==1:
        time.sleep(0.1)
        waitf = open("../../waitlist.crwl","r")
        waitl = waitf.read().split('\n')
        waitf.close()
        usef = open("../../inuse.crwl","r")
        inuse = int(usef.read()[0])
        usef.close()
        
        
    waitf = open("../../waitlist.crwl","r+")
    waitl = waitf.read().split('\n')
    
    waitl.pop(1)
    
    os.system("echo '1'>../../inuse.crwl")
    
    runf = open("../../running_"+model+".crwl","r")
    running = runf.read().split('\n')[0]
    runf.close()
    running = running.split()
    running[job.home-1] = '0'
    running = ' '.join(running[0])+'\n'
    runf = open("../../running_"+model+".crwl","w")
    runf.write(running)
    runf.close()
    
    tasksf=open("../../tasks.crwl","r")
    tasks = tasksf.read().split('\n')
    tasksf.close()
    for i in range(0,len(tasks)-1):
        if tasks[i]!='':
            tasks[i] = tasks[i].split()
            if tasks[i][0]!='#':
                if int(tasks[i][0])==pid:
                    tasks[i][3]="2"
                    tasks[i] = ' '.join(tasks[i])
                    break
            tasks[i] = ' '.join(tasks[i])
    tasks = '\n'.join(tasks)
    tasksf=open("../../tasks.crwl","w")
    tasksf.write(tasks)
    tasksf.close()
        
    os.system("cd ../../ && python crawler2.py")
    
    waitf.write('\n'.join(waitl))
    waitf.close()
    
    