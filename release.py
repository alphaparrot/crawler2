import os
import time
import numpy as np
from crawldefs import Job
import glob
import sys
from sets import Set

'''
With this new system, rather than each routine trying to edit a file, they just
create files for themselves in the 'waitlist' folder. They check to see if the 
main program is currently in use, and if it's not, then they go ahead, seize
control, clean up *everyone* currently in the waitlist, and then run the main
program. This should *in theory* reduce collisions. If anything, there is a risk
that occasionally the program might stall out, since nobody calls the main program.
'''

def joinwaitlist():
    job = np.load("job.npy").item()
    name = job.name
    pid = job.pid
    os.system("cp job.npy "+job.top+"/waitlist/job_"+name+"_"+pid+".npy")
    
def sweep(top):
    usef = open(top+"/inuse.crwl","r")
    inuse = int(usef.read()[0])
    usef.close()
    waitlist = glob.glob(top+"/waitlist/*.npy")
    if inuse==0: #If it's in use, we do nothing and just pass away
        name = np.load("job.npy").item().name
        os.system("echo 'Job "+name+" assuming control'>>"+top+"/inuse.log")
        os.system("echo '1'>"+top+"/inuse.crwl")
        for jobf in waitlist:
            print(jobf)
            job = np.load(jobf).item()
            sig = job.name
            jid = job.home
            pid = job.pid
            args = job.args
            fields = job.fields
            params = job.parameters
            model = job.model
            
            #runf = open(top+"/running_"+model+".crwl","r")
            #running = runf.read().split('\n')[0]
            #runf.close()
            #running = running.split()
            #running[job.home-1] = '0'
            #running = ' '.join(running)+'\n'
            #runf = open(top+"/running_"+model+".crwl","w")
            #runf.write(running)
            #runf.close()
            
            if int(pid)<0:
                tasksf=open(top+"/priority.crwl","r")
                tasks = tasksf.read().split('\n')
                tasksf.close()
                for i in range(0,len(tasks)-1):
                    if tasks[i]!='':
                        tasks[i] = tasks[i].split()
                        if tasks[i][0]!='#':
                            if int(tasks[i][0])==int(pid):
                                tasks[i][3]="2"
                                tasks[i] = ' '.join(tasks[i])
                                break
                        tasks[i] = ' '.join(tasks[i])
                tasks = '\n'.join(tasks)
                tasksf=open(top+"/priority.crwl","w")
                tasksf.write(tasks)
                tasksf.close()
                
            else:
                os.system("cp "+top+"/tasks.crwl "+top+"/tasks.crwl_bkup")
                tasksf=open(top+"/tasks.crwl","r")
                tasks = tasksf.read().split('\n')
                tasksf.close()
                for i in range(0,len(tasks)-1):
                    if tasks[i]!='':
                        tasks[i] = tasks[i].split()
                        if tasks[i][0]!='#':
                            if int(tasks[i][0])==int(pid):
                                tasks[i][3]="2"
                                tasks[i] = ' '.join(tasks[i])
                                break
                        tasks[i] = ' '.join(tasks[i])
                tasks = '\n'.join(tasks)
                if Set(tasks).issubset(Set(" \n")):
                    os.system("echo 'Something is horribly wrong!'>>"+top+"/crashlog.crwl")
                else:
                    tasksf=open(top+"/tasks.crwl","w")
                    tasksf.write(tasks)
                    tasksf.close()
                time.sleep(1.0)
                f=open(top+"/tasks.crwl","r")
                test=f.read()
                f.close()
                if Set(test).issubset(Set(" \n")): #File is somehow empty
                    os.system("cp "+top+"/tasks.crwl_bkup "+top+"/tasks.crwl")
                
            
            os.system("rm "+jobf)
        
        os.system("echo 'Job "+name+" attempting to run crawler2'>>"+top+"/inuse.log")
        os.system("cd "+top+" && python crawler2.py")    
            

if __name__=="__main__":
    if "MASTER" in sys.argv[:]:
        sweep(os.getcwd()) #ONLY CALL FROM TOP-LEVEL DIRECTORY
    else:
        joinwaitlist()
        job = np.load("job.npy").item()
        #os.system("cat "+job.top+"/inuse.crwl>>"+job.top+"/inuse.log")
        sweep(job.top)
        