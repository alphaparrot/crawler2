import os
import time

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
    
    while waitl[1].split()[0]!=(sig+"@"+pid):
        time.sleep(0.1)
        waitf = open("../../waitlist.crwl","r")
        waitl = waitf.read().split('\n')
        waitf.close()
        
    waitf = open("../../waitlist.crwl","rw")
    waitl = waitf.read().split('\n')
    
    waitl.pop(1)
    
    os.system("echo 1>../../inuse.crwl")
    
    runf = open("../../running_"+model+".crwl","rw")
    running = runf.read().split('\n')
    running[0] = running[0].split()
    running[0][job.home] = '0'
    running[0] = ' '.join(running[0])
    running = '\n'.join(running)
    runf.write(running)
    runf.close()
    
    tasksf=open("../../tasks.crwl","rw")
    tasks = tasksf.read().split('\n')
    for i in range(0,len(tasks)):
        tasks[i] = tasks[i].split()
        if tasks[i][0]!='#':
            if int(tasks[i][0])==pid:
                tasks[i][3]="2"
                tasks[i] = ' '.join(tasks[i])
                break
        tasks[i] = ' '.join(tasks[i])
    tasks = '\n'.join(tasks)
    tasksf.write(tasks)
    tasksf.close()
        
    os.system("cd ../../ && python crawler2.py")
    
    waitf.write('\n'.join(waitl))
    waitf.close()
    
    