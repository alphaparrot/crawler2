import numpy as np
import netCDF4 as nc
import os
import time
from batch_system import SUB, BATCHSCRIPT


def edit_data(datafile,namelist,keyword,value,form='scalar'):
    '''form can be scalar or list'''
    try:
        f=open(datafile,"r")
        datat = f.read().split('\n')
        f.close()
        found=False
        for l in range(0,len(datat)):
            if (keyword in datat[l]) and (datat[l][0]!='#'):
                if form=='scalar':
                    datat[l] = ' '+keyword+'='+str(value)+','
                elif form=='list':
                    datat[l] = ' '+keyword+'='+','.join(np.array(value).astype(str))+','
                found=True
                break
                
        if not found:
            for l in range(0,len(datat)):
                if (('&'+namelist) in datat[l]) or (namelist in datat[l]):
                    if form=='scalar':
                        datat.insert(l+1,' '+keyword+'='+str(value)+',')
                    elif form=='list':
                        datat.insert(l+1,' '+keyword+'='+','.join(np.array(value).astype(str))+',')
                    found=True
                    break
        datat = '\n'.join(datat)
        f=open(datafile,"w")
        f.write(datat)
        f.close()
        return found
    except:
        return False

    
def prep(job): 
  workdir = job.top+"/mitgcm/job"+str(job.home)
  if "source" in job.parameters:
    source = job.parameters["source"]
  else:
    source = "locked"
    
  NLEV=10
  if "levs" in job.parameters:
      NLEV = int(job.parameters["levs"])
  
  source+="/l%02d"%NLEV
  
  os.system("cp "+source+"/* "+workdir+"/")
  
  if "p0" in job.parameters: #in Pa
    p0 = float(job.parameters["p0"])
    if "linearpressure" in job.parameters:
        if int(job.parameters['linearpressure'])==1:
            dps = np.zeros(NLEV)+p0*0.1
        else:
            if "ptop" in job.parameters:
                ptop = float(job.parameters["ptop"])
            else:
                ptop = 10.0
            dps = -np.diff(np.geomspace(ptop,p0,num=NLEV+1)[::-1])
    if not edit_data(workdir+"/data","PARM01","atm_Po",p0):
        print "WARNING: Unable to set atm_Po in data/&PARM01!"
    if not edit_data(workdir+"/data","PARM04","delR",dps,format='list'):
        print "WARNING: Unable to set delR in data/&PARM04!"
    
  if "flux" in job.parameters: #in W/m^2
    solc = float(job.parameters["flux"])/1367.0 * 342.0 #Normalize to MITgcm standards
    if not edit_data(workdir+"/data.aimphys","AIM_PAR_FOR","SOLC",solc):
        print "WARNING: Unable to set SOLC in data.aimphys/&AIM_PAR_FOR!"
        
  if "obliq" in job.parameters: #in degrees
    obliq = float(job.parameters["obliq"])
    if not edit_data(workdir+"/data.aimphys","AIM_PAR_FOR","OBLIQ",obliq):
        print "WARNING: Unable to set OBLIQ in data.aimphys/&AIM_PAR_FOR!"
    
  if "pCO2" in job.parameters: #in bars
    pCO2 = float(job.parameters["pCO2"])
    if not edit_data(workdir+"/data.aimphys","AIM_PARAMS","aim_select_pCO2",1):
        print "WARNING: Unable to set CO2 type in data.aimphys/&AIM_PARAMS!"
    if not edit_data(workdir+"/data.aimphys","AIM_PARAMS","aim_fixed_pCO2",pCO2):
        print "WARNING: Unable to set pCO2 in data.aimphys/&AIM_PARAMS!"
    
  if "grav" in job.parameters: #in m/s^2
    grav = float(job.parameters["grav"])
    if not edit_data(workdir+"/data","PARM01","gravity",grav):
        print "WARNING: Unable to set gravity in data/&PARM01!"
    
  if "radius" in job.parameters: #in Rearth
    radius = float(job.parameters["radius"])*6370.0e3
    if not edit_data(workdir+"/data","PARM04","radius_fromHorizGrid",radius):
        print "WARNING: Unable to set radius_fromHorizGrid in data/&PARM04!"

  if "timestep" in job.parameters: #in seconds
    dt = float(job.parameters["timestep"])
    if not edit_data(workdir+"/data","PARM03","deltaT",dt):
        print "WARNING: Unable to set deltaT in data/&PARM03!"
  
  if "runtime" in job.parameters: #in days
    tt = float(job.parameters["runtime"])*86400.0 #convert to seconds
    if not edit_data(workdir+"/data","PARM03","endTime",tt):
        print "WARNING: Unable to set endTime in data/&PARM03!"
        
  
  if "rotation" in job.parameters: #in days
    rt = float(job.parameters["rotation"])*86400.0 #convert to seconds
    if not edit_data(workdir+"/data","PARM01","rotationPeriod",rt):
        print "WARNING: Unable to set rotationPeriod in data/&PARM01!"
    
  notify = 'a'
  if "notify" in job.parameters:
      notify = job.parameters["notify"]

  jobscript =(BATCHSCRIPT(job,notify)+
              "mkdir /mnt/node_scratch/paradise/mitgcm                              \n"+
              "mkdir /mnt/node_scratch/paradise/mitgcm/job"+str(job.home)+"         \n"+
              "rm /mnt/node_scratch/paradise/mitgcm/job"+str(job.home)+"/*          \n"+
              "cp -a * /mnt/node_scratch/paradise/mitgcm/job"+str(job.home)+"/      \n"+
              "cd /mnt/node_scratch/paradise/mitgcm/job"+str(job.home)+"            \n"+
              "                                                                     \n"+
              "module load intel/intel-18                                           \n"+
              "module load openmpi/3.0.0-intel-18                                   \n"+
              "load_netcdf_intel                                                    \n"+
              "                                                                     \n"+
              "mpirun -np 6 ./mitgcmuv > output.txt                                 \n"+
              "python gluemnc.py                                                    \n"+
              "rm *.t*.nc                                                           \n"+
              "mkdir output && mv *.nc output/                                      \n"+
              "                                                                     \n"+
              "cp -a * $PBS_O_WORKDIR/                                              \n"+
              "rm -rf *                                                             \n"+
              "cd $PBS_O_WORKDIR                                                    \n"+                  
              './release.sh "'+job.name+'"                                \n')
  

  rs = open(workdir+"/runmit","w")
  rs.write(jobscript)
  rs.close()
  
  os.system("cp "+job.top+"/crawldefs.py "+workdir+"/")
  os.system("cp "+job.top+"/identity.py "+workdir+"/")
  os.system("cp "+job.top+"/release.py "+workdir+"/")
  

def submit(job):
  workdir = job.top+"/mitgcm/job"+str(job.home)
  
  os.system("cd "+workdir+" && "+SUB+" runmit && cd "+job.top)