import numpy as np
import os
import sys
import netCDF4 as nc
import glob
import time

gplasim = True
TIMELIMIT = 1.44e5

#This version lets the model relax, ignores initial restart files, and freezes the model before
#commencing the experiment.

def edit_namelist(filename,arg,val): 
  f=open(filename,"r")
  fnl=f.read().split('\n')
  f.close()
  found=False
  for l in range(1,len(fnl)-2):
    fnl[l]=fnl[l].split(' ')
    if '=' in fnl[l]:
      mode='EQ'
    else:
      mode='CM'
    if arg in fnl[l]:
      fnl[l]=['',arg,'','=','',str(val),'']
      found=True
    elif (arg+'=') in fnl[l]:
      fnl[l]=['',arg+'=','',str(val),'',',']
      found=True
    fnl[l]=' '.join(fnl[l])
  if not found:
    if mode=='EQ':
      fnl.insert(-3,' '+arg+' = '+val+' ')
    else:
      fnl.insert(-3,' '+arg+'= '+val+' ,')
  f=open(filename,"w")
  f.write('\n'.join(fnl))
  f.close() 

def spatialmath(lt,ln,variable,mean=True,radius=6.371e6):
    lt1 = np.zeros(len(lt)+1)
    lt1[0] = 90
    for n in range(0,len(lt)-1):
        lt1[n+1] = 0.5*(lt[n]+lt[n+1])
    lt1[-1] = -90
    ln1 = np.zeros(len(ln)+1)
    ln1[0] = -2.8125
    for n in range(0,len(ln)-1):
        ln1[n+1] = 0.5*(ln[n]+ln[n+1])
    ln1[-1] = 360.0-2.8125
    
    lt1*=np.pi/180.0
    ln1*=np.pi/180.0
    
    darea = np.zeros((len(lt),len(ln)))
    for jlat in range(0,len(lt)):
        for jlon in range(0,len(ln)):
            dln = ln1[jlon+1]-ln1[jlon]
            darea[jlat,jlon] = (np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*dln
    
    svar = variable*darea
    if mean:
        outvar = np.sum(svar)/np.sum(darea)
    else:
        outvar = np.sum(svar) * radius**2
    
    return outvar

def isflat(key="ts",mean=True,radius=6.371e6,baseline=13,threshhold=0.05):
    #Key is the netCDF variable to evaluate, mean toggles whether to track the average or total,
    #radius is the planet radius in meters, and baseline is the number of years over which to measure
    #slope. Default is to track surface temperature. Threshhold is the maximum slope we'll allow.
  files = sorted(glob.glob("*.nc"))
  nfiles = len(files)
  prior=False
  if len(glob.glob("thistory.ps*"))>0:
      thistory = np.loadtxt("thistory.pso")
      nfiles += len(thistory)
      prior=True
  dd = np.zeros(nfiles)
  nstart=0
  if prior:
      dd[:len(thistory)] = thistory[:]
      nstart=len(thistory)
  if len(files) < baseline+2:
    return False
  else:
    for n in range(0,len(files)):
        ncd = nc.Dataset(files[n],"r")
        variable = ncd.variables[key][:]
        if len(variable.shape)>3:
            variable = variable[:,-1,:,:]
        for m in range(0,variable.shape[0]):
            dd[n+nstart] += spatialmath(ncd.variables['lat'][:],ncd.variables['lon'][:],variable[m,:,:],
                                 mean=mean,radius=radius)
        dd[n+nstart] /= variable.shape[0] #Monthly mean
        ncd.close()
    n=len(dd)-3
    tt=np.arange(baseline)+1
    linfits=[]
    for n in range(len(dd)-3,len(dd)):
      sample=dd[n-(baseline-1):n+1]
      linfit=np.polyfit(tt,sample,1)[0]
      linfits.append(abs(linfit))
      
    avglinfit = (linfits[-3]+linfits[-2]+linfits[-1])/3.0
    if avglinfit <= 0.05:
      return np.mean(dd[-(baseline-2):])
    else:
      return False
  
def gethistory(key="ts",mean=True,radius=6.371e6):
    files = sorted(glob.glob("*.nc"))
    dd=np.zeros(len(files))
    for n in range(0,len(files)):
        ncd = nc.Dataset(files[n],"r")
        variable = ncd.variables[key][:]
        if len(variable.shape)>3:
            variable = variable[:,-1,:,:]
        for m in range(0,variable.shape[0]):
            dd[n] += spatialmath(ncd.variables['lat'][:],ncd.variables['lon'][:],variable[m,:,:],
                                 mean=mean,radius=radius)
        dd[n] /= variable.shape[0] #Monthly mean
        ncd.close()
    return dd
  
def getsurftemp():
    files = sorted(glob.glob("*.nc"))
    if len(files)==0:
        return 1.0e10
    nfile = files[-1]
    ncd = nc.Dataset(nfile,"r")
    ts = np.mean(ncd.variables['ts'][:],axis=0) #Annual average
    lt = ncd.variables['lat'][:]
    ln = ncd.variables['lon'][:]
    tavg = spatialmath(lt,ln,ts)
    ncd.close()
    return tavg
  
def hasnans():
    files = sorted(glob.glob("*.nc"))
    ncd = nc.Dataset(files[-1],"r") #Should be most recent
    if np.isnan(np.amax(ncd.variables['ts'][-1,:])):
        return True
    return False

def getsol():
    f=open("planet_namelist","r")
    nml = f.read().split('\n')
    for k in nml:
        if k.split('=')[0].split()[0]=='GSOL0':
            sol = float(k.split('=')[1].split()[0])
            break
    return sol

if __name__=="__main__":
  EXP="MOST"
  NCPU=int(sys.argv[1])
  nlevs=int(sys.argv[2])
  #os.system("rm -f *.nc") #Clean up after old runs
  #os.system("rm -f plasim_restart") #Uncomment for a fresh run when you haven't cleaned up beforehand
  os.system("rm -f Abort_Message")
  tstart = time.clock()
  exfiles = glob.glob("*DIAG*")
  cyear = len(exfiles)
  
  maxyears=cyear+100
  stoprun=False
  
  if cyear==0:
    if gplasim:
        wf=open("weathering.pso","w")
        wf.write("     CO2       AVG SURF T   WEATHERING    OUTGASSING      DpCO2       NEW CO2\n")
        wf.close()
  
    os.system("echo 'New Cooldown'>cooldown.log")
    
    os.system("cp planet_namelist planet_namelist_wait")
    
    sol0 = getsol()
    sol = sol0-25.0
    
    edit_namelist("planet_namelist","GSOL0","%4.1f"%(sol)) #Turn the Sun dim enough to freeze over
    yearsfrozen = 0
    
    while getsurftemp()>255.0 or yearsfrozen<30:
        for nspool in range(0,3):
            os.system("echo 'cooldown year "+str(cyear)+"'>>cooldown.log")
            cyear += 1
            dataname=EXP+".%04d"%cyear
            snapname=EXP+"_SNAP.%04d"%cyear
            diagname=EXP+"_DIAG.%04d"%cyear
            restname=EXP+"_REST.%03d"%cyear
            snowname=EXP+"_SNOW_%1d"%(cyear%5)
            os.system("mpiexec -np "+str(NCPU)+" most_plasim_t21_l%d_p"%nlevs+str(NCPU)+".x")
            os.system("[ -e restart_dsnow ] && rm restart_dsnow")
            os.system("[ -e restart_xsnow ] && rm restart_xsnow")
            os.system("[ -e Abort_Message ] && exit 1")
            os.system("[ -e plasim_output ] && mv plasim_output "+dataname)
            os.system("[ -e plasim_snapshot ] && mv plasim_snapshot "+snapname)
            os.system("[ -e plasim_diag ] && mv plasim_diag "+diagname)
            os.system("[ -e plasim_status ] && cp plasim_status plasim_restart")
            os.system("[ -e plasim_status ] && mv plasim_status "+restname)
            os.system("[ -e restart_snow ] && mv restart_snow "+snowname)
            os.system("[ -e "+dataname+" ] && ./burn7.x -n <example.nl>burnout "+dataname+" "+dataname+".nc")
            os.system("[ -e "+snapname+" ] && ./burn7.x -n <snapshot.nl>snapout "+snapname+" "+snapname+".nc")
            os.system("[ -e "+dataname+" ] && cp "+dataname+" "+EXP+"_OUT.%04d"%cyear)
            os.system("[ -e "+dataname+".nc ] && rm "+dataname)
            os.system("[ -e "+snapname+".nc ] && rm "+snapname)
            os.system("[ -e "+snapname+".nc ] && mv "+snapname+".nc snapshots/")
            os.system("cp weathering.pso $PBS_O_WORKDIR/")
            os.system("cp "+diagname+" $PBS_O_WORKDIR/")
            if hasnans():
                break
            if getsurftemp()<=230.0:
                yearsfrozen+=1
        sol-=25.0
        edit_namelist("planet_namelist","GSOL0","%4.1f"%(sol)) #Turn the Sun dim enough to freeze over
        
    sol += 25.0
    
    while sol<sol0 and not hasnans():   #Ramp back up to target insolation.
        edit_namelist("planet_namelist","GSOL0",str(sol))
        for nt in range(0,3):
            cyear += 1
            dataname=EXP+".%04d"%cyear
            snapname=EXP+"_SNAP.%04d"%cyear
            diagname=EXP+"_DIAG.%04d"%cyear
            restname=EXP+"_REST.%03d"%cyear
            snowname=EXP+"_SNOW_%1d"%(cyear%5)
            os.system("mpiexec -np "+str(NCPU)+" most_plasim_t21_l%d_p"%nlevs+str(NCPU)+".x")
            os.system("[ -e restart_dsnow ] && rm restart_dsnow")
            os.system("[ -e restart_xsnow ] && rm restart_xsnow")
            os.system("[ -e Abort_Message ] && exit 1")
            os.system("[ -e plasim_output ] && mv plasim_output "+dataname)
            os.system("[ -e plasim_snapshot ] && mv plasim_snapshot "+snapname)
            os.system("[ -e plasim_diag ] && mv plasim_diag "+diagname)
            os.system("[ -e plasim_status ] && cp plasim_status plasim_restart")
            os.system("[ -e plasim_status ] && mv plasim_status "+restname)
            os.system("[ -e restart_snow ] && mv restart_snow "+snowname)
            os.system("[ -e "+dataname+" ] && ./burn7.x -n <example.nl>burnout "+dataname+" "+dataname+".nc")
            os.system("[ -e "+snapname+" ] && ./burn7.x -n <snapshot.nl>snapout "+snapname+" "+snapname+".nc")
            os.system("[ -e "+dataname+" ] && cp "+dataname+" "+EXP+"_OUT.%04d"%cyear)
            os.system("[ -e "+dataname+".nc ] && rm "+dataname)
            os.system("[ -e "+snapname+".nc ] && rm "+snapname)
            os.system("[ -e "+snapname+".nc ] && mv "+snapname+".nc snapshots/")
            os.system("cp weathering.pso $PBS_O_WORKDIR/")
            os.system("cp "+diagname+" $PBS_O_WORKDIR/")
            if hasnans():
                break
        sol+=25.0
    
    
    os.system("cp planet_namelist_wait planet_namelist") #Turn the Sun back up
    if abs(cyear-maxyears)<10 or (time.clock()-tstart)>=TIMELIMIT:
        stoprun=True
        os.system("touch keepgoing")
        temps = gethistory()
        with open("thistory.pso","a+") as f:
            text = '\n'+'\n'.join(temps.astype(str))
            f.write(text)
  
  if not stoprun:
    year=0
    minyears=25
    relaxed=False
    while (year < minyears or not relaxed) and (cyear < maxyears) and (time.clock()-tstart)<=TIMELIMIT:
        os.system("echo 'relaxation year "+str(year)+"'>>cooldown.log")
        year+=1
        cyear += 1
        dataname=EXP+".%04d"%cyear
        snapname=EXP+"_SNAP.%04d"%cyear
        diagname=EXP+"_DIAG.%04d"%cyear
        restname=EXP+"_REST.%03d"%cyear
        snowname=EXP+"_SNOW_%1d"%(cyear%5)
        os.system("mpiexec -np "+str(NCPU)+" most_plasim_t21_l%d_p"%nlevs+str(NCPU)+".x")
        os.system("[ -e restart_dsnow ] && rm restart_dsnow")
        os.system("[ -e restart_xsnow ] && rm restart_xsnow")
        os.system("[ -e Abort_Message ] && exit 1")
        os.system("[ -e plasim_output ] && mv plasim_output "+dataname)
        os.system("[ -e plasim_snapshot ] && mv plasim_snapshot "+snapname)
        os.system("[ -e plasim_diag ] && mv plasim_diag "+diagname)
        os.system("[ -e plasim_status ] && cp plasim_status plasim_restart")
        os.system("[ -e plasim_status ] && mv plasim_status "+restname)
        os.system("[ -e restart_snow ] && mv restart_snow "+snowname)
        os.system("[ -e "+dataname+" ] && ./burn7.x -n <example.nl>burnout "+dataname+" "+dataname+".nc")
        os.system("[ -e "+snapname+" ] && ./burn7.x -n <snapshot.nl>snapout "+snapname+" "+snapname+".nc")
        os.system("[ -e "+dataname+" ] && cp "+dataname+" "+EXP+"_OUT.%04d"%cyear)
        os.system("[ -e "+dataname+".nc ] && rm "+dataname)
        os.system("[ -e "+snapname+".nc ] && rm "+snapname)
        os.system("[ -e "+snapname+".nc ] && mv "+snapname+".nc snapshots/")
        os.system("cp weathering.pso $PBS_O_WORKDIR/")
        os.system("cp "+diagname+" $PBS_O_WORKDIR/")
        if hasnans():
            break
        relaxed=isflat(baseline=minyears)
        
    if not relaxed:
        os.system("touch keepgoing")
        temps = gethistory()
        with open("thistory.pso","a+") as f:
            text = '\n'+'\n'.join(temps.astype(str))
            f.write(text)
    
