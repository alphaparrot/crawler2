import numpy as np
import os
import sys
import netCDF4 as nc
import glob
import struct

gplasim = True

#This version lets the model relax, ignores initial restart files, and freezes the model before
#commencing the experiment, then allowing glaciers to adjust.


def getndco2(weatheringfile):
  f=open(weatheringfile,'r')
  r=f.read().split('\n')[1:-1]
  f.close()
  l1=r[1].split()
  l2=r[-1].split()
  sco2=float(l1[0])*1e6
  tco2=float(l2[5])*1e6
  return (sco2,tco2)

def getdco2(weatheringfile):
  f=open(weatheringfile,'r')
  r=f.read().split('\n')[1:-1]
  f.close()
  l=r[-1].split()
  dco2=float(l[4])*1e6
  return dco2

def getweathering(weatheringfile):
  f=open(weatheringfile,'r')
  r=f.read().split('\n')[1:-1]
  f.close()
  l=r[-1].split()
  ww=float(l[2])
  return ww

def changeCO2(pCO2):
  nl=open("radmod_namelist","r")
  nltxt=nl.read().split('\n')
  nl.close()
  l=0
  while l<len(nltxt):
    line=nltxt[l].split()
    if line[0]=='CO2':
      line[2]=str(pCO2)
      nltxt[l]=' '+' '.join(line)
      break
    elif line[0]=='CO2=':
      line[1]=str(pCO2)
      nltxt[l]=' '+' '.join(line)
      break
    l+=1
  nltxt='\n'.join(nltxt)
  nl=open('radmod_namelist','w')
  nl.write(nltxt)
  nl.close()
  
def changep(psurf):
  nl=open("plasim_namelist","r")
  nltxt=nl.read().split('\n')
  nl.close()
  l=0
  while l<len(nltxt)-1:
    line=nltxt[l].split()
    if line[0]=='PSURF':
      line[2]=str(psurf)
      nltxt[l]=' '+' '.join(line)
      break
    elif line[0]=='PSURF=':
      line[1]=str(psurf)
      nltxt[l]=' '+' '.join(line)
      break
    l+=1
  nltxt='\n'.join(nltxt)
  nl=open('plasim_namelist','w')
  nl.write(nltxt)
  nl.close()  

def getmaxdsnow(filename1,filename2):
  f1=open(filename1,'rb')
  r1=f1.read()
  f1.close()
  
  f2=open(filename2,'rb')
  r2=f2.read()
  f2.close()
  try:
    dd1=np.array(struct.unpack('2048d',r1[4:-4])).reshape((32,64))
    dd2=np.array(struct.unpack('2048d',r2[4:-4])).reshape((32,64))
  except:
    dd1=np.array(struct.unpack('2048f',r1[4:-4])).reshape((32,64))
    dd2=np.array(struct.unpack('2048f',r2[4:-4])).reshape((32,64))
    
  
  dsnow = dd2-dd1
  maxdsnow = np.amax(np.abs(dsnow))
  
  return maxdsnow

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

def isflat(key="ts",mean=True,radius=6.371e6,baseline=13,threshhold=0.05,prefix="*"):
    #Key is the netCDF variable to evaluate, mean toggles whether to track the average or total,
    #radius is the planet radius in meters, and baseline is the number of years over which to measure
    #slope. Default is to track surface temperature. Threshhold is the maximum slope we'll allow.
  files = sorted(glob.glob(prefix+"*.nc"))
  if len(files) < baseline+2:
    return False
  else:
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
  
def getsurftemp(prefix=''):
    files = sorted(glob.glob(prefix+"*.nc"))
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
  nlevs = int(sys.argv[2])
  #os.system("rm -f *.nc") #Clean up after old runs
  #os.system("rm -f plasim_restart") #Uncomment for a fresh run when you haven't cleaned up beforehand
  os.system("rm -f Abort_Message")
  
  
  os.system("cp planet_namelist planet_namelist_wait")
  
  sol0 = getsol()
  
  p0 = 1010670.0
  
  stoprun=False
  
  cyear = len(glob.glob("*MOST.*.nc"))
  cmaxyears = cyear+2000
  ntgl = len(glob.glob("*OUTGLAC*.nc"))
  if cyear==0:
    os.system("echo 'New Cooldown'>cooldown.log")
    
    edit_namelist("planet_namelist","GSOL0","400.0") #Turn the Sun dim enough to freeze over
    yearsfrozen = 0  
    if gplasim:
        wf=open("weathering.pso","w")
        wf.write("     CO2       AVG SURF T   WEATHERING    OUTGASSING      DpCO2       NEW CO2\n")
        wf.close()
    while getsurftemp()>255.0 or yearsfrozen<30:
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
        
    sol = 500.0
    while sol<sol0 and not hasnans():   #Ramp back up to target insolation.
        edit_namelist("planet_namelist","GSOL0",str(sol))
        for nt in range(0,2):
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
        sol+=100.0
    if abs(cyear-cmaxyears)<10:
        stoprun=True
        os.system("touch keepgoing")
    
  if not hasnans() and not stoprun:  
    
    os.system("cp planet_namelist_wait planet_namelist") #Turn the Sun back up
    
    pco2 = getndco2('weathering.pso')[1]
    pco2_old = pco2
    
    psurf = p0 + pco2
    
    os.system("cp restart_dsnow restart_dsnow_old")
    os.system("cp restart_xsnow restart_xsnow_old")
    
    os.system("cp newdsnow_old newdsnow")
    os.system("cp newxsnow_old newxsnow")
    os.system("cp plasim_restart_old plasim_restart")
    
    
  #  ntgl = 0
    glacrelaxed = False
    trelaxed = False
    longtime=False
    while ntgl < ntgl+90 and not glacrelaxed and not trelaxed and cyear<cmaxyears: #No more than 50 iterations
      os.system("rm "+EXP+"_OUT.*")
      os.system("cp -a "+EXP+"*.nc $PBS_O_WORKDIR/")
      os.system("rm "+EXP+"*.nc")
      os.system("cp -a "+EXP+"_REST.* $PBS_O_WORKDIR/")
      os.system("rm "+EXP+"_REST.*")
      os.system("cp -a "+EXP+"_OUTGLAC.*.nc $PBS_O_WORKDIR/")
      os.system("cp -a "+EXP+"_SNOW.* $PBS_O_WORKDIR/")
      os.system("rm "+EXP+"_SNOW.*")
      year=0
      minyears=25
      maxyears = 125
      relaxed=False
      while (year < minyears or not relaxed) and (year < maxyears) and cyear<cmaxyears:
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
        os.system("cp cyclelog.txt $PBS_O_WORKDIR/")
        if hasnans():
            break
        relaxed=isflat(baseline=minyears,prefix=EXP+".*")
      os.system("cp "+dataname+".nc "+EXP+"_OUTGLAC.%04d.nc"%ntgl)
      os.system("cp "+snowname+" "+EXP+"_SNOW.%04d"%ntgl)
      sfile1 = EXP+"_SNOW_0"
      sfile2 = EXP+"_SNOW_1"
      sfile3 = EXP+"_SNOW_2"
      sfile4 = EXP+"_SNOW_3"
      sfile5 = EXP+"_SNOW_4"
      
      ts = getsurftemp(prefix=EXP+".")
      if ts > 275: #Undo the last changes
        os.system("cp restart_dsnow_old restart_dsnow")
        os.system("cp restart_xsnow_old restart_xsnow")
        os.system("cp newdsnow_old newdsnow")
        os.system("cp newxsnow_old newxsnow")
        os.system("cp plasim_restart_old plasim_restart")
        pco2 = 0.5*(pco2+pco2_old) #bisection
        psurf = p0 + pco2
        changeCO2(pco2/psurf*1e6)
        changep(psurf*0.1)
        f=open('cyclelog.txt','a')
        f.write("Reverted.\n")
        f.close()
      else:
        
        if getweathering('weathering.pso')>1.0e-6:
            break
        
        os.system("cp restart_dsnow restart_dsnow_old")
        os.system("cp restart_xsnow restart_xsnow_old")
        os.system("cp newdsnow newdsnow_old")
        os.system("cp newxsnow newxsnow_old")
        os.system("cp plasim_restart plasim_restart_old")
          
        deltat = 500.0
        os.system("./newsnow.x "+sfile1+" "+sfile2+" "+sfile3+" "+sfile4+" "+sfile5+" "+str(deltat))
        maxdsnow = getmaxdsnow("newdsnow_old","newdsnow")
        
        #if maxdsnow>1.0e-4: #Just let the glaciers relax if things will change on 100 kyr timescales.
          #fraction=127.5/maxdsnow
          #deltat*=fraction
          ##deltat = min(deltat,1.0e5) #No more than 100kyr at a time
          ##if deltat==1.0e5:
          ##    longtime=True
          #os.system("cp newdsnow_old newdsnow")
          #os.system("cp newxsnow_old newxsnow")
          #os.system("./newsnow.x "+sfile1+" "+sfile2+" "+sfile3+" "+sfile4+" "+sfile5+" "+str(deltat))
          #f=open('cyclelog.txt','a')
          #maxdsnow = getmaxdsnow("newdsnow_old","newdsnow")
          #f.write("Advanced "+str(deltat)+" years. Maximum snow change was "+str(maxdsnow)+".\n")
          #f.close()
          #os.system("mv newdsnow restart_dsnow")
          #os.system("mv newxsnow restart_xsnow")
        
        #else: #Warm up
          ##longtime=False
          #os.system("cp newdsnow_old newdsnow")
          #os.system("cp newxsnow_old newxsnow")
          #os.system("mv newdsnow restart_dsnow")
          #os.system("mv newxsnow restart_xsnow")
          #pco2_old = pco2
          #pco2 = pco2 * 10**0.05 #20 per dex
          #f=open('cyclelog.txt','a')
          #f.write("Increased CO2 to "+str(pco2*1.0e-6)+" bars.\n")
          #f.close()
          #psurf = p0 + pco2
          #changeCO2(pco2/psurf*1e6)
          #changep(psurf*0.1)
          
        fraction=127.5/maxdsnow
        deltat*=fraction  
        pco2_old = pco2
        pco2 += 0.05*deltat #50 nbars/yr
        if maxdsnow<1.0e-15 or np.isnan(maxdsnow):
            pco2 = pco2_old*10**0.05
            dpco2 = pco2-pco2_old
            deltat = dpco2/(0.05)
        if pco2/pco2_old > pco2_old*10**0.05:
            pco2 = pco2_old*10**0.05
            dpco2 = pco2-pco2_old
            deltat = dpco2/0.05
        #deltat = min(deltat,1.0e5) #No more than 100kyr at a time
        #if deltat==1.0e5:
        #    longtime=True
        os.system("cp newdsnow_old newdsnow")
        os.system("cp newxsnow_old newxsnow")
        os.system("./newsnow.x "+sfile1+" "+sfile2+" "+sfile3+" "+sfile4+" "+sfile5+" "+str(deltat))
        f=open('cyclelog.txt','a')
        maxdsnow = getmaxdsnow("newdsnow_old","newdsnow")
        f.write("Advanced "+str(deltat)+" years. Maximum snow change was "+str(maxdsnow)+".\n")
        f.close()
        os.system("mv newdsnow restart_dsnow")
        os.system("mv newxsnow restart_xsnow")
        f=open('cyclelog.txt','a')
        f.write("Increased CO2 to "+str(pco2*1.0e-6)+" bars.\n")
        f.close()
        psurf = p0 + pco2
        changeCO2(pco2/psurf*1e6)
        changep(psurf*0.1)      
        
        
        glacrelaxed = isflat(baseline=5,key='icez',threshhold=10.0,prefix=EXP+"_OUTGLAC.*") #less than 10m average annual change over 5 cycles
        trelaxed = isflat(baseline=5,prefix=EXP+"_OUTGLAC.*") #Temperature should be flat too
        ntgl+=1
    if not glacrelaxed or not trelaxed:
        os.system("touch keepgoing")
