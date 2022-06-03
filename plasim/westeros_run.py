import numpy as np
import os
import sys
import netCDF4 as nc
import glob
import time

gplasim = True
TIMELIMIT = 1.44e5

#This version lets the model relax.

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

  
def hasnans():
    files = sorted(glob.glob("*.nc"))
    if len(files)<1:
        return False
    print("NetCDF  files:",files)
    if type(files)!=type([1,2,3]):
        files = [files,]
    ncd = nc.Dataset(files[-1],"r") #Should be most recent
    if np.nansum(1.0*np.isnan(ncd.variables['ts'][-1,:]))>0.5:
        return True
    return False


if __name__=="__main__":

  EXP="MOST"
  os.system("rm keepgoing")
  tstart = time.clock()
  NCPU=int(sys.argv[1])
  nlevs = int(sys.argv[2])
  runyears = int(sys.argv[3])
  #os.system("rm -f plasim_restart") #Uncomment for a fresh run when you haven't cleaned up beforehand
  os.system("rm -f Abort_Message")
  os.system("echo 'SURFACE      TOA'>balance.log")
  os.system("echo 'SURFACE      TOA'>slopes.log")
  exfiles = glob.glob("*DIAG*")
  year=len(exfiles)
  if gplasim and year==0:
    os.system("rm sitnikov.pso")
    wf=open("weathering.pso","w")
    wf.write("     CO2       AVG SURF T   WEATHERING    OUTGASSING      DpCO2       NEW CO2\n")
    wf.close()
  minyears=75
  maxyears=year+300
  relaxed=False
  while (year < minyears or year<maxyears) and (time.clock()-tstart)<=TIMELIMIT and year<runyears:
    year+=1
    dataname=EXP+".%04d"%year
    snapname=EXP+"_SNAP.%04d"%year
    diagname=EXP+"_DIAG.%04d"%year
    restname=EXP+"_REST.%03d"%year
    snowname=EXP+"_SNOW_%1d"%(year%5)
    os.system("mpiexec -np "+str(NCPU)+" most_plasim_t21_l"+str(nlevs)+"_p"+str(NCPU)+".x")
    os.system("[ -e restart_dsnow ] && rm restart_dsnow")
    os.system("[ -e restart_xsnow ] && rm restart_xsnow")
    os.system("[ -e Abort_Message ] && exit 1")
    os.system("[ -e plasim_output ] && mv plasim_output "+dataname)
    os.system("[ -e plasim_snapshot ] && mv plasim_snapshot "+snapname)
    os.system("[ -e plasim_diag ] && mv plasim_diag "+diagname)
    os.system("[ -e plasim_status ] && cp plasim_status plasim_restart")
    os.system("[ -e plasim_status ] && mv plasim_status "+restname)
    os.system("[ -e restart_snow ] && mv restart_snow "+snowname)
    os.system("[ -e "+dataname+" ] && ./burn7.x -n <essos.nl>essos.out "+dataname+" "+dataname+".nc")
    os.system("[ -e "+snapname+" ] && ./burn7.x -n <snapshot.nl>snapout "+snapname+" "+snapname+".nc")
    os.system("[ -e "+dataname+" ] && cp "+dataname+" "+EXP+"_OUT.%04d"%year)
    os.system("[ -e "+dataname+".nc ] && rm "+dataname)
    os.system("[ -e "+snapname+".nc ] && rm "+snapname)
    os.system("[ -e "+snapname+".nc ] && mv "+snapname+".nc snapshots/")
    if hasnans():
        os.system("echo 'NAN ENCOUNTERED'>>weathering.pso")
        break
  os.system("rm keepgoing")
  if not hasnans() and year<runyears:
    os.system("touch keepgoing")
        
