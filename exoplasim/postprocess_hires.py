import os
import glob
import numpy as np
import netCDF4 as nc


#This version lets the model relax.
def get_namelist(filename,arg): 
  f=open(filename,"r")
  fnl=f.read().split('\n')
  f.close()
  val = "-1"
  found=False
  for l in range(1,len(fnl)-2):
    if arg in fnl[l]:
        val = fnl[l].split('=')[1]
        found=True
        break
  return val



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
  
def hasnans():
    files = sorted(glob.glob("*.nc"))
    print("NetCDF  files:",files)
    if type(files)!=type([1,2,3]):
        files = [files,]
    ncd = nc.Dataset(files[-1],"r") #Should be most recent
    if np.sum(1.0*np.isnan(ncd.variables['ts'][-1,:]))>0.5:
        return True
    return False

def energybalanced(cyear,threshhold = 1.0e-4,baseline=50): #Takes an average of 200 years
    if cyear<baseline:
        return False
    files = sorted(glob.glob("*.nc"))
    nfiles = len(files)
    prior=False
    if len(glob.glob("toahistory.ps*"))>0:
        try:
            toahistory = np.loadtxt("toahistory.pso")
            nfiles+=len(toahistory)
            shistory = np.loadtxt("shistory.pso")
            prior=True
        except:
            pass
    sbalance = np.zeros(nfiles)
    toabalance=np.zeros(nfiles)
    nstart=0
    if prior:
        sbalance[:len(toahistory)] = shistory[:]
        toabalance[:len(toahistory)] = toahistory[:]
        nstart = len(toahistory)
    if len(files) < baseline: #Run for minimum of baseline years
        return False
    else:
        for n in range(0,len(files)):
            ncd = nc.Dataset(files[n],"r")
            ntr = ncd.variables['ntr'][:]
            hfns = ncd.variables['hfns'][:]
            lat = ncd.variables['lat'][:]
            lon = ncd.variables['lon'][:]
            ncd.close()
            ntimes = ntr.shape[0]
            topt = np.zeros(ntimes)
            bott = np.zeros(ntimes)
            for m in range(0,ntimes):
                topt[m] = spatialmath(lat,lon,ntr[m,:,:])
                bott[m] = spatialmath(lat,lon,hfns[m,:,:])
            sbalance[n+nstart] = np.mean(bott)
            toabalance[n+nstart] = np.mean(topt)
        savgs = []
        tavgs = []
        for n in range(9,len(sbalance)):
            savgs.append(abs(np.mean(sbalance[n-9:n+1]))) #10-year average energy balance
            tavgs.append(abs(np.mean(toabalance[n-9:n+1])))
        sslopes = []
        tslopes = []
        for n in range(4,len(savgs)): #5-baseline slopes in distance from energy balance
            sslopes.append(np.polyfit(np.arange(5)+1,savgs[n-4:n+1],1)[0])
            tslopes.append(np.polyfit(np.arange(5)+1,tavgs[n-4:n+1],1)[0])
        savgslope = abs(np.mean(sslopes[-30:])) #30-year average of 5-year slopes  
        tavgslope = abs(np.mean(tslopes[-30:]))
        os.system("echo '%02.8f  %02.8f'>>slopes.log"%(savgslope,tavgslope))
        print("%02.8f %02.8f"%(savgslope,tavgslope))
        if savgslope<threshhold and tavgslope<threshhold: #Both TOA and Surface are changing at average 
            return True                                  # of <0.1 mW/m^2/yr on 45-year baselines
        else:
            return False
        
def getbalance():
    files = sorted(glob.glob("*.nc"))
    ncd = nc.Dataset(files[-1],"r")
    ntr = ncd.variables['ntr'][:]
    hfns = ncd.variables['hfns'][:]
    lat = ncd.variables['lat'][:]
    lon = ncd.variables['lon'][:]
    ncd.close()
    ntimes = ntr.shape[0]
    topt = np.zeros(ntimes)
    bott = np.zeros(ntimes)
    for m in range(0,ntimes):
        topt[m] = spatialmath(lat,lon,ntr[m,:,:])
        bott[m] = spatialmath(lat,lon,hfns[m,:,:])
    return (np.mean(bott),np.mean(topt))
       

    
def rpvorticity(name,daylen):
    
    omega = 2*np.pi/(daylen*86400.0)
    
    ncd = nc.Dataset(name,"r")
    ln = ncd.variables['lon'][:]
    lt = ncd.variables['lat'][:]
    lev = ncd.variables['lev'][:]
    ta = ncd.variables['ta'][:]
    ps = ncd.variables['ps'][:]
    
    #ln,lt,ta = gt.parse(name,"ta")
    #ln,lt,lev = gt.parse(name,"lev")
    #ln,lt,ps = gt.parse(name,"ps")
    
    print("Bunch 1 loaded")
    
    nlevs = len(lev)
    ntimes = ta.shape[0]
    nlats = len(lt)
    nlons = len(ln)
    pa = np.zeros(ta.shape)
    #print nlevs
    for k in range(nlevs):
        #print k,type(k),ps.shape,pa.shape
        pa[:,k,:,:] = ps*100*lev[k]
    
    print("pa")
    
    pvap = 6.1078*10**(7.5*(ta-273.15)/(ta-273.15+237.3))
    
    print("psat")
    
    hur = ncd.variables['hur'][:]
    
    pvap = hur*pvap
    del hur
    
    print("pvap")
    
    pdry = pa-pvap
    rho = (pdry*0.0289654+pvap*0.018016)/(8.314*ta)
    del pdry
    del pvap
    
    print("rho",np.nanmin(rho),np.nanmax(rho))
    
    rr = -8.314*ta/9.80665*np.log(pa*0.01/ps[:,np.newaxis,:,:])
    rr = 6371e6+rr
    del ps
    
    print("rr",np.nanmin(rr),np.nanmax(rr))
    
    theta = ta*(1e5/pa)**0.286
    del pa
    del ta
    
    print("theta")
        
    #lt = 90.0-lt #polar angle starts at the north pole
        
    #ln,lt,zeta = gt.parse(name,"zeta")
    ln *= np.pi/180.0
    lt *= np.pi/180.0
    
    
    #ncd = nc.Dataset(name,"r")
    #ua = ncd.variables['ua'][:]
    #va = ncd.variables['va'][:]
    #wa = ncd.variables['wa'][:]
    #ncd.close()
    
    #drvdr = np.zeros(rr.shape)
    #drudr = np.zeros(rr.shape)
    #dthetadr = np.zeros(rr.shape)
    #for t in range(ntimes):
        #print "Time ",t
        #for j in range(nlats):
            #for k in range(nlons):
                #drvdr[t,:,j,k] = np.gradient(rr[t,:,j,k]*va[t,:,j,k], rr[t,:,j,k])
                #drudr[t,:,j,k] = np.gradient(rr[t,:,j,k]*ua[t,:,j,k], rr[t,:,j,k])
                #dthetadr[t,:,j,k] = np.gradient(theta[t,:,j,k],rr[t,:,j,k])
    
    #duslatdlat = np.gradient(np.sin(lt)[np.newaxis,np.newaxis,:,np.newaxis]*ua,lt,axis=2)
    #dvdlon = np.gradient(va,ln,axis=3)
    #dwdlon = np.gradient(wa,ln,axis=3)
    #dthetadlat = np.gradient(theta,lt,axis=2)
    #dthetadlon = np.gradient(theta,ln,axis=3)
    #dwdlat = np.gradient(wa,lt,axis=2)
    
    #rsinlat = rr*(np.sin(lt)[np.newaxis,np.newaxis,:,np.newaxis])
    
    #pvort = (dthetadr/rsinlat*(2*omega*np.cos(lt)[np.newaxis,np.newaxis,:,np.newaxis] +
                                #duslatdlat - dvdlon)
            #+dthetadlat/rr**2*(2*omega*np.sin(lt)[np.newaxis,np.newaxis,:,np.newaxis] +
                               #dwdlon/(np.sin(lt)[np.newaxis,np.newaxis,:,np.newaxis])-
                               #drudr)
            #+dthetadlon/(rr*rsinlat)*(drvdr - dwdlat))/rho
    
    clt = np.cos(lt)
    slt = np.sin(lt)
    
    #pv = 1/rho*((dwdlat/rr - drvdr/rr)*dthetadlon/(rr*clats) +
                #(2*omega*clt + drudr/rr - dwdlon/(rr*clats))*dthetadlat/rr +
                #(2*omega*slt + dvdlon/(rr*clats) - duclatdlat/(rr*clats))*dthetadr)

    wa = ncd.variables['wa'][:]
    dwdlat = np.zeros(wa.shape)
    
    for t in range(ntimes):
        for j in range(nlevs):
            for k in range(nlons):
                dwdlat[t,j,:,k] = np.gradient(wa[t,j,:,k],lt)
                
    dwdlon = np.zeros(wa.shape)
    for t in range(ntimes):
        for j in range(nlevs):
            for k in range(nlats):
                dwdlon[t,j,k,:] = np.gradient(wa[t,j,k,:],ln)
    del wa
    
    print("bunch 2")
    

    va = ncd.variables['va'][:]
    ua = ncd.variables['ua'][:]
    ncd.close()

    
    drvdr = np.zeros(rr.shape)
    drudr = np.zeros(rr.shape)
    dthetadr = np.zeros(rr.shape)
    print("Starting gradients")
    for t in range(ntimes):
        print("Time ",t)
        for j in range(nlats):
            for k in range(nlons):
                drvdr[t,:,j,k] = np.gradient(rr[t,:,j,k]*va[t,:,j,k], rr[t,:,j,k])
                drudr[t,:,j,k] = np.gradient(rr[t,:,j,k]*ua[t,:,j,k], rr[t,:,j,k])
                dthetadr[t,:,j,k] = np.gradient(theta[t,:,j,k],rr[t,:,j,k])

    print("3D Gradients finished")
    
    dvdlon = np.zeros(va.shape)
    for t in range(ntimes):
        for j in range(nlevs):
            for k in range(nlats):
                dvdlon[t,j,k,:] = np.gradient(va[t,j,k,:],ln)
    del va
    
    print("bunch 3")
    
    thing = dwdlat-drvdr
    del dwdlat
    del drvdr
    
    print("bunch 4")
    
    dthetadlon = np.zeros(theta.shape)
    for t in range(ntimes):
        for j in range(nlevs):
            for k in range(nlats):
                dthetadlon[t,j,k,:] = np.gradient(theta[t,j,k,:],ln)
    thing = thing*dthetadlon
    del dthetadlon
    
    print("bunch 5")
    
    clats = clt[np.newaxis,np.newaxis,:,np.newaxis]*np.ones(ua.shape)
    
    thing = thing/(clats*rr**2)
    #thing = 1/rho*((dwdlat-drvdr)/rr)*dthetadlon/(rr*clats)
    
    
    duclatdlat = np.zeros(ua.shape)
    for t in range(ntimes):
        for j in range(nlevs):
            for k in range(nlons):
                duclatdlat[t,j,:,k] = np.gradient(ua[t,j,:,k]*clt,lt)
    del ua
    print("bunch 6")
    
    dthetadlat = np.zeros(theta.shape)
    for t in range(ntimes):
        for j in range(nlevs):
            for k in range(nlons):
                dthetadlat[t,j,:,k] = np.gradient(theta[t,j,:,k],lt)
    del theta
    print("bunch 7")
    
    #thing = thing + 2*omega*clats
    thyng = dwdlon/clats
    del dwdlon
    thyng = drudr - thyng
    del drudr
    thyng = (2*omega*clats + thyng/rr)*dthetadlat/rr
    del dthetadlat
    thing = thing + thyng
    del thyng
    print("bunch 8")
    
    
    thyng = dvdlon - duclatdlat
    del dvdlon
    del duclatdlat
    thyng = thyng/(rr*clats)
    del rr
    del clats
    slats = slt[np.newaxis,np.newaxis,:,np.newaxis]*np.ones(thing.shape)
    thyng = (2*omega*slats + thyng)*dthetadr
    del dthetadr
    #del slats
    thing = thing + thyng
    del thyng
    thing = thing/rho
    del rho
    #thing = pvort
    
    print(np.nanmin(thing))
    print(np.nanmax(thing))
    
    latvort = np.nanmean(thing,axis=(0,3))
    
    print("potential vorticity")
    
    #print "De-trending...."
    #for k in range(nlevs):
    #    thing[:,k,:,:] = thing[:,k,:,:]/latvort[np.newaxis,k,:,np.newaxis]
    #    print "Level %d"%k
    thing = thing/slats
    del slats    
        
    ncd = nc.Dataset(name,"r+")
    print("Opening variable")
    try:
        pv = ncd.createVariable("pvort",thing.dtype,["time","lev","lat","lon"])
    except:
        pv = ncd.variables["pvort"]
    print("masking....")
    pv.set_auto_mask(False)
    print("Filling....")
    pv[:] = thing[:]#/latvort[np.newaxis,np.newaxis,:,np.newaxis]
    print("Syncing....")
    ncd.sync()
    print("Closing!")
    ncd.close()
                
    print("Modified original file in place; we have reached the end")
    

if __name__=="__main__":
    
    daylen = 1.0/float(get_namelist("planet_namelist","ROTSPD"))
    if daylen<0:
        daylen = 1.0
        
    os.system("[ -e MOST_HC ] && ./burn7.x -n <snapshot.nl>hcout MOST_HC highcadence/MOST_HC.nc")
    os.system("[ -e highcadence/MOST_HC.nc ] && rm MOST_HC")
    rpvorticity("highcadence/MOST_HC.nc",daylen)
    
    
    snapfiles = sorted(glob.glob("*_SNAP.*"))
    if len(snapfiles)>0:
        for snapf in snapfiles:
            try:
                os.system("[ -e "+snapf+" ] && ./burn7.x -n <snapshot.nl>snapout "+snapf+" snapshots/"+snapf+".nc")
                rpvorticity("snapshots/%s.nc"%snapf,daylen)
                os.system("[ -e snapshots/"+snapf+".nc ] && rm "+snapf)
            
            except:
                pass
    outfiles = sorted(glob.glob("*_OUT.????"))
    if len(outfiles)>0:
        for outf in outfiles:
            try:
                os.system("[ -e "+outf+" ] && ./burn7.x -n <example.nl>burnout "+outf+" "+outf+".nc")
                os.system("[ -e "+outf+".nc ] && rm "+outf)
            except:
                pass
    
    year = len(outfiles)
    
    keepon=True
    
    if hasnans():
        os.system("echo 'NAN ENCOUNTERED'>>weathering.pso")
        os.system("rm keepgoing")
        keepon=False
    if energybalanced(year,threshhold=4.0e-4):
        os.system("rm keepgoing")
        keepon=False
    sb,tb = getbalance()
    os.system("echo '%02.6f  %02.6f'>>balance.log"%(sb,tb))
    if keepon:
        bott = gethistory(key="hfns")
        topt = gethistory(key="ntr")
        with open("shistory.pso","a+") as f:
            text='\n'+'\n'.join(bott.astype(str))
            f.write(text)
        with open("toahistory.pso","a+") as f:
            text='\n'+'\n'.join(topt.astype(str))
            f.write(text)
            
        
        
        