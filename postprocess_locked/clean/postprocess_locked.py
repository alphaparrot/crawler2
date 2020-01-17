import numpy as np
import sys
import netCDF4 as nc
import time
import colormatch
  
def postprocess(lats,lons,times,angles,makemap=True,color=True):
  collect(lats,lons,times,angles)
  specs = nc.Dataset("spectra.nc","r")
  makephase(specs,times,angles=angles,makemap=makemap,color=color)

def postphases(lats,lons,times,angles,makemap=True,color=True):  
  specs = nc.Dataset("spectra.nc","r")
  makephase(specs,times,angles=angles,makemap=makemap,color=color)

def readradiance(filename):
  fl = open(filename,"r")
  dd = fl.read().split('\n')[:-1]
  typeid = dd[1]
  if typeid == '"tbf':
    try:
        nrecs = int(dd[2].split()[0])
    except:
        print filename
        print dd[2]
        raise
    dd = dd[3:]
    nzens = int(dd[1].split()[1])
    nphis = int(dd[1].split()[0])
    
    wvls = np.zeros(nrecs)
    insol = np.zeros(nrecs)
    rads = np.zeros((nphis,nzens,nrecs))
    
    inc=1
    if nzens>6:
        zens = dd[3] + dd[4]
        inc=2
    else:
        zens = dd[3]
    zens = zens.split()
    for z in range(0,len(zens)):
        zens[z] = float(zens[z])
    zens = np.array(zens)
    
    phis = dd[2].split()
    for p in range(0,len(phis)):
        phis[p] = float(phis[p])
    phis = np.array(phis)
    
    njump = 3+inc+len(zens)
    
    for n in range(0,nrecs):
      wvls[n] = float(dd[njump*n].split()[0])
      insol[n] = float(dd[njump*n].split()[2])
      for k in range(3+inc,njump):
        for a in range(0,nphis):
            try:
                rads[a,k-3-inc,n] = float(dd[njump*n+k].split()[a])
            except:
                print a,k-3-inc,n,njump,k,len(dd),rads.shape
                print njump*n+k
                print len(dd[njump*n+k].split())
                raise
    
    bundle = {"type":typeid[1:],
              "phis":phis,
              "zens":zens,
              "wvl":wvls,
              "radiance":rads,
              "input":insol}
    return bundle

  elif typeid == '"fzw':
    n=4
    altz = np.array(''.join(dd[7:9]).split()).astype(float)[::-1]
    wvls = []
    while n<len(dd): #first pass, just count records
      line = dd[n].split()
      if len(line)==1:
          wvls.append(float(line[0]))
          n+=18
      else:
          n+=1
    nrecs = len(wvls)
    wvls = np.array(wvls)
    ddir = np.zeros((nrecs,len(altz)))
    ddif = np.zeros((nrecs,len(altz)))
    dtot = np.zeros((nrecs,len(altz)))
    utot = np.zeros((nrecs,len(altz)))
    n=4
    iw=0
    while n<len(dd): #rewind and pass again
        line = dd[n].split()
        if len(line)==1:
            n+=5
            ddir[iw,:] = np.array(''.join(dd[n:n+2]).split()).astype(float)
            n+=3
            ddif[iw,:] = np.array(''.join(dd[n:n+2]).split()).astype(float)
            n+=3
            dtot[iw,:] = np.array(''.join(dd[n:n+2]).split()).astype(float)
            n+=3
            utot[iw,:] = np.array(''.join(dd[n:n+2]).split()).astype(float)
            n+=2
            iw+=1
        else:
            n+=1
            
    ddir = np.transpose(ddir)
    ddif = np.transpose(ddif)
    dtot = np.transpose(dtot)
    utot = np.transpose(utot)
    bundle = {"type":typeid[1:],
              "altz":altz,
              "fluxdr":ddir,
              "fluxdf":ddif,
              "fluxdt":dtot,
              "fluxut":utot,
              "wvl":wvls}
    return bundle
  else:
    return {"type":-1}

def collect(lats,lons,times,angles,modelfile="spectra.nc"):
  
  views = np.array([0.0,90.0,180.0,270.0])
  
  lon0 = views[times[0]]
  nlon0 = np.where(abs(lons-lon0) == np.amin(abs(lons-lon0)))[0][0]
  
  t1 = time.time()
  
  nlats = len(lats)
  nlons = len(lons)
  ntimes = max(len(times),times.max()+1)
  nangles = len(angles)
  
  print nlats,nlons
  
  bundle = readradiance("sbout.%02d_%02d_%1d_%s"%(nlats/2,nlon0,times[0],angles[0]))    #substellar
  if bundle["type"]=="tbf":
    phis = bundle["phis"]
    zens = bundle["zens"]
    wvls = bundle["wvl"]
    fluxes = bundle["radiance"]
    insol = bundle["input"]
    zflux=False
  elif bundle["type"]=="fzw":
    phis = np.array((0.,))
    zens = np.array((0.,))
    wvls = bundle["wvl"]
    altz = bundle["altz"]
    fluxes = bundle["fluxut"]
    insol = bundle["fluxdt"][-1,:]
    zflux=True
  else:
      return bundle["type"]
  nwvs = len(wvls)
  for nln in range(0,3*nlons/4):
      bundle = readradiance("sbout.%02d_%02d_%1d_%s"%(nlats/2,nln,times[0],angles[0]))
      if bundle["type"]=="tbf":
        p = bundle["phis"]
        z = bundle["zens"]
        w = bundle["wvl"]
        f = bundle["radiance"]
        i = bundle["input"]
      elif bundle["type"]=="fzw":
        p = np.array((0.,))
        z = np.array((0.,))
        w = bundle["wvl"]
        altz = bundle["altz"]
        f = bundle["fluxut"]
        i = bundle["fluxdt"][-1,:]
      if len(w)>nwvs:
          nwvs = len(w)
          wvls = w
          fluxes = f
          zens = z
          phis = p
          insol = i
          break
    
  model = nc.Dataset(modelfile, "w", format="NETCDF4")
  if zflux:
      altd = model.createDimension("alt",len(altz))
  lat = model.createDimension("lat",nlats)
  lon = model.createDimension("lon",nlons)
  #phi = model.createDimension("phi",len(phis))
  #zen = model.createDimension("zen",len(zens))
  ttime = model.createDimension("time",ntimes)
  phi = model.createDimension("phi",nangles)
  zen = model.createDimension("zen",nangles)
  wvl = model.createDimension("wavelength",len(wvls))
  
  latitudes   = model.createVariable("lat","f4",("lat",),zlib=True,least_significant_digit=6)
  longitudes  = model.createVariable("lon","f4",("lon",),zlib=True,least_significant_digit=6)
  azimuths    = model.createVariable("azm","f4",("lat","lon","time","phi",),
                                     zlib=True,least_significant_digit=6)
  zeniths     = model.createVariable("zen","f4",("lat","lon","time","zen",),
                                     zlib=True,least_significant_digit=6)
  if zflux:
      altitudes = model.createVariable("alt","f4",("alt",),zlib=True,least_significant_digit=6)
  wavelengths = model.createVariable("wvl","f4",("wavelength",),zlib=True,least_significant_digit=6)
  if not zflux:
    synthspec = model.createVariable("spectra","f4",("lat","lon","time","zen","wavelength",),
                                     zlib=True,least_significant_digit=6)
  else:
    synthspec = model.createVariable("spectra","f4",("alt","lat","lon","time","zen","wavelength",), 
                                     zlib=True,least_significant_digit=6)
  insolation = model.createVariable("input","f4",("lat","lon","time","zen","wavelength",),
                                     zlib=True,least_significant_digit=6)
  
  model.set_auto_mask(False)
  latitudes.set_auto_mask(False)  
  longitudes.set_auto_mask(False)   
  azimuths.set_auto_mask(False)     
  zeniths.set_auto_mask(False)      
  wavelengths.set_auto_mask(False)  
  synthspec.set_auto_mask(False)    
  insolation.set_auto_mask(False)
  if zflux:
      altitudes.set_auto_mask(False)
      altitudes.units = "km"
      altitudes[:] = altz.astype("float32")
  
  wavelengths.units = "microns"
  latitudes.units = "degrees"
  longitudes.units = "degrees"
  azimuths.units = "degrees"
  zeniths.units = "degrees"
  synthspec.units = "W/m^2/um/sr"
  insolation.units = "W/m^2/um"
  
  latitudes[:] = lats.astype("float32")
  longitudes[:] = lons.astype("float32")
  azimuths[nlats/2,0,times[0],0] = phis.astype("float32")
  zeniths[nlats/2,0,times[0],0] = zens.astype("float32")
  wavelengths[:] = wvls.astype("float32")
  
  synthspec[:] = 0.0
  
  if not zflux:
    synthspec[nlats/2,nln,times[0],0,:] = fluxes[0,0,:].astype("float32")
  else:
    synthspec[:,nlats/2,nln,times[0],0,:] = fluxes[:,:]
  
  insolation[:] = 0.0
  insolation[nlats/2,nln,times[0],0,:] = insol[:].astype("float32")
  
  if not zflux:
    fullspec = np.zeros((nlats,nlons,ntimes,len(angles),len(wvls)))
    fullspec[nlats/2,times[0],0,0,:] = fluxes[0,0,:]
  else:
    fullspec = np.zeros((len(altz),nlats,nlons,ntimes,len(angles),len(wvls)))
    fullspec[:,nlats/2,times[0],0,0,:] = fluxes[:,:]  
      
  print synthspec.mask
  
  for kt in times:
    for ka in range(0,nangles):
      for jlat in range(0,nlats):
        for jlon in range(0,nlons):
          print "Collecting lat %02d and lon %02d from view %s of time %1d"%(jlat,jlon,angles[ka],kt)
          try:
              phis,zens,wvls,fluxes,insol 
              bundle = readradiance("sbout.%02d_%02d_%1d_%s"%(jlat,jlon,kt,angles[ka]))
              if bundle["type"]=="tbf":
                phis = bundle["phis"]
                zens = bundle["zens"]
                wvls = bundle["wvl"]
                fluxes = bundle["radiance"]
                insol = bundle["input"]
              elif bundle["type"]=="fzw":
                phis = np.array((0.,))
                zens = np.array((0.,))
                wvls = bundle["wvl"]
                fluxes = bundle["fluxut"]
                insol = bundle["fluxdt"][-1,:]
          except:
              fluxes[:] = 0.0
              insol[:] = 0.0
              print "Spectrum %02d_%02d_%1d_%s is missing!"%(jlat,jlon,kt,angles[ka])
          n_init=0
          n_end =len(wavelengths[:])
          if len(wvls)<len(wavelengths[:]):
            n_init = n_end-len(wvls)
          if not zflux:
            print synthspec[jlat,jlon,:,:,n_init:n_end].shape,fluxes[:,:,:].shape
            synthspec[jlat,jlon,kt,ka,n_init:n_end] = fluxes[0,0,:].astype("float32")
            insolation[jlat,jlon,kt,ka,n_init:n_end] = insol[:].astype("float32")
            print n_init,n_end,np.mean(fluxes),np.mean(synthspec[jlat,jlon,:,:,:])
            fullspec[jlat,jlon,kt,ka,n_init:n_end] = fluxes[0,0,:]
            azimuths[jlat,jlon,kt,ka] = phis.astype("float32")
            zeniths[jlat,jlon,kt,ka] = zens.astype("float32")
          else:
            print synthspec[-1,jlat,jlon,:,:,n_init:n_end].shape,fluxes[:,:].shape
            #print insolation[jlat,jlon,kt,ka,n_init:n_end].shape,insol[:].shape
            synthspec[:,jlat,jlon,kt,ka,n_init:n_end] = fluxes[:,:].astype("float32")
            insolation[jlat,jlon,kt,ka,n_init:n_end] = insol[:].astype("float32")
            print n_init,n_end,np.mean(fluxes),np.mean(synthspec[-1,jlat,jlon,:,:,:])
            fullspec[:,jlat,jlon,kt,ka,n_init:n_end] = fluxes[:,:]
            azimuths[jlat,jlon,kt,ka] = phis.astype("float32")
            zeniths[jlat,jlon,kt,ka] = zens.astype("float32")
            
        
  model.sync()
  model.close() 

  t2 = time.time()
  print "Time elapsed: ",str(t2-t1),"seconds"

  #return fullspec

def lognorm(x):
    v = np.log10(np.maximum(x,1.0e-15))
    vmin = np.amin(v)
    v-=vmin
    vmax = np.amax(v)
    v/=vmax
    return v

def diskintegrate(spectra,ntime,ndir,planet_radius=1.0,makemap=True,colorize=True):
  lts = spectra.variables['lat'][:]
  lns = spectra.variables['lon'][:]
  
  nlats = len(lts)
  nlons = len(lns)
  
  faces = ['Z','E','W','N','S']
  decs = np.array([0.,0.,0.,90.,-90.])
  views = np.array([0,90,-90,0,0])+180.0
  view = views[ndir]
 
  darea = np.zeros((nlats,nlons))
  nlats = len(lts)
  lts1 = np.zeros(nlats)
  lts2 = np.zeros(nlats)
  lts1[0] = 0.5*np.pi
  lts1[nlats-1] = 0.5*(lts[nlats-2]+lts[nlats-1])*np.pi/180.0
  lts2[0] = 0.5*(lts[0]+lts[1])*np.pi/180.0
  lts2[nlats-1] = -0.5*np.pi
  for jlat in range(1,nlats-1):
      lts1[jlat] = 0.5*(lts[jlat-1]+lts[jlat])*np.pi/180.0
      lts2[jlat] = 0.5*(lts[jlat+1]+lts[jlat])*np.pi/180.0
  for jlat in range(0,nlats):
      darea[jlat,:] = 0.5/nlons*(np.sin(lts1[jlat])-np.sin(lts2[jlat]))
  
  #plarad = planet_radius * 6371000.0 #We are ignoring the fact that the light we're seeing is 20-100km up
  #darea *= 4*np.pi*plarad**2 #Turn areas into physical units (square meters)
  globalarea = np.sum(darea)
  
  #darea /= globalarea #Actually let's keep this in W/m^2
  
  if "alt" in spectra.variables:
      zflux=True
      nz = len(spectra.variables['alt'][:])
  else:
      zflux=False
  
  if zflux:
    diskspec = np.zeros((nz,len(spectra.variables['wvl'][:])))
  else:
    diskspec = np.zeros(len(spectra.variables['wvl'][:]))
  diskinspec = np.zeros(len(spectra.variables['wvl'][:]))
  if makemap:
    if zflux:
      phasemap = np.zeros((nz,nlats,nlons,len(spectra.variables['wvl'])))
    else:
      phasemap = np.zeros((nlats,nlons,len(spectra.variables['wvl'])))
    phaseinmap = np.zeros((nlats,nlons,len(spectra.variables['wvl'])))
  else:
    phasemap = None
    phaseinmap = None
  if colorize:
    if zflux:
      colors = np.zeros((nz,nlats,nlons,3))
      intens = np.zeros((nz,nlats,nlons,3))
    else:
      colors = np.zeros((nlats,nlons,3))
      intens = np.zeros((nlats,nlons,3))
  else:
    colors = None
  
  va = 180.0-view-lns #if view=0, lon=0 should be facing away, at 180 degrees.
                            #if view=120, lon=0 is 60 degrees from zenith.
                            #if view=40, lon=0 is 50 degrees past the terminator (140 degrees)
                            #if view=200, lon=0 is 20 degrees from zenith.
  va[va<0] = (va[va<0])%(-360)
  va[va>0] = (va[va>0])%(360)
  va[va<-180] += 360
  va[va>180]  -= 360  
  dec = decs[ndir] * np.pi/180.0
  
  ilat = lts*np.pi/180.0
  
  vas,ilats = np.meshgrid(va*np.pi/180.,ilat)
  
  coszenview = np.sin(dec)*np.sin(ilats)+np.cos(ilats)*np.cos(dec)*np.cos(vas)
      #This is the cosine of the zenith angle, right? Should we turn this back to degrees?
      #We will need to turn this into projected area, but we need to select the scattering
      #angle based on the viewing angle, which is the zenith angle for each point at that view.
      #We select the azimuth based on whether the point is on the left or right side of the planet.
      #If the point is on the left side, then the azimuth is 90 degrees; if on the right, it's 270.
      #Basically the azimuth is whether we're to the East or the West of that grid point.
  
  zenview = np.arccos(coszenview)
  sinzenview = np.sin(zenview)
  
  cosazs = (np.sin(dec)-coszenview*np.sin(ilats)) / (sinzenview*np.cos(ilats))
  cosazs[np.isnan(cosazs)] = 0.0
  uazs = np.arccos(np.minimum(np.maximum(cosazs,-1.0),1.0)) * 180.0/np.pi
  
  avguaz = np.mean(uazs)  
  
  projectionweight=np.maximum(np.cos(ilats)*np.cos(vas),0.0)
  
  norm = 1.0/np.sum(projectionweight*darea)
  
  for jlat in range(0,nlats):  
    
    for jlon in range(0,nlons):
  
      #if va[jlon] > 0:
        #avuaz = 360.0-avguaz
      #else:
        #avuaz = avguaz
        
      #zens = spectra.variables['zen'][jlat,jlon,:] #16 length
      #phis = spectra.variables['azm'][jlat,jlon,:] #2 length
      
      #iz = closest(zens,zenview[jlat,jlon]*180./np.pi)
      #ip = closest(phis,avuaz)
      
      if zflux:
        spec = spectra.variables['spectra'][:,jlat,jlon,ntime,ndir,:]
      else:
        spec = spectra.variables['spectra'][jlat,jlon,ntime,ndir,:]
          
      inspec = spectra.variables['input'][jlat,jlon,ntime,ndir,:]
      if makemap:
        if zflux:
          phasemap[:,jlat,jlon,:] = spec
        else:
          phasemap[jlat,jlon,:] = spec
        phaseinmap[jlat,jlon,:] = inspec
      if colorize:
        if not zflux:
          intens[jlat,jlon,:] = colormatch.makexyz(spectra.variables['wvl'][:]*1000.0,spec)
        else:
          for jz in range(0,nz):
            intens[jz,jlat,jlon,:] = colormatch.makexyz(spectra.variables['wvl'][:]*1000.0,spec[jz])
      
      spec *= projectionweight[jlat,jlon]
      
      spec *= darea[jlat,jlon]*norm
      
      inspec *= projectionweight[jlat,jlon]
      inspec *= darea[jlat,jlon]*norm
      
      
      diskspec += spec #return intensity in W/m^2/um
      diskinspec += inspec
  
  if colorize:
    if not zflux:
      norms = lognorm(intens[:,:,2]/np.amax(intens[:,:,2]))
      for j in range(0,nlats):
          for k in range(0,nlons):
              colors[j,k,:] = colormatch.xyz2rgb(intens[j,k,0],intens[j,k,1],
                                                 intens[j,k,2]*norms[j,k])
        
      colors /= np.amax(colors)
    else:
      for z in range(0,nz):
        norms = lognorm(intens[z,:,:,2]/np.amax(intens[z,:,:,2]))
        for j in range(0,nlats):
            for k in range(0,nlons):
                colors[z,j,k,:] = colormatch.xyz2rgb(intens[z,j,k,0],intens[z,j,k,1],
                                                   intens[z,j,k,2]*norms[j,k])
          
        colors[z,:,:] /= np.amax(colors[z,:,:])
  
  return diskspec,phasemap,colors,diskinspec,phaseinmap



def makephase(spectra,times,angles=['Z',],outfile="phases.nc",planet_radius=1,color=True,makemap=True):
  nlats = len(spectra.variables['lat'][:])
  nlons = len(spectra.variables['lon'][:])
  if "alt" in spectra.variables:
      nz = len(spectra.variables["alt"][:])
      altz = spectra.variables["alt"][:]
      zflux=True
  else:
      zflux=False
  ntimes = max(len(times),times.max()+1)
  nangles = len(angles)
  phasecurve = nc.Dataset(outfile,"w", format="NETCDF4")
  time = phasecurve.createDimension("time",ntimes)
  angle = phasecurve.createDimension("angle",nangles)
  if zflux:
    altsz = phasecurve.createDimension("alt",nz)
  wvl = phasecurve.createDimension("wavelength",len(spectra.variables['wvl'][:]))
  if makemap or color:
    lt = phasecurve.createDimension("latitude",nlats)
    ln = phasecurve.createDimension("longitude",nlons)
  cl = phasecurve.createDimension("color",3)
  
  wvls = spectra.variables['wvl'][:]
  if makemap or color:
    latitudes   = phasecurve.createVariable("lat","f4",("latitude",),zlib=True,least_significant_digit=6)
    longitudes  = phasecurve.createVariable("lon","f4",("longitude",),zlib=True,least_significant_digit=6)
  obstimes = phasecurve.createVariable("time","f4",("time",),zlib=True,least_significant_digit=6)
  phaseangles = phasecurve.createVariable("phase",str,("angle",))
  wavelengths = phasecurve.createVariable("wvl","f4",("wavelength",),zlib=True,least_significant_digit=6)
  if not zflux:
    synthspec   = phasecurve.createVariable("spectra","f4",("time","angle","wavelength",),
                                            zlib=True,least_significant_digit=6)
  else:
    altitude    = phasecurve.createVariable("alt","f4",("alt",),zlib=True,least_significant_digit=6)
    synthspec   = phasecurve.createVariable("spectra","f4",("time","angle","alt","wavelength",),
                                            zlib=True,least_significant_digit=6)
  inspec   = phasecurve.createVariable("input","f4",("time","angle","wavelength",),
                                       zlib=True,least_significant_digit=6)
  if makemap:
    if zflux:
      phasemaps = phasecurve.createVariable("map","f4",("time","angle","alt","latitude",
                                                        "longitude","wavelength",),
                                            zlib=True,least_significant_digit=6)
    else:
      phasemaps = phasecurve.createVariable("map","f4",("time","angle","latitude",
                                                        "longitude","wavelength",),
                                            zlib=True,least_significant_digit=6)
    phaseinmaps = phasecurve.createVariable("inmap","f4",("time","angle","latitude",
                                                          "longitude","wavelength",),
                                            zlib=True,least_significant_digit=6)
  if color:
    if zflux:
      phasecols = phasecurve.createVariable("colors","f4",("time","angle","alt","latitude",
                                                           "longitude","color",),
                                            zlib=True,least_significant_digit=6)
    else:
      phasecols = phasecurve.createVariable("colors","f4",("time","angle","latitude",
                                                           "longitude","color",),
                                            zlib=True,least_significant_digit=6)
  if makemap or color:
    latitudes.set_auto_mask(False)  
    longitudes.set_auto_mask(False) 
  phasecurve.set_auto_mask(False)
  wavelengths.set_auto_mask(False)
  phaseangles.set_auto_mask(False)
  obstimes.set_auto_mask(False)
  synthspec.set_auto_mask(False)
  inspec.set_auto_mask(False)
  if makemap:
    phasemaps.set_auto_mask(False)
    phaseinmaps.set_auto_mask(False)
  if color:
    phasecols.set_auto_mask(False)
  if zflux:
    altitude.set_auto_mask(False)
    altitude[:] = altz.astype("float32")
    altitude.units = "km"
  if makemap or color:
    latitudes[:] = spectra.variables['lat'][:].astype("float32")
    longitudes[:] =spectra.variables['lon'][:].astype("float32")
  synthspec[:] = 0.
  if makemap:
    phasemaps[:] = 0.
  if color:
    phasecols[:] = 0.
  
  if makemap or color:
    latitudes.units = "degrees"
    longitudes.units = "degrees"
  wavelengths.units = "microns"
  phaseangles.units = "direction"
  obstimes.units = "int"
  synthspec.units = "W/m^2/um"
  inspec.units = "W/m^2/um"
  if makemap:
    phasemaps.units = "W/m^2/um"
    phaseinmaps.units = "W/m^2/um"
  if color:
    phasecols.units = "1"
  
  obstimes[times] = np.array(times).astype("float32")
  if len(angles)>1:
    for nal in range(0,len(angles)):
      phaseangles[nal] = angles[nal]
  else:
    phaseangles[0] = angles[0]
  wavelengths[:] = wvls.astype("float32")
  
  for t in times:
    for d in range(0,nangles):
      di,pm,colors,ii,pim = diskintegrate(spectra,t,d,planet_radius=planet_radius,makemap=makemap,colorize=color)
      if zflux:
        synthspec[t,d,:,:] = di[:]
      else:
        synthspec[t,d,:] = di[:]
      inspec[t,d,:] = ii[:]
      if makemap:
          if zflux:
            phasemaps[t,d,:,:,:,:] = pm[:]
          else:
            phasemaps[t,d,:,:,:] = pm[:]
          phaseinmaps[t,d,:,:,:] = pim[:]
      if color:
          if zflux:
            phasecols[t,d,:,:,:,:] = colors[:]
          else:
            phasecols[t,d,:,:,:] = colors[:]
  phasecurve.sync()  
  phasecurve.close()
  

def closest(data,val):
  try:
    return np.where(abs(data-val)==np.amin(abs(data-val)))[0][0]
  except:
    return np.where(abs(data-val)==np.amin(abs(data-val)))[0]


if __name__=="__main__":
    lat = np.load("latitudes.npy")
    lon = np.load("longitudes.npy")
    
    times = sys.argv[1]
    angles = sys.argv[2]
    
    times = times.split('^')
    if times[-1]=='':
        times = times[:-1]
    times = np.array(times).astype(int)
    angles = angles.split('^')
    if angles[-1]=='':
        angles = angles[:-1]
    
    color=False
    makemap=False
    if "color" in sys.argv[:]:
        color=True
    if "map" in sys.argv[:]:
        makemap=True
    if "phases" in sys.argv[:]:
        postphases(lat,lon,times,angles,color=color,makemap=makemap)
    else:
        postprocess(lat,lon,times,angles,color=color,makemap=makemap)
