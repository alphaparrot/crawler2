import os
import numpy as np
import time
from batch_system import SUB, BATCHSCRIPT
from identity import USER, SCRATCH
from crawldefs import Job
import exoplasim as exo

# Options:
#   noutput
#   pCO2     (ubars)
#   pressure
#   flux
#   gravity
#   radius
#   year
#   restart
#   ncarbon
#   volcanCO2
#   wmax
#   nglacier
#   glacelim
#   sidday   (hours)
#   solday   (hours)
#   sidyear  (days)
#   solyear  (days)
#   rotspd
#   obliq
#   eccen
#   vernlon
#   radius
#   snowmax
#   soilh2o
#   naqua
#   script
#   extra
#   nfixorb
#   nwpd     (writes per day)
#   months   (runtime)
#   days     (runtime)
#   steps    (runtime)
#   timestep (minutes/timestep)
#   runyears
#   fixedlon (longitude of fixed substellar point)
#   source
#   notify



def prep(job):

  sig=job.name
  jid=str(job.home)
  pid=job.pid
  args=job.args
  fields=job.fields
  
  workdir = "exoplasim/job"+jid
  
  if "source" in job.parameters:
    source = job.parameters["source"]
  else:
    source = "clean"

  if "westeros" in job.parameters:
    source = "westeros"
    
  if "source" in job.parameters:
    source = job.parameters["source"]
  
  print "Setting stuff for job "+sig+" in plasim/job"+jid+" which is task number "+pid
  print "Arguments are:",fields[2:]
  
  notify = 'ae'
  scriptfile = "super_relax.py"
  
  #PlaSim
  #os.system("rm plasim/job"+jid+"/*.sh")
  #os.system("mkdir plasim/errorlogs/")
  #os.system("cp exoplasim/job"+jid+"/*.e* plasim/errorlogs/")
  #os.system("rm -rf plasim/job"+jid+"/*")
  #os.system("cp exoplasim/"+source+"/* plasim/job"+jid+"/")
  #os.system("rm plasim/job"+jid+"/plasim_restart")
  #os.system("cp exoplasim/"+scriptfile+" plasim/job"+jid+"/")
  
  nlevs = 10
  if "nlevs" in job.parameters:
      nlevs = int(job.parameters["nlevs"])
      
  mars = False    
  if "mars" in job.parameters:
      val = job.parameters["mars"]
      if val=="True" or val=="1":
          mars=True
          
  if "year_init" in job.parameters:
      yearini = int(job.parameters["year_init"])
  
  resolution = "21"
  nlats = 32
  
  if "resolution" in job.fields[2:]:
    val = job.parameters["resolution"]
    if val=="T21" or val=="t21" or val=="21" or val=="32":
        resolutions="21"
        nlats=32
    if val=="T42" or val=="t42" or val=="42" or val=="64":
        resolutions="42"
        nlats=64
    if val=="T63" or val=="t63" or val=="63" or val=="63":
        resolutions="63"
        nlats=96
    if val=="T85" or val=="t85" or val=="85" or val=="128":
        resolutions="85"
        nlats=128
    if val=="T127" or val=="t127" or val=="127" or val=="192":
        resolutions="127"
        nlats=192
    if val=="T170" or val=="t170" or val=="170" or val=="256":
        resolutions="170"
        nlats=256
  
  if "lockedyear" in job.parameters or "locked" in job.parameters:
    model = exo.TLmodel(workdir=workdir,ncpus=job.ncores,modelname=job.name,
                        layers=nlevs,source=source,mars=mars,inityear=yearini,
                        resolution="T%s"%resolutions)
  else:
    model = exo.Model(workdir=workdir,ncpus=job.ncores,modelname=job.name,layers=nlevs,
                    source=source,mars=mars,inityear=yearini,resolution="T%s"%resolutions)
  
  model.configure()
  
  os.system("cp exoplasim/synthoutput.py %s/"%model.workdir)
  os.system("cp crawldefs.py %s/"%model.workdir)
  os.system("cp identity.py %s/"%model.workdir)
  
  if "cleanup" in job.parameters:
    cleanup = job.parameters["cleanup"]
  else:
    cleanup = 'release-plasim.py'
  os.system("cp exoplasim/"+cleanup+" plasim/job"+jid+"/")
  
  tag = 'all'
  
  keeprs = False
  monitor = False
  weathrestart = False
  nsn = False
  nwesteros = False
  setgas=False
  setgasx=False
  maketransit = False
  makesbdart = False
  hcout = False
  highcadence = False
  stormclim = False
  prescgascon = False
  setglacier = False
  alloutput = False
  
  sbdvar = "earth"
  plarad = 1.0
  grav = 9.80665
  starrad = 1.0
  
  threshhold = 4.0e-4
  
  flux = 1367.0
  startemp = -1
  starfile = False
  
  #SI units:
  hh = 6.62607004e-34
  cc = 2.99792458e8
  kb = 1.38064852e-23
  rsun = 6.95510e8 #meters
  au = 1.496e11
  adjfac = 3.1011857558763545
  
  for name in job.fields[2:]:
    val = job.parameters[name]
    found=False    
        
    if name=="noutput":
      model.modify(noutput=bool(int(val)))
      found=True
      
    if name=="flux":
      model.modify(flux=float(val))
      flux = float(val)
      found=True
      
    if name=="startemp":
      model.modify(startemp=float(val))
      startemp = float(val)
      found=True
      
    if name=="starspec":
      model.modify(starspec=val)
      starfile = val[:-4]+"_hr.dat"
      found=True
      
    if name=="starrad": #radius of the star in solar radii--only used for transit calc
      found=True
      starrad = float(val)
    
    if name=="transit":
      maketransit=True
      found=True
      
    if name=="sbdart":
      makesbdart=True
      hcout=False
      found=True
      ntimes = val.split('|')[0]
      lviews = val.split('|')[1]
      
    if name=="sbdart_hc":
      makesbdart=True
      found=True
      hcout=True
      ntimes = "{0,}"
      lviews = val
      
      
    if name=='pH2u': #in ubars
      found=True
      model.modify(pH2=float(val)*1.0e-6
            
    if name=='pH2b': #in bars
      found=True
      model.modify(pH2=float(val)
      
    if name=='pHeu': #in ubars
      found=True
      model.modify(pHe=float(val)*1.0e-6
            
    if name=='pHeb': #in bars
      found=True
      model.modify(pHe=float(val)
          
    if name=='pN2u': #in ubars
      found=True
      model.modify(pN2=float(val)*1.0e-6
            
    if name=='pN2b': #in bars
      found=True
      model.modify(pN2=float(val)
      
    if name=='pO2u': #in ubars
      found=True
      model.modify(pO2=float(val)*1.0e-6
            
    if name=='pO2b': #in bars
      found=True
      model.modify(pO2=float(val)
    
    if name=='pCO2u': #in ubars
      found=True
      model.modify(pCO2=float(val)*1.0e-6
            
    if name=='pCO2b': #in bars
      found=True
      model.modify(pCO2=float(val)
      
    if name=='pAru': #in ubars
      found=True
      model.modify(pAr=float(val)*1.0e-6
            
    if name=='pArb': #in bars
      found=True
      model.modify(pAr=float(val)
      
    if name=='pNeu': #in ubars
      found=True
      model.modify(pNe=float(val)*1.0e-6
            
    if name=='pNeb': #in bars
      found=True
      model.modify(pNe=float(val)
      
    if name=='pKru': #in ubars
      found=True
      model.modify(pKr=float(val)*1.0e-6
            
    if name=='pKrb': #in bars
      found=True
      model.modify(pKr=float(val)
      
    if name=='pH2Ou': # NOTE this will ONLY change the mmw--no affect on moist processes
      found=True
      model.modify(pH2O=float(val)*1.0e-6
      
    if name=='pH2Ob': # NOTE this will ONLY change the mmw--no affect on moist processes
      found=True
      model.modify(pH2O=float(val)
    
    #if name=='xH2': # mass fraction
      #found=True
      #setgasx=True
      #gasesx['H2'] = eval(val)
            
      
    #if name=='xHe': # mass fraction
      #found=True
      #setgasx=True
      #gasesx['He'] = eval(val)
            
          
    #if name=='xN2': # mass fraction
      #found=True
      #setgasx=True
      #gasesx['N2'] = eval(val)
            
      
    #if name=='xO2': # mass fraction
      #found=True
      #setgasx=True
      #gasesx['O2'] = eval(val)
    
    #if name=='xCO2': # mass fraction
      #found=True
      #setgasx=True
      #gasesx['CO2'] = eval(val)
            
      
    #if name=='xAr': # mass fraction
      #found=True
      #setgasx=True
      #gasesx['Ar'] = eval(val)
            
      
    #if name=='xNe': # mass fraction
      #found=True
      #setgasx=True
      #gasesx['Ne'] = eval(val)
            
      
    #if name=='xKr': # mass fraction
      #found=True
      #setgasx=True
      #gasesx['Kr'] = eval(val)
      
    #if name=='xH2O': # NOTE this will ONLY change the mmw--no affect on moist processes
      #found=True
      #setgasx=True
      #gasesx['H2O'] = eval(val)
      
    #if name=="pCO2":
      #found=True
      #gotpress=False
      #if "pressure" in fields:
        #p0 = float(job.parameters["pressure"])
        #gotpress=True
      #if not gotpress:
        #p0 += float(val)
        #p0 *= 1.0e-6
        #edit_namelist(jid,"plasim_namelist","PSURF",str(p0*1.0e5))
    
      #pCO2 = float(val)/(p0*1.0e6)*1.0e6 #ppmv
      #gases['pCO2'] = float(val)*1.0e-6
      #gasesx['CO2'] = pCO2*1e-6     #volumefraction
      #edit_namelist(jid,"radmod_namelist","CO2",str(pCO2))      
      
    if name=="pressure": #in bars
      found=True
      model.modify(pressure=float(val))
       
    if name=="alloutput":
      found=True
      alloutput=True
      
    if name=="keeprestart":
      found=True
      if int(job.parameters["keeprestart"])==1:
          keeprs = True
        
    if name=="monitor":
      found=True
      if int(job.parameters["monitor"])==1:
          monitor = True
          
    if name=="pressurebroaden":
      found=True
      model.modify(pressurebroaden=bool(int(val)))
      
    if name=="year": #In days
      found=True
      model.modify(year=float(val))
      
    if name=="vtype":
      found=True
      model.modify(vtype=int(val))
     
    if name=="rotationperiod":
      found=True
      model.modify(rotationperiod=float(val))
      
    if name=="lockedyear": #Year length in days for a tidally-locked planet_namelist, and lon0.
      found=True
      sbdvar="locked"
      val0 = val.split('/')[0]
      val1 = val.split('/')[1]
      model.modify(rotationperiod=float(val0))
      model.modify(substellarlon=float(val1))
      
    if name=="sbdart_type": #Whether to use SBDART in tidally-locked or Earth configuration
      found=True
      sbdvar = val
      
    if name=="restart":
      found=True
      model.modify(restartfile=val)
        
    if name=="hweathering":
      found=True
      os.system("cp hopper/"+val+" plasim/job"+jid+"/")
      wfile = val
      weathrestart = True
      
    if name=="balance":
      found=True
      model.modify(threshold=float(val))
     
    if name=="gravity":
      found=True
      grav = float(val)
      model.modify(gravity=float(val))
       
    if name=="radius":
      found=True
      plarad=float(val)
      model.modify(radius=float(val))
     
    if name=="eccen":
      found=True
      model.modify(eccentricity=float(val))
     
    if name=="obliq":
      found=True
      model.modify(obliquity=float(val))
       
    if name=="vernlon":
      found=True
      model.modify(lonvernaleq=float(val))
           
    if name=="nfixorb":
      found=True
      model.modify(fixedorbit=bool(int(val)))
      
    if name=="oroscale":
      found=True
      model.modify(orography=float(val))
       
    if name=="nglacier":
      found=True
      setglacier=True
      glaciers["toggle"] = bool(int(val))
      
    if name=="glacelim":
      found=True
      setglacier=True
      glaciers["mindepth"] = float(val)
      
    if name=="icesheets":
      found=True
      setglacier=True
      glaciers["initialh"] = float(val)
    
    if name=="lowerANTGRN":
      found=True
      if int(val)==1:
        os.system("cp exoplasim/glacsra/*.sra %s/"%model.workdir)
    
    if name=="nseaice": #Toggle whether sea ice and snow are allowed (1=yes,0=no, 1 is default)
      found=True
      model.modify(seaice=bool(int(val))
    
    if name=="ncarbon":
      found=True
      model.modify(co2weathering=bool(int(val)))
      
    if name=="co2evolve":
      found=True
      model.modify(evolveco2=bool(int(val)))
          
    if name=="filterpower": ##will work with physfilter and 'exp'
      found=True
      model.modify(filterpower = int(val))
        
    if name=="filterkappa": ##will work with physfilter and 'exp'
      found=True
      model.modify(filterkappa = float(val))
        
    if name=="filterLHN0": ##will work with physfilter and 'lh'
      found=True
      model.modify(filterLHN0 = float(val))
        
    if name=="physfilter":
      found=True
      model.modify(physicsfilter = val)
        
    if name=="frictionmod":
      found=True
      model.modify(otherargs={"FRCMOD@plasim_namelist":val})
      
    if name=="qdiff": #spec humiditiy diffusion timescale in days (default=0.1)
      found=True
      model.modify(qdiffusion=float(val))
        
    if name=="tdiff": #temperature diffusion timescale in days (default=5.6)
      found=True
      model.modify(tdiffusion=float(val))
        
    if name=="zdiff": #vorticity diffusion timescale in days (default=1.1)
      found=True
      model.modify(zdiffusion=float(val))
        
    if name=="ddiff": #divergence diffusion timescale in days (default=0.2)
      found=True
      model.modify(ddiffusion=float(val))
      
    if name=="diffpower": #power for wavelength cutoff diffusion (default=2 for T21, 4 for T42)
      found=True
      model.modify(diffusionpower=int(val))
      
    if name=="nhdiff": #Critical wavenumber for diffusion
      found=True
      model.modify(diffusionwaven=int(val))
        
    if name=="nsupply" or name=="wmax":
      found=True
      model.modify(erosionsupplylimit=float(val)) 
      
    if name=="volcanCO2":
      found=True
      model.modify(outgassing=float(val))
    
    if name=="snowicealbedo":
      found=True
      model.modify(snowicealbedo=float(val))
      
    if name=="nbandalbedo":
      found=True
      model.modify(twobandalbedo=bool(int(val)))
      
    if name=="snowmax":
      found=True
      model.modify(maxsnow=float(val))
      
    if name=="soilalbedo":
      found=True
      model.modify(soilalbedo=float(val))
      
    if name=="oceanalbedo":
      found=True
      model.modify(oceanalbedo=float(val))
      
    if name=="oceanalbzad": #zenith-angle dependence for the ocean direct beam reflectivity
      found=True
      model.modify(oceanzenith=val)
          
    if name=="wetsoil":
      found=True
      model.modify(wetsoil=bool(int(val)))
      
    if name=="soilh2o":
      found=True
      model.modify(soilwatercap=float(val))
      
    if name=="naqua":
      found=True
      model.modify(aquaplanet=bool(int(val)))
      
    if name=="ndesert":
      found=True
      model.modify(desertplanet=bool(int(val)))
      
    if name=="soilwetness":
      found=True
      model.modify(soilsaturation=float(val))
      
    if name=="drycore":
      found=True
      model.modify(drycore=bool(int(val)))
        
    if name=="ozone":
      found=True
      model.modify(ozone=bool(int(val)))
      
    if name=="soilheat": # 0.001 J/m^3/K = 1 J/kg/K. Ocean = 4180 J/kg/K = 4.18e6 J/m^3/K
      found=True
      model.modify(cpsoil=float(val))
      
    if name=="soildepth":
      found=True
      model.modify(soildepth=float(val))
      
    if name=="mixedlayer":
      found=True
      model.modify(mldepth=float(val))
      
    if name=="nwpd":
      found=True
      model.modify(writefrequency=int(val))
      
    if name=="months":
      found=True
      model.modify(otherargs={"N_RUN_MONTHS@plasim_namelist":val})
      
    if name=="rhum_c":
      found=True
      model.modify(otherargs={"RCRITMOD@rainmod_namelist":val})
      
    if name=="rhum_m":
      found=True
      model.modify(otherargs={"RCRITSLOPE@rainmod_namelist":val})
      #if negative, cloud formation will be suppressed at high altitudes
      #if positive, cloud formation will be encouraged at high altitudes
      
    if name=="rcrit":
      found=True
      model.modify(ortherargs={"RCRIT@rainmod_namelist":val})
      
    if name=="days":
      found=True
      model.modify(otherargs={"N_RUN_DAYS@plasim_namelist":val})
      
    if name=="steps":
      found=True
      model.modify(otherargs={"N_RUN_STEPS@plasim_namelist":val})
      
    if name=="o3":
      found=True
      model.modify(ozone=bool(int(val)))
      
    if name=="ptop":
      found=True
      model.modify(modeltop=float(val))
      
    if name=="vlin_top":
      found=True
      model.modify(vtype=4,modeltop=float(val))
      
    if name=="vlog_top":
      found=True
      model.modify(vtype=3,modeltop=float(val))
      
    if name=="stratosphere":
      found=True
      ptop2 = val.split('|')
      ptop2 = np.array(ptop2).astype(float)
      ptop = min(ptop2)
      ptop2 = max(ptop2)
      model.modify(stratosphere=True,tropopause=float(ptop),modeltop=float(ptop2))
      
    if name=="timestep":
      found=True
      model.modify(timestep=float(val))
      
    if name=="rayleigh":
      found=True
      model.modify(otherargs={"NEWRSC@radmod_namelist":val})
      
    if name=="script":
      found=True
      scriptfile=val
      os.system("cp exoplasim/"+val+" %s/"%model.workdir)
      
    if name=="postprocessor":
      found=True
      os.system("cp exoplasim/"+source+"/"+val+" %s/example.nl"%model.workdir)
      
    if name=="notify":
      found=True
      notify = val

    if name=="highcadence":
      found=True
      highcadence=True
      vals = val.split('|')
      hcstart = vals[0]
      hcend   = vals[1]
      hcint   = vals[2]
      hc = {"toggle":1,"start":hcstart,"end":hcend,"interval":hcint}
      model.modify(highcadence=hc)
      
    if name=="columnmode":
      found=True
      if val=="-":
          val=None
      model.modify(columnmode=val)
      
    if name=="snapshots":
      found=True
      model.modify(snapshots=int(val))
      nsn = True
      
    if name=="extra":
      found=True
      model.modify(resources=[val,])
      
    if name=="fixedlon":
      found=True
      model.modify(sychronous=True,substellarlon=float(val))
      
    if name=="source":
      found=True #We already took care of it
      
    if name=='cleanup':
      found=True #We already took care of it
      
    #if name=="westeros":
      #found=True
      #edit_namelist(jid,"plasim_namelist","NWESTEROS",val)
      #edit_namelist(jid,"planet_namelist","NFIXORB","1") 
      #if "runyears" in job.parameters:
        #nrunyears = int(job.parameters["runyears"])
      #else:
        #nrunyears = 100
      #nwesteros=True
      
    #if name=="binarystar":
      #found=True
      #edit_namelist(jid,"radmod_namelist","NDOUBLESTAR",val)
      
    #if name=="flux2":
      #found=True
      #edit_namelist(jid,"planet_namelist","GSOL1",val)
      
    #if name=="initialmzv":
      #found=True
      #parts = val.split('/')
      #manom0 = parts[0]
      #x0 = parts[1]
      #v0 = parts[2]
      #edit_namelist(jid,"planet_namelist","MEANANOM0",manom0)
      #edit_namelist(jid,"planet_namelist","ZWZZ",x0)
      #edit_namelist(jid,"planet_namelist","ZWVV",v0)
      
    #if name=="semimajor":
      #found=True
      #edit_namelist(jid,"planet_namelist","SEMIMAJOR",val)
      
    #if name=="maxflux":
      #found=True
      #edit_namelist(jid,"radmod_namelist","SOLMAX",val)
      
    if name=="gascon":
      found=True
      model.modify(gascon=float(val))
      
    if name=="storms":
      found=True
      stormclim=True
      if val=="True" or val=="true" or val=="1":
          model.modify(stormclim=True)
      else:
          model.modify(stormclim=False)
    
    if name=="nstorms":
      found=True
      model.modify(nstorms=int(val))
      
    if name=="stormtrigger":
      found=True
      vals = val.split('|')
      stormdict={"toggle":1}
      for v in vals:
          keys = v.split("=")
          stormdict[keys[0]]=float(keys[1])
      model.modify(stormcapture=stormdict)
    
    if name=="topomap":
      found=True
      
    if not found: #Catchall for options we didn't include
      found=True
      args = name.split('@')
      if len(args)>1:
        model.modify(otherargs={name:val})
      else:
        print "Unknown parameter "+name+"! Submit unsupported parameters as KEY@NAMELIST in the header!"
    
  if setglacier:
      model.modify(glaciers=glaciers)
  if "landmap" in job.fields:
      lmapname = job.parameters["landmap"]
      tmapname = None
      if "topomap" in job.fields:
        tmapname = job.parameters["topomap"]
        model.modify(landmap=lmapname,topomap=tmapname)
      
  print "Arguments set"
  
  gasesvx = {"H2":0,"He":0,"CO2":0,"N2":0,"O2":0}
  for gas in model.pgases:
      if gas in gasesvx:
          gasesvx[gas] = model.pgases[gas]/model.pressure
  
  transittag = ''
  if maketransit:
      transitparams = '%f %f %f '%(plarad,grav,starrad)
      transitparams+="%f "%(gasesvx['H2']*model.mmw/exo.smws['mH2'])
      transitparams+="%f "%(gasesvx['He']*model.mmw/exo.smws['mHe'])
      transitparams+="%f "%(gasesvx['CO2']*model.mmw/exo.smws['mCO2'])
      transitparams+="%f "%(gasesvx['N2']*model.mmw/exo.smws['mN2'])
      transitparams+="%f "%(gasesvx['O2']*model.mmw/exo.smws['mO2'])
      transitparams+="%f"%model.gascon
      
      
      print "CONFIGURING POST-RUN TRANSIT SPECTROSCOPY..."
      txs = transitparams.split()[3:]
      txsn = ['H2','He','CO2',"N2",'O2']
      for ntt in range(len(txsn)):
          print "\t... %s mass fraction = \t %s"%(txsn[ntt],txs[ntt])
      if len(txs)>len(txsn):
          print "\t... gas constant R manually set to %s"%(txs[-1])
        
      transitjob = Job("# PID MODEL JOBNAME STATE NCORES QUEUE","%s transit transit_%s 0 1 sandyq"%(job.pid,job.name),-1)
      transitfw =(BATCHSCRIPT(transitjob,"abe")+
                   "mkdir "+SCRATCH+"/transit_%s   \n"%job.name+
                   "cp %s/exoplasim/output/%s_snapshot.nc /mnt/node_scratch/"%(job.top,job.name)+USER+"/transit_%s/  \n"%job.name+
                   "cp %s/exoplasim/plasimtransit.py /mnt/node_scratch/"%job.top+USER+"/transit_%s/   \n"%job.name+
                   "cd "+SCRATCH+"/transit_%s   \n"%job.name+
                   "python plasimtransit.py %s_snapshot.nc %s   \n"%(job.name,transitparams)+
                   "cp *.png %s/exoplasim/output/   \n"%job.top+
                   "cp *.pdf %s/exoplasim/output/   \n"%job.top+
                   "cp *.npy %s/exoplasim/output/   \n"%job.top+
                   "cd ../    \n"+
                   "rm -rf transit_%s/        \n"%job.name+
                   "cd $PBS_O_WORKDIR    \n")
      with open(workdir+"/runtransit_%s"%job.name,"w") as transitf:
          transitf.write(transitfw)
      transittag = ("   cp runtransit_%s %s/exoplasim/output/ \n"%(job.name,job.top)+
                    "   cd %s/exoplasim/output/ \n"%job.top+
                    "   qsub runtransit_%s\n"%job.name)
  
  sbdarttag = ''
  if makesbdart:
      
      #ensure that we only output in high-cadence mode if both flags are enabled
      if not highcadence:
          hcout=False
      if not hcout:
          highcadence=False
    
      sbdparams = [job.name,job.ncores]
      if "pCO2" in model.pgases:
          _pco2 = model.pgases['pCO2']*1.0e3
      else:
          _pco2 = 0.0
      sbdparams.append(_pco2)
      sbdparams.append(flux)
      print model.pressure*1.0e3, _pco2
      sbdparams.append(model.pressure*1.0e3-_pco2)
      sbdparams.append(grav)
      sbdparams.append('^'.join(ntimes[1:-1].split(',')))
      sbdparams.append('^'.join(lviews[1:-1].split(',')))
      sbdparams.append(str(highcadence))
      
      print tuple(sbdparams),sbdvar
      os.system("mkdir "+job.top+"/sbdart_%s/%s"%(sbdvar,job.name))
      
      sbd_specfile = ""
      
      if startemp>0:
          wvl_icm = np.arange(100.0,49980.0,20.0) # 20 cm-1 resolution
          wvl_cm = 1.0/wvl_icm #centimeters
          wvl = wvl_cm[::-1]*1.0e4 #microns
          wvl_m = wvl*1.0e-6 #meters
          bf = 2*hh*cc**2/wvl_m**5 * 1.0/(np.exp(hh*cc/(wvl_m*kb*startemp))-1) #Planck function
          Lsun = 4*np.pi*rsun**2 * bf #solar luminosity (W/m)
          Fsun = Lsun / (4*np.pi*au**2) #solar flux at 1 AU (W/m^2/m)
          flux = Fsun * 1367.0/np.trapz(Fsun,x=wvl_m) #adjust to equal 1367 W/m^2
          #It's okay that we're not using user-specified insolation here, since that's set
          #separately in SBDART with solfac--but that's a modifier relative to 1367 W/m^2.
          flux *= 1.0e-6 #convert to W/m^2/um
          spectxt = ''
          for nw in range(len(wvl)-1):
              spectxt += "%s %s \n"%(wvl[nw],flux[nw])
          spectxt += "%s %s"%(wvl[-1],flux[-1])
          with open(job.top+"/sbdart_%s/%s/sbdart_%d.dat"%(sbdvar,job.name,startemp),"w") as specf:
              specf.write(spectxt)
          sbd_specfile = "sbdart_%d.dat"%startemp 
      
      if starfile: #This will always supersede a specified blackbody temperature
          with open(job.top+"/exoplasim/clean/"+starfile,"r") as starf:
              spectxt = starf.read().split('\n')
          if spectxt[-1]=='':
              spectxt = spectxt[:-1]
          with open(job.top+"/sbdart_%s/%s/sbdart_%s"%(sbdvar,job.name,starfile),"w") as specf:
              specf.write('\n'.join(spectxt[1:]))
          sbd_specfile = "sbdart_%s"%starfile
      
      
      sbdtag = sbd_specfile
      
      os.system("cp "+job.top+"/identity.py "+job.top+"/sbdart_%s/"%sbdvar)
      sbdjob = Job("# PID MODEL JOBNAME STATE NCORES QUEUE","%s sbdart_%s sbd_%s 0 %d sandyq"%(job.pid,sbdvar,job.name,job.ncores),-1)
      sbdfw = (BATCHSCRIPT(sbdjob,"abe")+
               "cd "+job.top+"/sbdart_%s    \n"%sbdvar+
               "python buildsbdart.py %s %d sandyq %f %f %f %f %s %s %s "%tuple(sbdparams)+sbdtag+"    \n")
      with open(job.top+"/sbdart_%s/"%sbdvar+job.name+"/autosbdart","w") as sbdf:
          sbdf.write(sbdfw)
      sbdarttag = ("   cd "+job.top+"/sbdart_%s/"%sbdvar+job.name+"/   \n"+
                   "   qsub autosbdart               \n"+
                   "   cp %s/job.npy "%model.workdir+job.top+"/sbdart_%s/"%sbdvar+job.name+"/    \n"+
                   "   cd $PBS_O_WORKDIR             \n")
      os.system("cp "+job.top+"/release.py "+job.top+"/sbdart_%s/"%sbdvar+job.name+"/")

  histargs = ''
  #if weathrestart:
      #histargs = wfile+" "+str(yearini)
  #if nwesteros:
      #histargs = str(nrunyears)
  
  model.save()
  model.exportcfg()
  
  runprefix="./"
  if scriptfile[-3:]==".py":
      runprefix = "python "
      
  
  # You may have to change this part
  jobscript =(BATCHSCRIPT(job,notify)+
              "rm keepgoing                                                     \n"+
              "mkdir "+SCRATCH+"/job"+jid+"            \n")
  jobscript+=("mkdir "+SCRATCH+"/job"+jid+"/snapshots         \n"+
              "mkdir "+SCRATCH+"/job"+jid+"/highcadence         \n"+
                  "tar cvzf stuff.tar.gz --exclude='*_OUT*' --exclude='*_REST*' --exclude='*.nc' --exclude='snapshots/' --exclude='highcadence/' ./* \n")
  jobscript +=("rsync -avzhur stuff.tar.gz "+SCRATCH+"/job"+jid+"/         \n"+
              "rm -rf stuff.tar.gz                     \n"+
              "cd "+SCRATCH+"/job"+jid+"/              \n"+
              "tar xvzf stuff.tar.gz                   \n"+
              "rm stuff.tar.gz          \n"+
              runprefix+scriptfile+" "+str(job.ncores)+" "+str(nlevs)+" "+histargs+"   \n"+
              "rm ice_output ocean_output    \n"+
              "tar cvzf stuff.tar.gz *                                          \n"+
              "rsync -avzhur stuff.tar.gz $PBS_O_WORKDIR/                                          \n"+
              "cd ../ \n"+
              "rm -rf job"+jid+"                                                         \n"+
              "cd $PBS_O_WORKDIR                                                \n"+
              "tar xvzf stuff.tar.gz                                \n"+
              "rm stuff.tar.gz          \n"+
              "if [ -e keepgoing ]                                              \n"+
              "then                                                            \n"+
              "      qsub runexoplasim > %s/%s_new.id                \n"%(job.top,job.pid)+
              "      cd "+job.top+"                           \n"+
              "      python updatedependency.py %s            \n"%job.pid+
              "      cd $PBS_O_WORKDIR                        \n"+
              "else        \n")
  if keeprs:
      jobscript+= "cp plasim_restart ../output/"+job.name+"_restart            \n"
      tag+=" restart"
  if monitor:
      jobscript+= "python monitor_balance.py                                   \n"
      
  jobscript += ("   [ -e balance.log ] && cp balance.log ../output/"+job.name+"_balance.log    \n"+
                "   [ -e slopes.log ] && cp slopes.log ../output/"+job.name+"_slopes.log    \n"+
                "   python synthoutput.py MOST 1                                      \n"+
                "   cp MOST_history.npy ../output/"+job.name+"_history.npy            \n"+
                '   python '+cleanup+' '+job.name+' '+tag+'                           \n'+
                transittag+
                sbdarttag+
                "fi \n")
  
  rs = open(workdir+"/runexoplasim","w")
  rs.write(jobscript)
  rs.close()
  
def submit(job):
  if "DEPENDENCIES" in job.parameters:
     dlist = job.parameters["DEPENDENCIES"].split(',')
     priorjobs = []
     for d in dlist:
         with open(d+".id","r") as f:
             priorjobs.append(f.read().split('\n')[0].split()[0])
         os.system("echo %s >> "%job.pid+d+".id") #indicate that we depend on this job
     os.system("cd exoplasim/job"+str(job.home)+" && "+HOLD(priorjobs)+" runexoplasim > %s/%s.id && cd "%(job.top,job.pid)+job.top)

  else:
     os.system("cd exoplasim/job"+str(job.home)+" && "+SUB+" runexoplasim > %s/%s.id && cd "%(job.top,job.pid)+job.top)
  time.sleep(1.0)
  tag = job.getID()
  job.write()