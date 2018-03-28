import numpy as np
import netCDF4 as nc
import os
import time


def write_input(workdir,cszenith,azimuth,latitude,surface,pCO2,p0,tsurf,altz,flux,clouds=True,abs_sc=True,sic=0.0):
  template =("&INPUT                                                     \n"+ #0
             " IDATM=          4,                                        \n"+ #1
             " AMIX= -1.0000000000000000     ,                           \n"+ #2
             " ISAT=          -4,                                        \n"+ #3
             " WLINF= 10.0000000000000     ,                             \n"+ #4  
             " WLSUP= 4.7250000000000     ,                              \n"+ #5  
             " WLINC=  20  ,                                           \n"+ #6 #set to -0.01 for lower res
             " SZA=  0.0000000000000000     ,                            \n"+ #7
             " CSZA= -1.0000000000000000     ,                           \n"+ #8   --
             " SOLFAC=  1.0000000000000000     ,                         \n"+ #9   --
             " NF=          3,                                           \n"+ #10
             " IDAY=          0,                                         \n"+ #11
             " TIME=  16.000000000000000     ,                           \n"+ #12
             " ALAT= -64.766998291015625     ,                           \n"+ #13
             " ALON= -64.067001342773438     ,                           \n"+ #14
             " ZPRES= -1.0000000000000000     ,                          \n"+ #15  -- (alt)
             " PBAR= -1.0000000000000000     ,                           \n"+ #16  --
             " SCLH2O= -1.0000000000000000     ,                         \n"+ #17
             " UW= -1.0000000000000000     ,                             \n"+ #18
             " UO3= -1.0000000000000000     ,                            \n"+ #19
             " O3TRP= -1.0000000000000000     ,                          \n"+ #20
             " ZTRP=  0.0000000000000000     ,                           \n"+ #21
             " XRSC=  1.0000000000000000     ,                           \n"+ #22  --
             " XN2= -1.0000000000000000     ,                            \n"+ #23  --
             " XO2= -1.0000000000000000     ,                            \n"+ #24  --
             " XCO2= -1.0000000000000000     ,                           \n"+ #25  --
             " XCH4= -1.0000000000000000     ,                           \n"+ #26  --
             " XN2O= -1.0000000000000000     ,                           \n"+ #27  --
             " XCO= -1.0000000000000000     ,                            \n"+ #28  --
             " XNO2= -1.0000000000000000     ,                           \n"+ #29  --
             " XSO2= -1.0000000000000000     ,                           \n"+ #30  --
             " XNH3= -1.0000000000000000     ,                           \n"+ #31  --
             " XNO= -1.0000000000000000     ,                            \n"+ #32  --
             " XHNO3= -1.0000000000000000     ,                          \n"+ #33  --
             " XO4=  1.0000000000000000     ,                            \n"+ #34  --
             " ISALB=          4,                                        \n"+ #35  --
             " ALBCON=  0.0000000000000000     ,                         \n"+ #36
             " SC= 5*3.4028234663852886E+038  ,                          \n"+ #37  --
             " ZCLOUD= 5*0.0000000000000000       ,                      \n"+ #38
             " TCLOUD= 5*0.0000000000000000       ,                      \n"+ #39
             " LWP= 5*0.0000000000000000       ,                         \n"+ #40
             " NRE= 5*8.0000000000000000       ,                         \n"+ #41  --
             " RHCLD= -1.0000000000000000     ,                          \n"+ #42
             " KRHCLR=          0,                                       \n"+ #43
             " JAER= 5*0          ,                                      \n"+ #44
             " ZAER= 5*0.0000000000000000       ,                        \n"+ #45
             " TAERST= 5*0.0000000000000000       ,                      \n"+ #46
             " IAER=          0,                                         \n"+ #47
             " VIS= -1.0000000000000000     ,                            \n"+ #48
             " RHAER= -1.0000000000000000     ,                          \n"+ #49
             " TBAER= -1.0000000000000000     ,                          \n"+ #50
             " WLBAER= 150*-1.0000000000000000      ,                    \n"+ #51
             " QBAER= 150*-1.0000000000000000      ,                     \n"+ #52
             " ABAER=  0.0000000000000000     ,                          \n"+ #53
             " WBAER= 150*-1.0000000000000000      ,                     \n"+ #54
             " GBAER= 150*-1.0000000000000000      ,                     \n"+ #55
             " PMAER= 44850*-1.0000000000000000      ,                   \n"+ #56
             " ZBAER= 65*-1.0000000000000000      ,                      \n"+ #57
             " DBAER= 65*-1.0000000000000000      ,                      \n"+ #58
             " NOTHRM=         -1,                                       \n"+ #59
             " NOSCT=          0,                                        \n"+ #60
             " KDIST=          3,                                        \n"+ #61
             " ZGRID1=  1.0000000000000000     ,                         \n"+ #62
             " ZGRID2=  30.000000000000000     ,                         \n"+ #63
             " NGRID=          0,                                        \n"+ #64
             " IDB= 20*0          ,                                      \n"+ #65
             " ZOUT=  0.0000000000000000     ,  100.00000000000000     , \n"+ #66
             " IOUT=         5,                                          \n"+ #67
             " PRNT= 7*F,                                                \n"+ #68
             " TEMIS=  0.0000000000000000     ,                          \n"+ #69
             " NSTR=          0,                                         \n"+ #70
             " NZEN=          0,                                         \n"+ #71   --
             " UZEN=        0.0,84.375,                                  \n"+ #72   --
             " NPHI=          0,                                         \n"+ #73   --
             " PHI=          90,                                         \n"+ #74   --
             " SAZA=  180.00000000000000     ,                           \n"+ #75   --
             " IMOMC=          3,                                        \n"+ #76
             " IMOMA=          3,                                        \n"+ #77
             " TTEMP= -1.0000000000000000     ,                          \n"+ #78
             " BTEMP= 285.501     ,                                      \n"+ #79   --
             " CORINT=F,                                                 \n"+ #80
             " SPOWDER=F,                                                \n"+ #81
             " /                                                         \n")
# We have a problem in that if we consider polar regions, and us (the observer) to be a star near the 
# equator of the planet's celestial sphere, our star should never rise a few degrees above the horizon
# near the poles. And this means that the zenith angle will only reach 0 degrees at the equator. The 
# unfortunate consequence of this is that we cannot simply specify 2 azimuth angles (E & W). At any
# latitude other than the equator, the azimuth angle must also vary with zenith angle. But this is 
# a problem, because SBDART outputs NPHI*NZEN angles--we desire 16 zenith angles and 32 azimuth 
# angles (to account for forward-scattering), for 32 viewing angles total. But giving SBDART both
# lists of angles will give us 512 viewing angles, which is atrocious. We could run SBDART once for
# every viewing angle, resulting in 65,536 SBART runs per planet. But that only runs 3.5 times as fast,
# and we must run SBART 32 times more often. This is more than a factor of 10 increase in computation
# time, and therefore untenable. Doing it the brute force way, however, will likely cost us a factor of 
# 10 or so as well. It is therefore likely that we will have to modify SBDART itself if possible, to 
# accept a linear series of viewing angles rather than a regular az-zen grid.
#
# SBDART has radiance angles hard-coded as a regular grid. Not only that, but it's not SBDART that did
# that--it's DISORT. Remedying this will involve a not-insignificant overhaul of SBDART/DISORT, meriting
# a fork on github and a new name--perhaps SBDARTx, for SBDART Extended. Or DISORTx, along the same lines.
# Given that this bit of programming folly on the part of the DISORT programmers turns an O(N) operation
# into O(N^2), this may well be worth it, although it may set us back at least 2 weeks. Hrumph. So it goes.
  
 
  #Default mixing ratios
  xn2   = 781000.00 * 1.0e-6
  xo2   = 209000.00 * 1.0e-6
  xco2  = 360.00    * 1.0e-6
  xch4  = 1.74      * 1.0e-6
  xn2o  = 0.32      * 1.0e-6
  xco   = 0.15      * 1.0e-6
  xnh3  = 5.0e-4    * 1.0e-6
  xso2  = 3.0e-4    * 1.0e-6
  xno   = 3.0e-4    * 1.0e-6
  xhno3 = 5.0e-5    * 1.0e-6
  xno2  = 2.3e-5    * 1.0e-6
 
  pn2  = xn2   * p0 #mbars (hPa)
  po2  = xo2   * p0
  pch4 = xch4  * p0
  pn2o = xn2o  * p0
  pco  = xco   * p0
  pnh3 = xnh3  * p0
  pso2 = xso2  * p0
  pno  = xno   * p0
  phno3= xhno3 * p0
  pno2 = xno2  * p0
  
  # Ground types
  #-1 user
  # 0 uniform
  # 1 snow
  # 2 clearwater
  # 3 lakewater
  # 4 seawater
  # 5 sand
  # 6 vegetation
  
  input_text = template.split('\n')
  
  if cszenith <= 0.0: #There is no incident sunlight, so we only need the thermal range >2.5 microns
    input_text[4] = " WLINF= 10.972857000  ,    "   # as opposed to >0.55 microns. This saves us time.
    input_text[5] = " WLSUP= 4.23857210000 ,    "
                #The spectral binning if we do this isn't perfectly aligned with the full range, so
                #we'll need to devise a binning scheme for putting them all in the right bins.
  
#UZEN
  rlat = latitude*np.pi/180.0
  hours = np.linspace(0,84.375,num=16)*np.pi/180.0
  cuz = np.cos(rlat)*np.cos(hours)
  uz = np.arccos(cuz)*180.0/np.pi
  suz = np.sin(uz*np.pi/180.0)
  
#Azimuths
  cosazs = -cuz*np.sin(rlat) / (suz*np.cos(rlat))
  cosazs[np.isnan(cosazs)] = 0.0
  uazs = np.arccos(np.minimum(np.maximum(cosazs,-1.0),1.0)) * 180.0/np.pi
  
  avguaz = np.mean(uazs)
  
  input_text[72] = " UZEN= "
  for zz in uz:
    input_text[72]+= "%.16f,"%zz
  input_text[72] = input_text[72][:-1] #Strip off that last comma
  input_text[72]+="   , "
  
  input_text[74] = " PHI= %.16f,%.16f  , "%(avguaz,360.0-avguaz)
  
  input_text[8] = " CSZA= %.16f   ,     "%cszenith
  input_text[75] = " SAZA= %3.16f   ,   "%azimuth
  input_text[79] = " BTEMP= %3.3f   ,   "%tsurf
  input_text[9] = " SOLFAC= %.16f   ,   "%(flux/1367.0)
  input_text[15] = " ZPRES= %.16f   ,   "%altz
  if not abs_sc:
    input_text[16] = " PBAR= %.16f   ,   "%(0.0)
  if clouds:
    input_text[41] = " NRE = 5*%.16f  ,  "%(0.0)
    
  psurf = p0 + pCO2
  
  xco2  = pCO2  / psurf * 1.0e6
  xn2   = pn2   / psurf * 1.0e6
  xo2   = po2   / psurf * 1.0e6
  xch4  = pch4  / psurf * 1.0e6
  xn2o  = pn2o  / psurf * 1.0e6
  xco   = pco   / psurf * 1.0e6
  xnh3  = pnh3  / psurf * 1.0e6
  xso2  = pso2  / psurf * 1.0e6
  xno   = pno   / psurf * 1.0e6
  xhno3 = phno3 / psurf * 1.0e6
  xno2  = pno2  / psurf * 1.0e6
  
  input_text[23] = " XN2=   %.16f   ,  "%xn2 
  input_text[24] = " XO2=   %.16f   ,  "%xo2  
  input_text[25] = " XCO2=  %.16f   ,  "%xco2  
  input_text[26] = " XCH4=  %.16f   ,  "%xch4 
  input_text[27] = " XN2O=  %.16f   ,  "%xn2o 
  input_text[28] = " XCO=   %.16f   ,  "%xco  
  input_text[29] = " XNO2=  %.16f   ,  "%xnh3 
  input_text[30] = " XSO2=  %.16f   ,  "%xso2 
  input_text[31] = " XNH3=  %.16f   ,  "%xno  
  input_text[32] = " XNO=   %.16f   ,  "%xhno3
  input_text[33] = " XHNO3= %.16f   ,  "%xno2 
  
  if surface == "uniform":
    input_text[35] = " ISALB= 0  , "
  elif surface == "snow":
    input_text[35] = " ISALB= 1  , "
  elif surface == "clearwater":
    input_text[35] = " ISALB= 2  , "
  elif surface == "lakewater":
    input_text[35] = " ISALB= 3  , "
  elif surface == "seawater":
    input_text[35] = " ISALB= 4  , "
  elif surface == "sand":
    input_text[35] = " ISALB= 5  , "
  elif surface == "vegetation":
    input_text[35] = " ISALB= 6  , "
  elif surface == "seamix":
    input_text[35] = " ISALB= 10 , "
    input_text[37] = " SC= %.10f,%.10f,0.000,0.000   , "%(sic,1-sic)
  else:
    input_text[35] = " ISALB= -1  , "
    
  inputtext = '\n'.join(input_text)
  
  inp = open(workdir+"/INPUT","w")
  inp.write(inputtext)
  inp.close()
  
  
def writecolumn(z,p,t,q,o,lon,lat,workdir,name='atms.dat'):
    dat = ''
    nlev = len(z[:,lat,lon])
    dat+=str(nlev)+'\n'
    for k in range(0,nlev):
        layer = []
        layer.append(str(z[k,lat,lon]))
        layer.append(str(p[k,lat,lon]))
        layer.append(str(t[k,lat,lon]))
        layer.append(str(q[k,lat,lon]))
        layer.append(str(o[k,lat,lon]))
        line = ' '.join(layer)
        dat+=line+'\n'
    f=open(workdir+"/"+name,"w")
    f.write(dat)
    f.close()  

    
def cloudcolumn(dql,cld,lon,lat,workdir,name='usrcld.dat'):
    dat = ''
    nlev = len(cld[:,lat,lon])
    for k in range(nlev-1,-1,-1):
        layer = []
        layer.append(str(dql[k,lat,lon]))
        layer.append("8")
        layer.append("-1")
        layer.append("-1")
        layer.append(str(cld[k,lat,lon]))
        line = ' '.join(layer)
        dat+=line+'\n'
    f=open(workdir+'/'+name,"w")
    f.write(dat)
    f.close()
  
def writecolumn_single_plasim(z,p,t,q,o,workdir,name='atms.dat'):
    dat = ''
    nlev = len(z)
    dat+=str(nlev)+'\n'
    for k in range(0,nlev):
        layer = []
        layer.append(str(z[k]))
        layer.append(str(p[k]))
        layer.append(str(t[k]))
        layer.append(str(q[k]))
        layer.append(str(o[k]))
        line = ' '.join(layer)
        dat+=line+'\n'
    f=open(workdir+'/'+name,"w")
    f.write(dat)
    f.close()
    
def writecolumn_single_lmdz(z,p,t,q,o,workdir,name='atms.dat'):
    dat = ''
    nlev = len(z)
    dat+=str(nlev)+'\n'
    for k in range(nlev-1,-1,-1):
        layer = []
        layer.append(str(z[k]*1.0e-3))
        layer.append(str(p[k]*0.01))
        layer.append(str(t[k]))
        layer.append(str(q[k]))
        layer.append(str(o[k]))
        line = ' '.join(layer)
        dat+=line+'\n'
    f=open(workdir+'/'+name,"w")
    f.write(dat)
    f.close()    
    
def cloudcolumn_single(dql,cld,workdir,name='usrcld.dat'):
    dat = ''
    nlev = len(cld)
    for k in range(nlev-1,-1,-1):
        layer = []
        layer.append(str(dql[k]))
        layer.append("8")
        layer.append("-1")
        layer.append("-1")
        layer.append(str(cld[k]))
        line = ' '.join(layer)
        dat+=line+'\n'
    f=open(workdir+'/'+name,"w")
    f.write(dat)
    f.close()  
    
def cloudcolumn_single_lmdz(dql,dqi,rel,rei,cld,workdir,name='usrcld.dat'):
    dat = ''
    nlev = len(cld)
    for k in range(0,nlev):
        layer = []
        layer.append(str(dql[k])) #Liquid water path (g/m^2)
        layer.append(str(rel[k])) #Water droplet radius (um)
        layer.append(str(dqi[k])) #Frozen water path (g/m^2)
        layer.append(str(rei[k])) #Ice particle radius (um)
        layer.append(str(cld[k])) #Cloud fraction
        line = ' '.join(layer)
        dat+=line+'\n'
    f=open(workdir+'/'+name,"w")
    f.write(dat)
    f.close()  
    
def getalt(ta,lev,nlats=32,nlons=64,grav=9.80665,gascon=287.0):
    zdh = np.zeros((nlats,nlons))
    zzh = np.zeros(zdh.shape)
    zzf = np.zeros((len(lev),nlats,nlons))
    
    i=len(lev)-1
    while i>0:
        zdh[:] = -ta[i,:,:]*gascon/grav*np.log(lev[i-1]/lev[i])
        zzf[i,:,:] = zzh + 0.5*zdh
        zzh[:] += zdh
        i-=1
    zdh[:] = -ta[0,:,:]*gascon/grav*np.log(0.5)*0.5
    zzf[0,:,:] = zzh + zdh
    
    return zzf
  
def getalt_single(ta,lev,grav=9.80665,gascon=287.0):
    zdh = 0.0
    zzh = 0.0
    zzf = np.zeros(len(lev))
    
    i=len(lev)-1
    while i>0:
        zdh = -ta[i]*gascon/grav*np.log(lev[i-1]/lev[i])
        zzf[i] = zzh + 0.5*zdh
        zzh += zdh
        i-=1
    zdh = -ta[0]*gascon/grav*np.log(0.5)*0.5
    zzf[0] = zzh + zdh
    
    return zzf  


def analyzecell_plasim(data,lat,lon,workdir,grav=9.80665):
  #cszenith,azimuth,surface,pCO2,p0,tsurf,altz
  
  surf = "uniform"
  
  lsm = data.variables['lsm'][-1,lat,lon]
  if lsm < 0.5: #sea
    sic = data.variables['sic'][-1,lat,lon]
    if sic == 1.0:
      surf = "snow"
    elif sic == 0.0:
      surf = "seawater"
    else:
      surf = "seamix"
  else:
    surf = "sand"
    #We could put in an actual algorithm to get albedo given soil wetness, or combine sand, veg, etc
    
  tsurf = np.mean(data.variables['ts'][:,lat,lon])
  
  altz = data.variables["netz"][-1,lat,lon] / grav #meters
  
  ta = np.mean(data.variables["ta"][:,:,lat,lon],axis=0)
  
  lvs = data.variables['lev']
  #lns = data.variables['lon']
  #lts = data.variables['lat']
  
  zs = getalt_single(ta,lvs,grav=grav,gascon=287.0)
  
  ps = np.mean(data.variables['ps'][:,lat,lon])
  pa = ps * lvs
  
  gascon0 = 8.3144598
  mmair = 0.028964
  mmvap = 0.018016
  
  satvap = 610.78*10.0**(7.5*(ta-273.15)/(ta+0.15))
  
  relhum = np.mean(data.variables['hur'][:,:,lat,lon],axis=0)
  
  pvap = relhum*satvap
  
  rvap = 461.495
  rdry = 287.058
  
  rhohum = ((pa*100.0-pvap)*mmair+pvap*mmvap)/(gascon0*ta)
  
  hus = np.mean(data.variables['hus'][:,:,lat,lon],axis=0)
  rhoh2o = hus*rhohum
  
  print "Writing atms.dat"
  #write atms.dat, which contains the atmosphere profile--height, pressure, temp, water vapor, and ozone
  writecolumn_single_lmdz(zs,pa,ta,rhoh2o*1000.0,np.zeros(10),workdir)
  
  dsigma = np.zeros(len(lvs))
  dsigma[0] = lvs[0]
  for i in range(1,len(lvs)):
    dsigma[i] = lvs[i]-lvs[i-1]
    
  clw = np.mean(data.variables['clw'][:,:,lat,lon],axis=0)
  
  dql = np.zeros(len(lvs))
  
  dql = dsigma*1000.0*ps*clw/grav #g/m^2 #Cloud liquid water path
  
  cld = np.mean(data.variables['cl'][:,:,lat,lon],axis=0)
  
  print "Writing usrcld.dat"
  #Write usrcld.dat, which has level data on cloud water content and coverage fraction
  cloudcolumn_single(dql,cld,workdir)
  
  csz = data.variables["czen"][-1,lat,lon]
  
  nlcsz = data.variables['czen'][-1,lat,lon-1]
  
  if nlcsz>csz:
    direction = 'W' #sun is west of us, since it's brighter to the west
  else:
    direction = 'E' #sun is east of us, since it's brighter to the east
    
  ssz = np.sin(np.arccos(csz))
  clat = data.variables['lat'][lat] * np.pi/180.0
  
  cosaz = -csz*np.sin(clat) / (ssz*np.cos(clat))
  if np.isnan(cosaz):
    cosaz = 0.0
  if np.isinf(cosaz):
    cosaz = 0.0
  #print cosaz
  azm = np.arccos(np.minimum(np.maximum(cosaz,-1.0),1.0)) * 180.0/np.pi
  
  if direction=='W':
    azm = 360.0 - azm #Make it between 180 and 360
  
  return csz,azm,surf,sic,tsurf,altz
    
def cosine_zenith(latitude,longitude,sol_dec,sol_lon):
  sinlat = np.sin(latitude*np.pi/180.)
  sindec = np.sin(sol_dec*np.pi/180.)
  coslat = np.cos(latitude*np.pi/180.)
  cosdec = np.cos(sol_dec*np.pi/180.)
  coslon = np.cos((longitude-sol_lon)*np.pi/180.0)
  coszen = sinlat*sindec + coslat*cosdec*coslon
  
  #coszen[abs(longitude-sol_lon)>90] = 0.0
  
  return coszen

def analyzecell_lmdz(data,lat,lon,workdir,grav=9.80665,sol_dec=0.0,sol_lon=0.0):
  #cszenith,azimuth,surface,pCO2,p0,tsurf,altz
  
  surf = "uniform"
  
  #lsm = data['lsm'][-1,lat,lon]
  lsm = 0.0
  if lsm < 0.5: #sea
    sic = data['pctsrf_sic'][lat,lon]
    if sic == 1.0:
      surf = "snow"
    elif sic == 0.0:
      surf = "seawater"
    else:
      surf = "seamix"
  else:
    surf = "sand"
    #We could put in an actual algorithm to get albedo given soil wetness, or combine sand, veg, etc
    
  tsurf = data['tsurf'][lat,lon]
  
  zs = data["altitude"]*1000.0 #meters
  
  ta = data["temp"][:,lat,lon]
  
  #altz = data["phisinit"][lat,lon]
  altz=0.0
  
  #lvs = data['lev']
  #lns = data['lon']
  #lts = data['lat']
  
  #zs = getalt_single(ta,lvs,grav=grav,gascon=287.0)
  
  ps = data['ps'][lat,lon]
  pa = data['p'][:,lat,lon]
  
  gascon0 = 8.3144598
  mmair = 0.028964
  mmvap = 0.018016
  
  satvap = 610.78*10.0**(7.5*(ta-273.15)/(ta+0.15))
  
  relhum = data['RH'][:,lat,lon]
  
  pvap = relhum*satvap
  
  rvap = 461.495
  rdry = 287.058
  
  rhohum = ((pa*100.0-pvap)*mmair+pvap*mmvap)/(gascon0*ta)
  
  hus = np.maximum(data['h2o_vap'][:,lat,lon],0.0) #kg/kg?
  rhoh2o = hus*rhohum
  
  print "Writing atms.dat"
  #write atms.dat, which contains the atmosphere profile--height, pressure, temp, water vapor, and ozone
  writecolumn_single_lmdz(zs,pa,ta,rhoh2o*1000.0,np.zeros(len(zs)),workdir)
  
  #dpress = np.zeros(len(zs))
  #dpress[0] = ps-pa[0]*100.0
  #for i in range(1,len(zs)):
    #dpress[i] = pa[i-1]-pa[i]
    
  dpress = abs(np.diff(data['press_inter']))  
    
  #dsigma = np.diff(data['bps'])
  
  #dpress = np.diff(data['press_inter'][:,lat,lon])*100.0 #Pa
  
  clw = hus*1000 #g/kg
  
  dql = np.maximum(dpress*clw/grav,0.0) #g/m^2 #Water mass in cell
  
  cld = data['CLF'][:,lat,lon]
  
  dqi = np.maximum(dpress*data["h2o_ice"][:,lat,lon]*1000/grav,0.0)
  
  rei = np.minimum(data["H2Oice_reff"][:,lat,lon]*1.0e6,128.0) #um
  
  print "Writing usrcld.dat"
  #Write usrcld.dat, which has level data on cloud water content and coverage fraction
  cloudcolumn_single_lmdz(dql,dqi,np.zeros(len(cld))+8,rei,cld,workdir)
  
  latitude = data['latitude'][lat]
  longitude = data['longitude'][lon]
  plon = data['longitude'][lon-1]
  
  csz = cosine_zenith(latitude,longitude,sol_dec,sol_lon)  #still need--may have to compute by hand
  
  nlcsz = cosine_zenith(latitude,plon,sol_dec,sol_lon)
  
  if nlcsz>csz:
    direction = 'W' #sun is west of us, since it's brighter to the west
  else:
    direction = 'E' #sun is east of us, since it's brighter to the east
    
  ssz = np.sin(np.arccos(csz))
  clat = latitude * np.pi/180.0
  
  cosaz = -csz*np.sin(clat) / (ssz*np.cos(clat))
  if np.isnan(cosaz):
    cosaz = 0.0
  if np.isinf(cosaz):
    cosaz = 0.0
  #print cosaz
  azm = np.arccos(np.minimum(np.maximum(cosaz,-1.0),1.0)) * 180.0/np.pi
  
  if direction=='W':
    azm = 360.0 - azm #Make it between 180 and 360
  
  return csz,azm,surf,sic,tsurf,altz
    
def prep(job):
  if "type" not in job.parameters:
      print "Warning: need to specify which GCM produced this data!"
  else:
      if job.parameters["type"]=="plasim":
          _prep_plasim(job)
      elif job.parameters["type"]=="lmdz":
          _prep_lmdz(job)
      elif job.parameters["type"]=="fms":
          print "FMS-SBDART interface not yet implemented."
      else:
          print "Model not recognized."
    
def _prep_lmdz(job):    
  workdir = "sbdart/job"+str(job.home)
  if "source" in job.parameters:
    source = job.parameters["source"]
  else:
    source = "clean"
    
  if "lats" in job.parameters:
    lats = job.parameters["lats"].split(',')
    lats[0] = int(lats[0])
    if len(lats)>1:
      lats[1] = int(lats[1])
    else:
      lats.append(lats[0])
  else:
    lats = [0,48]
  
  if "lons" in job.parameters:
    lons = job.parameters["lons"].split(',')
    lons[0] = int(lons[0])
    if len(lons)>1:
      lons[1] = int(lons[1])
    else:
      lons.append(lons[0])
  else:
    lons = [0,64]
    
  data = np.load("hopper/"+job.parameters["gcm"]).item()  
  
  if "pCO2" in job.parameters:
    pCO2 = float(job.parameters["pCO2"])
  else:
    pCO2 = 0.360
  
  if "p0" in job.parameters:
    p0 = float(job.parameters["p0"])
  else:
    p0 = 1011.0
    
  if "flux" in job.parameters:
    flux = float(job.parameters["flux"])
  else:
    flux = 1367.0
    
  if "grav" in job.parameters:
    grav = float(job.parameters["grav"])
  else:
    grav = 9.80665
    
  if "outhopper" in job.parameters:
    dest = "../"+job.parameters["outhopper"]
    os.system("mkdir sbdart/"+job.parameters["outhopper"])
  else:
    dest = "../output"
    
  nlats = len(data['latitude'][:])
  for jlat in range(lats[0],lats[1]):
    for jlon in range(lons[0],lons[1]):
      print "Lat %02d Lon %02d"%(jlat,jlon)
      os.system("mkdir "+workdir+"/sbdart-%02d_%02d"%(jlat,jlon))
      os.system("cp -r sbdart/"+source+"/* "+workdir+"/sbdart-%02d_%02d/"%(jlat,jlon))
  
  jobscript =("#!/bin/bash -l                                                  \n"+
              "#PBS -l nodes=1:ppn=1                                            \n"+
              "#PBS -q workq                                                    \n"+
              "#PBS -r n                                                        \n"+
              "#PBS -l walltime=48:00:00                                        \n"+
              "#PBS -N "+job.name+"                                             \n"
              "# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE"+
              " nodes,ppn,walltime and my_job_name VALUES                       \n"+
              "cd $PBS_O_WORKDIR                                                \n"+
              "module unload gcc/4.9.1                                          \n"+
              "module unload python/2.7.9                                       \n"+
              "module load intel/intel-17                                       \n"+
              "module load openmpi/2.0.1-intel-17                               \n"+
              "module load python                                              \n"+
              "for jl in {%02d..%02d};                                  \n"%(lats[0],lats[1]-1)+
              "do \n"+
              "     for il in {%02d..%02d};                        \n"%(lons[0],lons[1]-1)+
              "     do \n"+
              "          ILAT=`printf '%02d' $(( 10#$jl ))`           \n"+
              "          ILON=`printf '%02d' $(( 10#$il ))`           \n"+
              "          echo $ILAT $ILON              \n"+
              "          TAG=${ILAT}_${ILON}                    \n"+
              "          cd sbdart-$TAG                          \n"+
              "          ./sbdart > ../sbout.$TAG                \n"+
              "          cd ../                                  \n"+
              "     done                                    \n"+
              "done \n"+
              './release.sh "'+dest+'"                                \n')
  
  rs = open(workdir+"/runsbdart","w")
  rs.write(jobscript)
  rs.close()
  
  os.system("cp -r sbdart/release.sh "+workdir+"/")
  
  for jlon in range(lons[0],lons[1]):
    for jlat in range(lats[0],lats[1]):
      csz,azm,surf,sic,tsurf,altz = analyzecell_lmdz(data,jlat,jlon,
                                                     workdir+"/sbdart-%02d_%02d"%(jlat,jlon),
                                                     grav=grav)
      latitude = data['latitude'][jlat]
      write_input(workdir+"/sbdart-%02d_%02d"%(jlat,jlon),csz,azm,latitude,surf,pCO2,
                  p0,tsurf,altz,flux,clouds=True,abs_sc=True,sic=sic)
      print "Prepped lat %02d lon %02d"%(jlat,jlon)

  
def _prep_plasim(job): #data,lats,lons,pCO2,p0,flux,grav=9.80665
  workdir = "sbdart/job"+str(job.home)
  if "source" in job.parameters:
    source = job.parameters["source"]
  else:
    source = "clean"
    
  if "lats" in job.parameters:
    lats = job.parameters["lats"].split(',')
    lats[0] = int(lats[0])
    if len(lats)>1:
      lats[1] = int(lats[1])
    else:
      lats.append(lats[0])
  else:
    lats = [0,32]
  
  if "lons" in job.parameters:
    lons = job.parameters["lons"].split(',')
    lons[0] = int(lons[0])
    if len(lons)>1:
      lons[1] = int(lons[1])
    else:
      lons.append(lons[0])
  else:
    lons = [0,64]
    
  data = nc.Dataset(job.parameters["gcm"],"r")  
  
  if "pCO2" in job.parameters:
    pCO2 = float(job.parameters["pCO2"])
  else:
    pCO2 = 0.360
  
  if "p0" in job.parameters:
    p0 = float(job.parameters["p0"])
  else:
    p0 = 1011.0
    
  if "flux" in job.parameters:
    flux = float(job.parameters["flux"])
  else:
    flux = 1367.0
    
  if "grav" in job.parameters:
    grav = float(job.parameters["grav"])
  else:
    grav = 9.80665
    
  nlats = len(data.variables['lat'][:])
  for jlat in range(lats[0],lats[1]):
    for jlon in range(lons[0],lons[1]):
      print "Lat %02d Lon %02d"%(jlat,jlon)
      os.system("mkdir sbdart/sbdart-%02d_%02d"%(jlat,jlon))
      os.system("cp -r sbdart/"+source+"/* sbdart/sbdart-%02d_%02d/"%(jlat,jlon))
  
  jobscript =("#!/bin/bash -l                                                  \n"+
              "#PBS -l nodes=1:ppn=1                                            \n"+
              "#PBS -q workq                                                    \n"+
              "#PBS -r n                                                        \n"+
              "#PBS -l walltime=48:00:00                                        \n"+
              "#PBS -N "+job.name+"                                             \n"
              "# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE"+
              " nodes,ppn,walltime and my_job_name VALUES                       \n"+
              "cd $PBS_O_WORKDIR                                                \n"+
              "module unload gcc/4.8.0                                          \n"+
              "module unload openmpi/1.6.4-gcc-4.8.0                            \n"+
              "module load intel/intel-17                                       \n"+
              "module load openmpi/2.0.1-intel-17                               \n"+
              "for jl in {%2d..%2d};                                  \n"%(lats[0],lats[1])+
              "do \n"+
              "     for il in {%2d..%2d};                        \n"%(lons[0],lons[1])+
              "     do \n"+
              "          ILAT=`printf '%02d' $(( 10#$jl ))`           \n"+
              "          ILON=`printf '%02d' $(( 10#$il ))`           \n"+
              "          TAG=${ILAT}_${ILON}                    \n"+
              "          cd sbdart-$TAG                          \n"+
              "          ./sbdart > ../sbout.$TAG                \n"+
              "          cd ../                                  \n"+
              "     done                                    \n"+
              "done \n")
  
  rs = open(workdir+"/runsbdart","w")
  rs.write(jobscript)
  rs.close()
  
  for jlon in range(lons[0],lons[1]):
    for jlat in range(lats[0],lats[1]):
      csz,azm,surf,sic,tsurf,altz = analyzecell(data,jlat,jlon,workdir,grav=grav)
      latitude = data.variables['lat'][jlat]
      write_input(workdir,csz,azm,latitude,surf,pCO2,p0,tsurf,altz,flux,clouds=True,abs_sc=True,sic=sic)
      print "Prepped lat %02d lon %02d"%(jlat,jlon)

def submit(job):
  workdir = "sbdart/job"+str(job.home)
  
  os.system("cd "+workdir+" && qsub runsbdart && cd ../../")
