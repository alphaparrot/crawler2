import numpy as np
import netCDF4 as nc
import os, sys, glob
import time
#from identity import ide.USER, ide.SCRATCH
import identity as ide
#from batch_system import bts.BATCHSCRIPT, bts.SUB, bts.HOLD
import batch_system as bts
#from crawldefs import crd.Job
import crawldefs as crd

def write_input_locked(workdir,nview,ntime,cszenith,azimuth,latitude,longitude,surface,pCO2,p0,
                tsurf,altz,flux,wmin=0.55,wmax=19.0,albedo=0.35,smooth=False,flat=True,sic=0.0,
                spec=False,iout=5,zout=(0,100),icefile='seaice.dat',waterfile='seawater.dat'):
  template =("&INPUT                                                     \n"+ #0
             " IDATM=          0,                                        \n"+ #1
             " AMIX= -1.0000000000000000     ,                           \n"+ #2
             " ISAT=          -4,                                        \n"+ #3
             " WLINF= 10.0000000000000     ,                             \n"+ #4  
             " WLSUP= 4.7250000000000     ,                              \n"+ #5  
             " WLINC=  20  ,                                           \n"+ #6 #set to -0.01 for lower res
             " SZA=  0.0000000000000000     ,                            \n"+ #7
             " CSZA= -1.0000000000000000     ,                           \n"+ #8   --
             " SOLFAC=  1.0000000000000000     ,                         \n"+ #9   --
             " NF=          3,                                           \n"+ #10  -- Star type
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
             " ALBCON=  0.0000000000000000     ,                         \n"+ #36  --
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
             " ZOUT=  %3.8f     ,  %3.8f     , \n"%zout+ #66
             " IOUT=         "+str(iout)+",                              \n"+ #67   
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
  
  wc = 0.5*(wmax+wmin)
  input_text[4] = " WLINF= %2.9f ,   "%wc
  input_text[5] = " WLSUP= %2.9f ,   "%(0.5*(wc-wmin))
  
  if cszenith <= 0.0: #There is no incident sunlight, so we only need the thermal range >2.5 microns
    input_text[4] = " WLINF= 10.972857000  ,    "   # as opposed to >0.55 microns. This saves us time.
    input_text[5] = " WLSUP= 4.23857210000 ,    "
                #The spectral binning if we do this isn't perfectly aligned with the full range, so
                #we'll need to devise a binning scheme for putting them all in the right bins.
#STAR TYPE
  if spec:
      input_text[10] = " NF=      -1,    "
    
#UZEN
#For Earth-like cases, we want to see the full phase, and quarter phases for East, West, North,
#and South. We also want to consider *four* solar zeniths: mostly ocean, mostly land, and two that
#are mixed. 

  rlat = latitude
  #rlat = abs(np.array([rlat,rlat,rlat,90.0-rlat,rlat+90.0]))
  #lcoeff = np.array([1,1,1,-1,1])
  loff = np.array([0,0,0,90,-90])
  rdec = loff[nview]*np.pi/180.0
  #rlat = lcoeff[nview]*rlat + loff[nview]
  #hours = np.linspace(0,84.375,num=16)*np.pi/180.0
  # cos(p)cos(l) = cos(z)
  
  faces = [0.0,90.0,180.0,270.0]
  
  # 12pm, E, W, N, S
  hours = longitude - faces[ntime] + np.array([0,-90,90,0,0])[nview]
  if hours<0:
      hours+=360
  elif hours>360:
      hours-=360
  #hours[hours<0]+=360
  #hours[hours>360] -=360
  
  rlat*=(np.pi/180.0)
  hours*=(np.pi/180.0)
  
  cuz = np.sin(rlat)*np.sin(rdec)+np.cos(rdec)*np.cos(rlat)*np.cos(hours)
  uz = np.arccos(cuz)*180.0/np.pi
  suz = np.sin(uz*np.pi/180.0)
  
#Azimuths
  # if longitude = pi, and lat = 0, cos(phi) is NaN.
  cosazs = (np.sin(rdec)-cuz*np.sin(rlat)) / (suz*np.cos(rlat))
  # -cos(z)sin(p)/(sin(z)cos(p)) = -tan(p)/tan(z)
  #cosazs[np.isnan(cosazs)] = 0.0
  if np.isnan(cosazs):
      cosazs = 0.0
  uazs = np.arccos(np.minimum(np.maximum(cosazs,-1.0),1.0)) * 180.0/np.pi
  
  avguaz = np.mean(uazs)
  
  input_text[72] = " UZEN= "
  input_text[72]+= "%.16f"%uz
  input_text[72]+="   , "
  
  input_text[74] = " PHI=  "
  input_text[74]+= "%.16f"%uazs
  input_text[74]+="   , "
  
  input_text[8] = " CSZA= %.16f   ,     "%cszenith
  if cszenith < 0:
      input_text[7] = " SZA= %.16f    ,     "%(np.arccos(cszenith)*180.0/np.pi)
      input_text[8] = " CSZA= -1     ,      "
  input_text[75] = " SAZA= %3.16f   ,   "%azimuth
  input_text[79] = " BTEMP= %3.3f   ,   "%tsurf
  input_text[9] = " SOLFAC= %.16f   ,   "%(flux/1367.0)
  input_text[15] = " ZPRES= %.16f   ,   "%altz
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
  
  if smooth:
    input_text[23] = " XN2=   %.16f   ,  "%1.0e6 
    input_text[24] = " XO2=   %.16f   ,  "%0.0
    input_text[25] = " XCO2=  %.16f   ,  "%0.0
    input_text[26] = " XCH4=  %.16f   ,  "%0.0
    input_text[27] = " XN2O=  %.16f   ,  "%0.0
    input_text[28] = " XCO=   %.16f   ,  "%0.0
    input_text[29] = " XNO2=  %.16f   ,  "%0.0
    input_text[30] = " XSO2=  %.16f   ,  "%0.0
    input_text[31] = " XNH3=  %.16f   ,  "%0.0
    input_text[32] = " XNO=   %.16f   ,  "%0.0
    input_text[33] = " XHNO3= %.16f   ,  "%0.0
      
  
  if surface == "uniform":
    input_text[35] = " ISALB= 0  , "
    input_text[36] = " ALBCON= %.16f ,"%albedo
  elif surface == "snow":
    #input_text[35] = " ISALB= 1  , "
    input_text[35] = " ISALB= -1  , "
    os.system("cp %s/%s %s/albedo.dat"%(workdir,icefile,workdir))
  elif surface == "clearwater":
    input_text[35] = " ISALB= 2  , "
  elif surface == "lakewater":
    input_text[35] = " ISALB= 3  , "
  elif surface == "seawater":
    #input_text[35] = " ISALB= 4  , "
    input_text[35] = " ISALB= -1  , "
    os.system("cp %s/%s %s/albedo.dat"%(workdir,waterfile,workdir))
  elif surface == "sand":
    input_text[35] = " ISALB= 5  , "
  elif surface == "soil":
    input_text[35] = " ISALB = 10, "
    input_text[37] = " SC = 0.000,0.1,0.85,0.05  , "
  elif surface == "vegetation":
    input_text[35] = " ISALB= 6  , "
  elif surface == "seamix":
    input_text[35] = " ISALB= 10 , "
    input_text[37] = " SC= %.10f,%.10f,0.000,0.000   , "%(sic,1-sic)
  else:
    input_text[35] = " ISALB= -1  , "

  if flat:
      input_text[35] = " ISALB= 0  , "
      input_text[10] = " NF=       0,    "
      input_text[69] = " NOTHRM=   1 , "
      input_text[36] = " ALBCON= %.16f ,"%1.0 #mirror

  inputtext = '\n'.join(input_text)
  
  inp = open(workdir+"/INPUT","w")
  inp.write(inputtext)
  inp.close()

def write_input(workdir,cszenith,azimuth,latitude,surface,pCO2,p0,
                tsurf,altz,flux,icefile='seaice.dat',waterfile='seawater.dat',
                wmin=0.55,albedo=0.35,smooth=False,flat=True,sic=0.0,spec=False):
  template =("&INPUT                                                     \n"+ #0
             " IDATM=          0,                                        \n"+ #1
             " AMIX= -1.0000000000000000     ,                           \n"+ #2
             " ISAT=          -4,                                        \n"+ #3
             " WLINF= 10.0000000000000     ,                             \n"+ #4  
             " WLSUP= 4.7250000000000     ,                              \n"+ #5  
             " WLINC=  20  ,                                           \n"+ #6 #set to -0.01 for lower res
             " SZA=  0.0000000000000000     ,                            \n"+ #7
             " CSZA= -1.0000000000000000     ,                           \n"+ #8   --
             " SOLFAC=  1.0000000000000000     ,                         \n"+ #9   --
             " NF=          3,                                           \n"+ #10  -- Star type
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
             " ALBCON=  0.0000000000000000     ,                         \n"+ #36  --
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
  
  wc = 0.5*(19.0+wmin)
  input_text[4] = " WLINF= %2.9f ,   "%wc
  input_text[5] = " WLSUP= %2.9f ,   "%(0.5*(wc-wmin))
  
  if cszenith <= 0.0: #There is no incident sunlight, so we only need the thermal range >2.5 microns
    input_text[4] = " WLINF= 10.972857000  ,    "   # as opposed to >0.55 microns. This saves us time.
    input_text[5] = " WLSUP= 4.23857210000 ,    "
                #The spectral binning if we do this isn't perfectly aligned with the full range, so
                #we'll need to devise a binning scheme for putting them all in the right bins.
#STAR TYPE
  if spec:
      input_text[10] = " NF=      -1,    "
    
#UZEN
  rlat = latitude*np.pi/180.0
  hours = np.linspace(0,84.375,num=16)*np.pi/180.0
  # cos(p)cos(l) = cos(z)
  
  # if longitude = pi, and lat = 0, cos(z) = -1
  cuz = np.cos(rlat)*np.cos(hours)
  # if longitude = pi, and lat = 0, z = pi
  uz = np.arccos(cuz)*180.0/np.pi
  # if longitude = pi, and lat = 0, sin(z) = 0
  suz = np.sin(uz*np.pi/180.0)
  
#Azimuths
  # if longitude = pi, and lat = 0, cos(phi) is NaN.
  cosazs = -cuz*np.sin(rlat) / (suz*np.cos(rlat))
  # -cos(z)sin(p)/(sin(z)cos(p)) = -tan(p)/tan(z)
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
  if cszenith < 0:
      input_text[7] = " SZA= %.16f    ,     "%(np.arccos(cszenith)*180.0/np.pi)
      input_text[8] = " CSZA= -1     ,      "
  input_text[75] = " SAZA= %3.16f   ,   "%azimuth
  input_text[79] = " BTEMP= %3.3f   ,   "%tsurf
  input_text[9] = " SOLFAC= %.16f   ,   "%(flux/1367.0)
  input_text[15] = " ZPRES= %.16f   ,   "%altz
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
  
  if smooth:
    input_text[23] = " XN2=   %.16f   ,  "%1.0e6 
    input_text[24] = " XO2=   %.16f   ,  "%0.0
    input_text[25] = " XCO2=  %.16f   ,  "%0.0
    input_text[26] = " XCH4=  %.16f   ,  "%0.0
    input_text[27] = " XN2O=  %.16f   ,  "%0.0
    input_text[28] = " XCO=   %.16f   ,  "%0.0
    input_text[29] = " XNO2=  %.16f   ,  "%0.0
    input_text[30] = " XSO2=  %.16f   ,  "%0.0
    input_text[31] = " XNH3=  %.16f   ,  "%0.0
    input_text[32] = " XNO=   %.16f   ,  "%0.0
    input_text[33] = " XHNO3= %.16f   ,  "%0.0
      
  
  if surface == "uniform":
    input_text[35] = " ISALB= 0  , "
    input_text[36] = " ALBCON= %.16f ,"%albedo
  elif surface == "snow":
    #input_text[35] = " ISALB= 1  , "
    input_text[35] = " ISALB= -1  , "
    os.system("cp %s/%s %s/albedo.dat"%(workdir,icefile,workdir))
  elif surface == "clearwater":
    input_text[35] = " ISALB= 2  , "
  elif surface == "lakewater":
    input_text[35] = " ISALB= 3  , "
  elif surface == "seawater":
    #input_text[35] = " ISALB= 4  , "
    input_text[35] = " ISALB= -1  , "
    os.system("cp %s/%s %s/albedo.dat"%(workdir,waterfile,workdir))
  elif surface == "sand":
    input_text[35] = " ISALB= 5  , "
  elif surface == "vegetation":
    input_text[35] = " ISALB= 6  , "
  elif surface == "seamix":
    input_text[35] = " ISALB= 10 , "
    input_text[37] = " SC= %.10f,%.10f,0.000,0.000   , "%(sic,1-sic)
  else:
    input_text[35] = " ISALB= -1  , "

  if flat:
      input_text[35] = " ISALB= 0  , "
      input_text[10] = " NF=       0,    "
      input_text[69] = " NOTHRM=   1 , "
      input_text[36] = " ALBCON= %.16f ,"%1.0 #mirror

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
        layer.append(str(z[k]*1.0e-3)) #km
        layer.append(str(p[k]*0.01))   #hPa
        layer.append(str(t[k]))
        layer.append(str(q[k]))
        layer.append(str(o[k]))
        line = ' '.join(layer)
        dat+=line+'\n'
    f=open(workdir+'/'+name,"w")
    f.write(dat)
    f.close()    
    
def cloudcolumn_single(dql,cld,workdir,name='usrcld.dat',dqi=np.array([]),cirrus=False):
    dat = ''
    nlev = len(cld)
    for k in range(nlev-1,-1,-1):
        layer = []
        layer.append(str(dql[k]))
        layer.append("8")
        if cirrus:
          layer.append(str(dqi[k]))
        else:
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
        if ta[i]!=ta[i-1]:
            zdh = gascon*(ta[i-1]-ta[i])/grav*(np.log(lev[i]/lev[i-1])/np.log(ta[i-1]/ta[i]))
        else: #isothermal
            zdh = ta[i]*gascon/grav*(np.log(lev[i])-np.log(lev[i-1]))
        zzf[i] = zzh + 0.5*zdh
        zzh += zdh
        i-=1
    zdh = -ta[0]*gascon/grav*np.log(0.5)*0.5
    zzf[0] = zzh + zdh
    zzf*=1.0e-3
    
    return zzf  

def analyzecell_plasim_locked(data,views,lat,lon,workdir,grav=9.80665,sol_dec=0.0,
                     sol_lon=0.0,smooth=False,clouds=True,istep=-1,cirrus=False):
  #cszenith,azimuth,surface,pCO2,p0,tsurf,altz
  
  surf = "uniform"
  
  lsm = data.variables['lsm'][-1,lat,lon]
  ts = data.variables['ts'][:]
  istep = min(istep,ts.shape[0]-1)
  
  if lsm < 0.5: #sea
    #sic = np.mean(data.variables['sic'][:,lat,lon])
    sic = data.variables['sic'][istep,lat,lon]
    if sic == 1.0:
      surf = "snow"
    elif sic == 0.0:
      surf = "seawater"
    else:
      surf = "seamix"
  else:
    #if (np.mean(data.variables['snd'][:,lat,lon])+np.mean(data.variables['glac'][:,lat,lon]))>=0.7:
    if ((data.variables['snd'][istep,lat,lon]+data.variables['glac'][istep,lat,lon])>=0.005 and 
        data.variables['as'][istep,lat,lon]>0.2):
        surf = "snow"
    else:
        surf = "soil"
    sic = 0.0
    #We could put in an actual algorithm to get albedo given soil wetness, or combine sand, veg, etc
    
  #tsurf = np.mean(data.variables['ts'][:,lat,lon])
  tsurf = data.variables['ts'][istep,lat,lon]
  
  altz = max(data.variables["sg"][-1,lat,lon]*1.0e-3,0.0) / grav #kilometers
  
  #ta = np.mean(data.variables["ta"][:,:,lat,lon],axis=0)
  ta = data.variables["ta"][istep,:,lat,lon]
  
  lvs = data.variables['lev']
  #lns = data.variables['lon']
  #lts = data.variables['lat']
  
  zs = getalt_single(ta,lvs,grav=grav,gascon=287.0)
  
  #ps = np.mean(data.variables['ps'][:,lat,lon])
  ps = data.variables['ps'][istep,lat,lon]
  pa = ps * lvs
  
  gascon0 = 8.3144598
  mmair = 0.028964
  mmvap = 0.018016
  
  satvap = 610.78*10.0**(7.5*(ta-273.15)/(ta+0.15))
  
  #relhum = np.mean(data.variables['hur'][:,:,lat,lon],axis=0)*0.01
  relhum = data.variables['hur'][istep,:,lat,lon]*0.01
  
  pvap = relhum*satvap
  
  rvap = 461.495
  rdry = 287.058
  
  rhohum = ((pa*100.0-pvap)*mmair+pvap*mmvap)/(gascon0*ta)
  
  #hus = np.mean(data.variables['hus'][:,:,lat,lon],axis=0)
  hus = data.variables['hus'][istep,:,lat,lon]
  rhoh2o = hus*rhohum
  
  print "Writing atms.dat"
  #write atms.dat, which contains the atmosphere profile--height, pressure, temp, water vapor, and ozone
  for vw in views:
    writecolumn_single_plasim(zs,pa,ta,rhoh2o*1000.0,np.zeros(len(lvs)),workdir+"_%s"%vw)
  
  dsigma = np.zeros(len(lvs))
  dsigma[0] = lvs[0]
  for i in range(1,len(lvs)):
    dsigma[i] = lvs[i]-lvs[i-1]
    
  #clw = np.mean(data.variables['clw'][:,:,lat,lon],axis=0)
  clw = data.variables['clw'][istep,:,lat,lon]
  
  dql = np.zeros(len(lvs))
  cld = np.zeros(len(lvs))
  
  if not smooth:
    dql = dsigma*1000.0*ps*hus/grav #g/m^2 #Cloud liquid water path #was clw not hus
  
    #cld = np.mean(data.variables['cl'][:,:,lat,lon],axis=0)
    cld = data.variables['cl'][istep,:,lat,lon]
  
  if not clouds:
    cld[:] = 0.0
    
  print "Writing usrcld.dat"
  #Need to implement a cloud-free option
  #Write usrcld.dat, which has level data on cloud water content and coverage fraction
  
  if cirrus:
    dqi = -np.ones(dql.shape)
    for k in range(len(lvs)):
        if ta[k]<273.15:
            dqi[k] = dql[k]
            dql[k] = -1
  else:
    dqi = np.array([])
           
  
  for vw in views:
    cloudcolumn_single(dql,cld,workdir+"_%s"%vw,dqi=dqi,cirrus=cirrus)
  
  latitude = data.variables['lat'][lat]
  longitude = data.variables['lon'][lon]
  plon = data.variables['lon'][lon-1]
              
  csz = cosine_zenith(latitude,longitude,sol_dec,sol_lon)
  
  nlcsz = cosine_zenith(latitude,plon,sol_dec,sol_lon)
  
  if nlcsz>csz:
    direction = 'W' #sun is west of us, since it's brighter to the west
  else:
    direction = 'E' #sun is east of us, since it's brighter to the east
    
  ssz = np.sin(np.arccos(csz))
  clat = data.variables['lat'][lat] * np.pi/180.0
  
  cosaz = -csz*np.sin(clat) / (ssz*np.cos(clat))
  if np.isnan(cosaz):
    cosaz = 1.0
    if np.sin(clat)<0:
        cosaz = -1.0
  if np.isinf(cosaz):
    cosaz = 1.0
    if np.sin(clat)<0:
        cosaz = -1.0
  #print cosaz
  azm = np.arccos(np.minimum(np.maximum(cosaz,-1.0),1.0)) * 180.0/np.pi
  
  if direction=='W':
    azm = 360.0 - azm #Make it between 180 and 360
  
  return csz,azm,surf,sic,tsurf,altz

def analyzecell_plasim(data,lat,lon,workdir,grav=9.80665,smooth=False,clouds=True):
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
  
  relhum = np.mean(data.variables['hur'][:,:,lat,lon],axis=0)*0.01
  
  pvap = relhum*satvap
  
  rvap = 461.495
  rdry = 287.058
  
  rhohum = ((pa*100.0-pvap)*mmair+pvap*mmvap)/(gascon0*ta)
  
  hus = np.mean(data.variables['hus'][:,:,lat,lon],axis=0)
  rhoh2o = hus*rhohum
  
  print "Writing atms.dat"
  #write atms.dat, which contains the atmosphere profile--height, pressure, temp, water vapor, and ozone
  writecolumn_single_plasim(zs,pa,ta,rhoh2o*1000.0,np.zeros(len(lvs)),workdir)
  
  dsigma = np.zeros(len(lvs))
  dsigma[0] = lvs[0]
  for i in range(1,len(lvs)):
    dsigma[i] = lvs[i]-lvs[i-1]
    
  clw = np.mean(data.variables['clw'][:,:,lat,lon],axis=0)
  
  dql = np.zeros(len(lvs))
  cld = np.zeros(len(lvs))
  
  if not smooth:
    dql = dsigma*1000.0*ps*clw/grav #g/m^2 #Cloud liquid water path
  
    cld = np.mean(data.variables['cl'][:,:,lat,lon],axis=0)
  
  if not clouds:
    cld[:] = 0.0
    dql[:] = 0.0
    
  print "Writing usrcld.dat"
  #Need to implement a cloud-free option
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
    cosaz = 1.0
    if np.sin(clat)<0:
        cosaz = -1.0
  if np.isinf(cosaz):
    cosaz = 1.0
    if np.sin(clat)<0:
        cosaz = -1.0
  #print cosaz
  azm = np.arccos(np.minimum(np.maximum(cosaz,-1.0),1.0)) * 180.0/np.pi
  
  if direction=='W':
    azm = 360.0 - azm #Make it between 180 and 360
  
  return csz,azm,surf,sic,tsurf,altz
    
def cosine_zenith(latitude,longitude,sol_dec,sol_lon):
    #If latitude = 0 and longitude = 180:
  sinlat = np.sin(latitude*np.pi/180.)
    #sinlat = 0
  sindec = np.sin(sol_dec*np.pi/180.)
    #sindec = 0
  coslat = np.cos(latitude*np.pi/180.)
    #coslat = 1
  cosdec = np.cos(sol_dec*np.pi/180.)
    #cosdec = 1
  coslon = np.cos((longitude-sol_lon)*np.pi/180.0)
   #coslon = -1
  coszen = sinlat*sindec + coslat*cosdec*coslon
   #coszen = -1
  #coszen[abs(longitude-sol_lon)>90] = 0.0
  
  if type(coszen)==type(-0.5) or type(coszen)==np.float64:
      if coszen<-0.99:
          coszen = -0.99 #Just in case sbdart is tripping over cos(z)=-1
  else:
      coszen[coszen<-0.99] = -0.99
  
  return coszen

def analyzecell_lmdz(data,lat,lon,workdir,grav=9.80665,sol_dec=0.0,
                     sol_lon=0.0,smooth=False,clouds=True):
  #cszenith,azimuth,surface,pCO2,p0,tsurf,altz
  
  surf = "uniform"
  
  #lsm = data['lsm'][-1,lat,lon]
  lsm = 0.0
  if lsm < 0.5: #sea
    if 'pctsrf_sic' in data.keys():
      sic = data['pctsrf_sic'][lat,lon]
    else:
      sic = min(data['h2o_ice_surf'][lat,lon],1.0)
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
  
  rhohum = ((pa-pvap)*mmair+pvap*mmvap)/(gascon0*ta)
  
  hus = np.maximum(data['h2o_vap'][:,lat,lon],0.0) #kg/kg?
  rhoh2o = hus*rhohum
  if smooth:
      rhoh2o = np.zeros(hus.shape)
  
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
  
  dql = np.zeros(hus.shape)
  cld = np.zeros(hus.shape)
  dqi = np.zeros(hus.shape)
  rei = np.zeros(hus.shape)
  
  if not smooth:
    dql = np.maximum(dpress*clw/grav,0.0) #g/m^2 #Water mass in cell
    
    cld = data['CLF'][:,lat,lon]
    
    dqi = np.maximum(dpress*data["h2o_ice"][:,lat,lon]*1000/grav,0.0)
    
    rei = np.minimum(data["H2Oice_reff"][:,lat,lon]*1.0e6,128.0) #um
  
  if not clouds:
      cld[:] = 0.0
  
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
    cosaz = 1.0
    if np.sin(clat)<0:
        cosaz = -1.0
  if np.isinf(cosaz):
    cosaz = 1.0
    if np.sin(clat)<0:
        cosaz = -1.0
  #print cosaz
  azm = np.arccos(np.minimum(np.maximum(cosaz,-1.0),1.0)) * 180.0/np.pi
  
  if direction=='W':
    azm = 360.0 - azm #Make it between 180 and 360
  
  return csz,azm,surf,sic,tsurf,altz
    
def prep(job):
  if "type" not in job.parameters:
      print "Warning: need to specify which GCM produced this data!"
  else:
      if "ROLE" in job.parameters:
          if "outhopper" in job.parameters:
              dest = job.top+"/sbdart_locked/"+job.parameters["outhopper"]
          else:
              dest = job.top+"/sbdart_locked/output"
          os.system("mkdir "+dest)
          os.system("mkdir "+dest+"/running")
          os.system("mkdir "+dest+"/finished")
      if job.parameters["type"]=="plasim":
          os.system("touch "+dest+"/running/token"+job.pid+".crwl")
          _prep_plasim(job)
      elif job.parameters["type"]=="plasim_locked":
          os.system("touch "+dest+"/running/token"+job.pid+".crwl")
          _prep_plasim_locked(job)
      elif job.parameters["type"]=="lmdz":
          os.system("touch "+dest+"/running/token"+job.pid+".crwl")
          _prep_lmdz(job)
      elif job.parameters["type"]=="mit":
          print "MITGCM-SBDART interface not yet implemented."
      else:
          print "Model not recognized."
          
    
def _prep_lmdz(job):    
  workdir = job.top+"/sbdart_locked/job"+str(job.home)
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
  
  if "outhopper" in job.parameters:
    dest = job.top+"/sbdart_locked/"+job.parameters["outhopper"]
    os.system("mkdir "+job.top+"/sbdart_locked/"+job.parameters["outhopper"])
  else:
    dest = job.top+"/sbdart_locked/output"
    
  token_name = "token"+job.pid+"_%02d-%02d-%02d-%02d.crwl"%(lats[0],lats[1],lons[0],lons[1])
  
  os.system("mv "+dest+"/running/token"+job.pid+".crwl "+
            dest+"/running/"+token_name)
  
  
    
  data = np.load(job.top+"/hopper/"+job.parameters["gcm"]).item()  
  
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
    
  if "domainlatlon" in job.parameters:
      domain=job.parameters["domainlatlon"].split('/')
      lat1=domain[0]
      lat2=domain[1]
      lon1=domain[2]
      lon2=domain[3]
  else:
      lat1=0
      lat2=48
      lon1=0
      lon2=64
    
  star = False
    
  if "starspec" in job.parameters:
      star = job.parameters["starspec"]
    
  smooth=False  
  if "smooth" in job.parameters:
      if int(job.parameters["smooth"])==1:
          smooth=True
    
  flat=False
  if "flat" in job.parameters:
      if int(job.parameters["flat"])==1:
          flat=True
    
  clouds=True
  if "clouds" in job.parameters:
      if int(job.parameters["clouds"])==0:
          clouds=False
    
  wmin=0.55 #microns
  if "wmin" in job.parameters:
      wmin = float(job.parameters["wmin"])
    
  uniform=False     
  unialb = 0.35
  if "albedo" in job.parameters:
      uniform=True
      unialb = float(job.parameters["unialb"])
      
  notify = 'a'
  if "notify" in job.parameters:
      notify = job.parameters["notify"]
      
  nlats = len(data['latitude'][:])
  for jlat in range(lats[0],lats[1]):
    for jlon in range(lons[0],lons[1]):
      print "Lat %02d Lon %02d"%(jlat,jlon)
      print "mkdir "+workdir+"/sbdart-%02d_%02d"%(jlat,jlon)
      os.system("mkdir "+workdir+"/sbdart-%02d_%02d"%(jlat,jlon))
      os.system("cp -r "+job.top+"/sbdart_locked/"+source+"/* "+workdir+"/sbdart-%02d_%02d/"%(jlat,jlon))
      if star:
          #os.system("cp -r "+job.top+"/sbdart_locked/"+star+" "+workdir+"/sbdart-%02d_%02d/solar.dat"%(jlat,jlon))
          w,f = np.loadtxt(job.top+"/sbdart_locked/"+star,unpack=True)
          ffac = flux/np.trapz(f,x=w)
          f*=ffac
          data = np.array([w,f])
          data = data.T
         np.savetxt(workdir+"/sbdart-%02d_%02d/solar.dat"%(jlat,jlon),data,fmt=['%f','%f'])
  
  role = "dom"
  if "ROLE" in job.parameters:
      role = job.parameters["ROLE"] #"dom" or "sub" for dominant/submissive
  
  
  color=False
  if "color" in job.parameters:
      if job.parameters["color"]=="True":
          color=True
      else:
          color=False
  
  makemap=False
  if "map" in job.parameters:
      if job.parameters["map"]=="True":
          makemap=True
      else:
          makemap=False
  
  
  tag = ''
  if color:
      tag+="color "
  if makemap:
      tag+="map "  
  
  if role == "sub":
      jobscript =(bts.BATCHSCRIPT(job,notify)+
                  "module load gcc/4.9.1                                          \n"+
                  "module load python/2.7.9                                       \n"+
                  "for jl in {%02d..%02d};                                  \n"%(lats[0],lats[1]-1)+
                  "do \n"+
                  "     for il in {%02d..%02d};                        \n"%(lons[0],lons[1]-1)+
                  "     do \n"+
                  "          ILAT=`printf '%02d' $(( 10#$jl ))`           \n"+
                  "          ILON=`printf '%02d' $(( 10#$il ))`           \n"+
                  "          echo $ILAT $ILON              \n"+
                  "          TAG=${ILAT}_${ILON}                    \n"+
                  "          cd "+workdir+"/sbdart-$TAG                          \n"+
                  "          ./sbdart > "+workdir+"/sbout.$TAG                \n"+
                  "          cd "+workdir+"                                  \n"+
                  "     done                                    \n"+
                  "done \n"+
                  "cp "+dest+"/running/"+token_name+" "+dest+"/finished/ \n"+
                  './release.sh "'+dest+'"                                \n')
  else:
      jobscript =(bts.BATCHSCRIPT(job,notify)+
                  "module load gcc/4.9.1                                          \n"+
                  "module load python/2.7.9                                       \n"+
                  "for jl in {%02d..%02d};                                  \n"%(lats[0],lats[1]-1)+
                  "do \n"+
                  "     for il in {%02d..%02d};                        \n"%(lons[0],lons[1]-1)+
                  "     do \n"+
                  "          ILAT=`printf '%02d' $(( 10#$jl ))`           \n"+
                  "          ILON=`printf '%02d' $(( 10#$il ))`           \n"+
                  "          echo $ILAT $ILON              \n"+
                  "          TAG=${ILAT}_${ILON}                    \n"+
                  "          cd "+workdir+"/sbdart-$TAG                          \n"+
                  "          ./sbdart > "+workdir+"/sbout.$TAG                \n"+
                  "          cd "+workdir+"                                  \n"+
                  "     done                                    \n"+
                  "done \n"+
                  "cp "+dest+"/running/"+token_name+" "+dest+"/finished/ \n"+
                  './release.sh "'+dest+'"                                \n'+
                  "python -B checkprogress.py "+dest+" "+lat1+" "+lat2+" "+lon1+" "+lon2+" 32 "+job.top+
                  " "+job.parameters["type"]+" "+job.parameters["gcm"]+" "+tag+"         \n")
      
  
  rs = open(workdir+"/runsbdart","w")
  rs.write(jobscript)
  rs.close()
  
  os.system("cp -r "+job.top+"/sbdart_locked/release.sh "+workdir+"/")
  os.system("cp -r "+job.top+"/sbdart_locked/checkprogress.py "+workdir+"/")
  os.system("cp "+job.top+"/crawldefs.py "+workdir+"/")
  os.system("cp "+job.top+"/identity.py "+workdir+"/")
  
  for jlon in range(lons[0],lons[1]):
    for jlat in range(lats[0],lats[1]):
      csz,azm,surf,sic,tsurf,altz = analyzecell_lmdz(data,jlat,jlon,
                                                     workdir+"/sbdart-%02d_%02d"%(jlat,jlon),
                                                     grav=grav,smooth=smooth,clouds=clouds)
      if uniform:
          surf='uniform'
      latitude = data['latitude'][jlat]
      write_input(workdir+"/sbdart-%02d_%02d"%(jlat,jlon),csz,azm,latitude,surf,pCO2,p0,tsurf,altz,
                  flux,wmin=wmin,albedo=unialb,flat=flat,sic=sic,spec=star,smooth=smooth)
      print "Prepped lat %02d lon %02d"%(jlat,jlon)

    
def _prep_plasim_locked(job): #data,lats,lons,pCO2,p0,flux,grav=9.80665
  os.system("rm "+job.top+"/sbdart_locked/job"+str(job.home)+"/*.e*")
  os.system("rm "+job.top+"/sbdart_locked/job"+str(job.home)+"/*.o*")
  workdir = ide.SCRATCH+"/sbdart_locked_job"+str(job.home)
  os.system("mkdir "+workdir)
  os.system("rm "+workdir+"/*.e* "+workdir+"/*.o*")
  if "source" in job.parameters:
    source = job.parameters["source"]
  else:
    source = "clean"
    
  ntimes = '{0,1,2,3}'
  lviews = '{Z,E,W,N,S}'
  
  iout=5
  
  if "iout" in job.parameters:
      iout = int(job.parameters["iout"])
  
  if "angles" in job.parameters:
      ntimes = job.parameters["angles"]
    
  if "views" in job.parameters:
      lviews = job.parameters["views"]
      
  vws = lviews[1:-1].split(',')
  if vws[-1]=='':
      vws = vws[:-1]
  itimes = ntimes[1:-1].split(',')
  if itimes[-1]=='':
      itimes = itimes[:-1]
  itimes = np.array(itimes).astype(int)
      
  
  viewdict = {'Z':0,
              'E':1,
              'W':2,
              'N':3,
              'S':4}    
      
  print "Setting up for times ",itimes,"and views",vws
    
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
    
  tempdest = workdir+"/output/"
  os.system("mkdir "+tempdest)
    
  if "outhopper" in job.parameters:
    finaldest = job.top+"/sbdart_locked/"+job.parameters["outhopper"]
    os.system("mkdir "+job.top+"/sbdart_locked/"+job.parameters["outhopper"])
  else:
    finaldest = job.top+"/sbdart_locked/output"
    
  token_name = "token"+job.pid+"_%02d-%02d-%02d-%02d.crwl"%(lats[0],lats[1],lons[0],lons[1])
  
  os.system("mv "+finaldest+"/running/token"+job.pid+".crwl "+
            finaldest+"/running/"+token_name)
  
  os.system("cp "+job.top+"/hopper/"+job.parameters["gcm"]+" "+workdir+"/"+job.parameters["gcm"])
  data = nc.Dataset(workdir+"/"+job.parameters["gcm"],"r")  
  
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
    
  if "domainlatlon" in job.parameters:
      domain=job.parameters["domainlatlon"].split('/')
      lat1=domain[0]
      lat2=domain[1]
      lon1=domain[2]
      lon2=domain[3]
  else:
      lat1=0
      lat2=48
      lon1=0
      lon2=64
    
  if job.ncores>1:
      mode="batch"
      if (lats[1]-lats[0])%job.ncores==0:
          latpairs = []
          dlat = (lats[1]-lats[0])/job.ncores
          for ll in range(lats[0],lats[1],dlat):
              latpairs.append([ll,ll+dlat])
          lonpairs = [lons,]
      elif (lats[1]-lats[0]==1) and (lons[1]-lons[0])%job.ncores==0:
          latpairs = [lats,]
          lonpairs = []
          dlon = (lons[1]-lons[0])/job.ncores
          for ll in range(lons[0],lons[1],dlon):
              lonpairs.append([ll,ll+dlon])
      elif ((lats[1]-lats[0])*(lons[1]-lons[0]))%job.ncores == 0:
          if (lats[1]-lats[0])<job.ncores:
            nlonchunks = job.ncores/(lats[1]-lats[0])
            dlons = (lons[1]-lons[0])/nlonchunks
            latpairs = []
            lonpairs = []
            for ll in range(lats[0],lats[1]):
                latpairs.append([ll,ll+1])
            for ll in range(lons[0],lons[1],dlons):
                lonpairs.append([ll,ll+dlons])
          else:
            nlatchunks = job.ncores/(lons[1]-lons[0])
            dlats = (lats[1]-lats[0])/nlatchunks
            latpairs = []
            lonpairs = []
            for ll in range(lats[0],lats[1],dlats):
                latpairs.append([ll,ll+dlats])
            for ll in range(lons[0],lons[1]):
                lonpairs.append([ll,ll+1])
      else:
          mode="single"
      if mode!="single":
          for n in range(len(latpairs)):
              for l in range(len(latpairs[n])):
                  latpairs[n][l] = min(31,latpairs[n][l])
          for n in range(len(lonpairs)):
              for l in range(len(lonpairs[n])):
                  lonpairs[n][l] = min(63,lonpairs[n][l])
  else:
      mode="single"
    
  star=False
  if "starspec" in job.parameters:
      star = job.parameters["starspec"]
  
  smooth=False  
  if "smooth" in job.parameters:
      if int(job.parameters["smooth"])==1:
          smooth=True
  
  flat=False
  if "flat" in job.parameters:
      if int(job.parameters["flat"])==1:
          flat=True
    
  clouds=True
  if "clouds" in job.parameters:
      if int(job.parameters["clouds"])==0:
          clouds=False
    
  notify = 'a'
  if "notify" in job.parameters:
      notify = job.parameters["notify"]
      
  istep=12
  if "istep" in job.parameters:
      istep = int(job.parameters["istep"])
      
  uniform=False     
  unialb = 0.35
  if "albedo" in job.parameters:
      uniform=True
      unialb = float(job.parameters["albedo"])
      

  wmin=0.55 #microns
  if "wmin" in job.parameters:
      wmin = float(job.parameters["wmin"])
  
  wmax = 19.0 #microns
  if "wmax" in job.parameters:
      wmax = float(job.parameters["wmax"])
  
  
  zout = (0.,100.)
  if "zout" in job.parameters:
      zout = (0.,float(job.parameters["zout"]))
          
      
  nlats = len(data.variables['lat'][:])
  
  os.system("mkdir "+workdir+"/source/")
  #os.system("tar cvzf "+job.top+"/sbdart_locked/"+source+"/source.tar.gz "+
            #job.top+"/sbdart_locked/"+source+"/* ")
  os.system("rsync "+job.top+"/sbdart_locked/"+source+"/source.tar.gz "+workdir+"/source/")
  #os.system("rm "+job.top+"/sbdart_locked/"+source+"/source.tar.gz ")
  os.system("tar xvzf "+workdir+"/source/source.tar.gz -C "+workdir+"/source/")
  os.system("rm "+workdir+"/source/source.tar.gz")
  if star:
      #os.system("cp -r "+job.top+"/sbdart_locked/"+star+" "+workdir+"/source/solar.dat")
      w,f = np.loadtxt(job.top+"/sbdart_locked/"+star,unpack=True)
      ffac = flux/np.trapz(f,x=w)
      f*=ffac
      data = np.array([w,f])
      data = data.T
      np.savetxt(workdir+"/source/solar.dat",data,fmt=['%f','%f'])
  
  for jlat in range(lats[0],lats[1]):
    for jlon in range(lons[0],lons[1]):
      for nang in itimes:
        for vw in vws:
            #print "Lat %02d Lon %02d Angle %1d View %s"%(jlat,jlon,nang,vw)
            #print "mkdir "+workdir+"/sbdart-%02d_%02d_%1d_%s"%(jlat,jlon,nang,vw)
            os.system("mkdir "+workdir+"/sbdart-%02d_%02d_%1d_%s"%(jlat,jlon,nang,vw))
            os.system("cp -r "+workdir+"/source/* "+workdir+"/sbdart-%02d_%02d_%1d_%s/"%(jlat,jlon,nang,vw))
        
  
  role = "dom"
  if "ROLE" in job.parameters:
      role = job.parameters["ROLE"] #"dom" or "sub" for dominant/submissive
  
   
  color=False
  if "color" in job.parameters:
      if job.parameters["color"]=="True":
          color=True
      else:
          color=False
  
  makemap=False
  if "map" in job.parameters:
      if job.parameters["map"]=="True":
          makemap=True
      else:
          makemap=False

  cirrus=False
  if "cirrus" in job.parameters:
      if job.parameters["cirrus"]=="1":
          cirrus=True
      else:
          cirrus=False  
  
  icefile='seaice.dat'
  waterfile='seawater.dat'
  if "icespec" in job.parameters:
      icefile=job.parameters["icespec"]
  if "waterspec" in job.parameters:
      waterfile=job.parameters["waterspec"]
      
  tag = ''
  if color:
      tag+="color "
  if makemap:
      tag+="map "
  
  if mode!="single":
      role="dom"
  
  if role=="sub":
      if mode=="single":
        jobscript =(bts.BATCHSCRIPT(job,notify)+
                  "module load gcc/4.9.1                                          \n"+
                  "module load python/2.7.9                                       \n"+
                  "for al in "+ntimes+";                         \n"+
                  "do \n"+
                  "     for vw in "+lviews+";                \n"+
                  "     do \n"+     
                  "          for jl in {%02d..%02d};                       \n"%(lats[0],lats[1]-1)+
                  "          do \n"+
                  "               for il in {%02d..%02d};                  \n"%(lons[0],lons[1]-1)+
                  "               do \n"+
                  "                    ILAT=`printf '%02d' $(( 10#$jl ))`           \n"+
                  "                    ILON=`printf '%02d' $(( 10#$il ))`           \n"+
                  "                    IANG=`printf '%1d' $al`           \n"+
                  "                    echo $ILAT $ILON $IANG $vw              \n"+
                  "                    TAG=${ILAT}_${ILON}_${IANG}_$vw                    \n"+
                  "                    cd "+workdir+"/sbdart-$TAG                          \n"+
                  "                    ./sbdart > "+workdir+"/output/sbout.$TAG            \n"+
                  "                    cd "+workdir+"                                  \n"+
                  "               done                       \n"+
                  "          done                         \n"+
                  "          mv "+workdir+"/output/sbout.* "+finaldest+"/                   \n"+
                  "     done                                    \n"+
                  "done \n"+
                  "cp "+finaldest+"/running/"+token_name+" "+finaldest+"/finished/ \n"+
                  './release.sh "'+finaldest+'"                                \n'+
                  "rm -rf "+workdir+"/*/                                  \n")
      else:
          runscript = ("#!/bin/bash \n\n"+
                      "cd "+workdir+" \n"+
                  "for al in "+ntimes+";                         \n"+
                  "do \n"+
                  "     for vw in "+lviews+";                \n"+
                  "     do \n"+     
                  "          for ((jl=$1; jl<=$2; jl++));                       \n"+
                  "          do \n"+
                  "               for ((il=$3; il<=$4; il++));                  \n"+
                  "               do \n"+
                  "                    ILAT=`printf '%02d' $(( 10#$jl ))`           \n"+
                  "                    ILON=`printf '%02d' $(( 10#$il ))`           \n"+
                  "                    IANG=`printf '%1d' $al`           \n"+
                  "                    echo $ILAT $ILON $IANG $vw              \n"+
                  "                    TAG=${ILAT}_${ILON}_${IANG}_$vw                    \n"+
                  "                    cd "+workdir+"/sbdart-$TAG                          \n"+
                  "                    ./sbdart > "+workdir+"/output/sbout.$TAG            \n"+
                  "                    cd "+workdir+"                                  \n"+
                  "               done                       \n"+
                  "          done                         \n"+
                  "          mv "+workdir+"/output/sbout.* "+finaldest+"/                   \n"+
                  "     done                                    \n"+
                  "done \n")
          with open(workdir+"/sbdart_batch.sh","w") as rs:
              rs.write(runscript)
          os.system("chmod a+x "+workdir+"/sbdart_batch.sh")
          jobscript =(bts.BATCHSCRIPT(job,notify)+
                      "module load gcc/4.9.1 \n"+
                      "module load python/2.7.9 \n"+
                      "cd "+workdir+" \n"+
                      "#ls -d * \n"+
                      "cat sbdart_batch.sh  \n\n")
          for ltp in latpairs:
              for lnp in lonpairs:
                  jobscript += "./sbdart_batch.sh %2d %2d %2d %2d &    \n\n"%(ltp[0],ltp[1],lnp[0],lnp[1])
          jobscript += ("wait     \n"+
                        "cp "+finaldest+"/running/"+token_name+" "+finaldest+"/finished/ \n"+
                  './release.sh "'+finaldest+'"                                \n'+
                  "rm -rf "+workdir+"/*/                                  \n")
          
      
  else:
      if mode=="single":
        jobscript =(bts.BATCHSCRIPT(job,notify)+
                  "module load gcc/4.9.1                                          \n"+
                  "module load python/2.7.9                                       \n"+
                  "for al in "+ntimes+";                          \n"+
                  "do \n"+
                  "     for vw in "+lviews+";                \n"+
                  "     do \n"+     
                  "          for jl in {%02d..%02d};                       \n"%(lats[0],lats[1]-1)+
                  "          do \n"+
                  "               for il in {%02d..%02d};                  \n"%(lons[0],lons[1]-1)+
                  "               do \n"+
                  "                    ILAT=`printf '%02d' $(( 10#$jl ))`           \n"+
                  "                    ILON=`printf '%02d' $(( 10#$il ))`           \n"+
                  "                    IANG=`printf '%1d' $al`           \n"+
                  "                    echo $ILAT $ILON $IANG $vw              \n"+
                  "                    TAG=${ILAT}_${ILON}_${IANG}_$vw                    \n"+
                  "                    cd "+workdir+"/sbdart-$TAG                          \n"+
                  "                    ./sbdart > "+workdir+"/output/sbout.$TAG            \n"+
                  "                    cd "+workdir+"                                  \n"+
                  "               done                       \n"+
                  "          done                         \n"+
                  "          mv "+workdir+"/output/sbout.* "+finaldest+"/                   \n"+
                  "     done                                    \n"+
                  "done \n"+
                  "cp "+finaldest+"/running/"+token_name+" "+finaldest+"/finished/ \n"+
                  './release.sh "'+finaldest+'"                                \n'+
                  "rm -rf "+workdir+"/*/                                  \n"+
                  "#python -B checkprogress_locked.py "+finaldest+" "+lat1+" "+lat2+" "+lon1+" "+lon2+" 0 "+job.top+
                  " "+job.parameters["type"]+" "+job.parameters["gcm"]+" "+tag+"         \n")
      else:
          runscript = ("#!/bin/bash \n\n"+
                      "cd "+workdir+" \n"+
                  "for al in "+ntimes+";                         \n"+
                  "do \n"+
                  "     for vw in "+lviews+";                \n"+
                  "     do \n"+     
                  "          for ((jl=$1; jl<=$2; jl++));                       \n"+
                  "          do \n"+
                  "               for ((il=$3; il<=$4; il++));                  \n"+
                  "               do \n"+
                  "                    ILAT=`printf '%02d' $(( 10#$jl ))`           \n"+
                  "                    ILON=`printf '%02d' $(( 10#$il ))`           \n"+
                  "                    IANG=`printf '%1d' $al`           \n"+
                  "                    echo $ILAT $ILON $IANG $vw              \n"+
                  "                    TAG=${ILAT}_${ILON}_${IANG}_$vw                    \n"+
                  "                    cd "+workdir+"/sbdart-$TAG                          \n"+
                  "                    ./sbdart > "+workdir+"/output/sbout.$TAG            \n"+
                  "                    cd "+workdir+"                                  \n"+
                  "               done                       \n"+
                  "          done                         \n"+
                  "          mv "+workdir+"/output/sbout.* "+finaldest+"/                   \n"+
                  "     done                                    \n"+
                  "done \n")
          with open(workdir+"/sbdart_batch.sh","w") as rs:
              rs.write(runscript)
          os.system("chmod a+x "+workdir+"/sbdart_batch.sh")
          jobscript =(bts.BATCHSCRIPT(job,notify)+
                      "module load gcc/4.9.1 \n"+
                      "module load python/2.7.9 \n"+
                      "cd "+workdir+" \n"+
                      "#ls -d * \n"+
                      "cat sbdart_batch.sh  \n\n")
          for ltp in latpairs:
              for lnp in lonpairs:
                  jobscript += "./sbdart_batch.sh %2d %2d %2d %2d &    \n\n"%(ltp[0],ltp[1],lnp[0],lnp[1])
          jobscript += ("wait     \n"+
                  "cp "+finaldest+"/running/"+token_name+" "+finaldest+"/finished/ \n"+
                  './release.sh "'+finaldest+'"                                \n'+
                  "rm -rf "+workdir+"/*/                                  \n"+
                  "#python -B checkprogress_locked.py "+finaldest+" "+lat1+" "+lat2+" "+lon1+" "+lon2+" 0 "+job.top+
                  " "+job.parameters["type"]+" "+job.parameters["gcm"]+" "+tag+"         \n"+
                  "cd $PBS_O_WORKDIR   \n"+
                  "python -B release.py   \n")
  
  
  #print jobscript
  #print workdir
  rs = open(workdir+"/runsbdart","w")
  rs.write(jobscript)
  rs.close()
  os.system("chmod a+x "+workdir+"/runsbdart")
  os.system("cp -r "+job.top+"/sbdart_locked/release.sh "+workdir+"/")
  os.system("cp -r "+job.top+"/sbdart_locked/checkprogress_locked.py "+workdir+"/")
  os.system("cp "+job.top+"/crawldefs.py "+workdir+"/")
  os.system("cp "+job.top+"/identity.py "+workdir+"/")
  
  #views = [28.125,118.125,196.875,275.625]
  views = [0.,90.,180.,270.]
  
  for jlon in range(lons[0],lons[1]):
    for jlat in range(lats[0],lats[1]):
      for nang in itimes:
        csz,azm,surf,sic,tsurf,altz = analyzecell_plasim_locked(data,vws,jlat,jlon,
                                                         workdir+"/sbdart-%02d_%02d_%1d"%(jlat,jlon,nang),
                                                         sol_lon = views[nang],
                                                         grav=grav,smooth=smooth,clouds=clouds,
                                                         istep=istep,cirrus=cirrus)
        if uniform:
            surf='uniform'
        latitude = data.variables['lat'][jlat]
        longitude = data.variables['lon'][jlon]
        #vws = ['Z','N','E','S','W']
        for vv in range(0,len(vws)):
            nv = viewdict[vws[vv]]
            write_input_locked(workdir+"/sbdart-%02d_%02d_%1d_%s"%(jlat,jlon,nang,vws[vv]),
                              nv,nang,csz,azm,latitude,longitude,surf,pCO2,p0,tsurf,altz,flux,
                              wmin=wmin,wmax=wmax,albedo=unialb,flat=flat,sic=sic,spec=star,
                              smooth=smooth,iout=iout,zout=zout,waterfile=waterfile,
                              icefile=icefile)
        #print "Prepped lat %02d lon %02d Angle %1d View %s"%(jlat,jlon,nang,vws[vv])
      
def _prep_plasim(job): #data,lats,lons,pCO2,p0,flux,grav=9.80665
  workdir = job.top+"/sbdart_locked/job"+str(job.home)
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
    
  if "outhopper" in job.parameters:
    dest = job.top+"/sbdart_locked/"+job.parameters["outhopper"]
    os.system("mkdir "+job.top+"/sbdart_locked/"+job.parameters["outhopper"])
  else:
    dest = job.top+"/sbdart_locked/output"
    
  token_name = "token"+job.pid+"_%02d-%02d-%02d-%02d.crwl"%(lats[0],lats[1],lons[0],lons[1])
  
  os.system("mv "+dest+"/running/token"+job.pid+".crwl "+
            dest+"/running/"+token_name)
  
  data = nc.Dataset(job.top+"/hopper/"+job.parameters["gcm"],"r")  
  
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
    
  if "domainlatlon" in job.parameters:
      domain=job.parameters["domainlatlon"].split('/')
      lat1=domain[0]
      lat2=domain[1]
      lon1=domain[2]
      lon2=domain[3]
  else:
      lat1=0
      lat2=48
      lon1=0
      lon2=64
    
    
  star=False
  if "starspec" in job.parameters:
      star = job.parameters["starspec"]
  
  smooth=False  
  if "smooth" in job.parameters:
      if int(job.parameters["smooth"])==1:
          smooth=True
  
  flat=False
  if "flat" in job.parameters:
      if int(job.parameters["flat"])==1:
          flat=True
    
  clouds=True
  if "clouds" in job.parameters:
      if int(job.parameters["clouds"])==0:
          clouds=False
    
  notify = 'a'
  if "notify" in job.parameters:
      notify = job.parameters["notify"]
      
  uniform=False     
  unialb = 0.35
  if "albedo" in job.parameters:
      uniform=True
      unialb = float(job.parameters["unialb"])
      

  wmin=0.55 #microns
  if "wmin" in job.parameters:
      wmin = float(job.parameters["wmin"])
          
      
  nlats = len(data.variables['lat'][:])
  for jlat in range(lats[0],lats[1]):
    for jlon in range(lons[0],lons[1]):
      print "Lat %02d Lon %02d"%(jlat,jlon)
      print "mkdir "+workdir+"/sbdart-%02d_%02d"%(jlat,jlon)
      os.system("mkdir "+workdir+"/sbdart-%02d_%02d"%(jlat,jlon))
      os.system("cp -r "+job.top+"/sbdart_locked/"+source+"/* "+workdir+"/sbdart-%02d_%02d/"%(jlat,jlon))
      if star:
          #os.system("cp -r "+job.top+"/sbdart_locked/"+star+" "+workdir+"/sbdart-%02d_%02d/solar.dat"%(jlat,jlon))
          w,f = np.loadtxt(ob.top+"/sbdart_locked/"+star,unpack=True)
          ffac = flux/np.trapz(f,x=w)
          f*=ffac
          data = np.array([w,f])
          data = data.T
          np.savetxt(workdir+"/sbdart-%02d_%02d/solar.dat"%(jlat,jlon),data,fmt=['%f','%f'])
  
  
  role = "dom"
  if "ROLE" in job.parameters:
      role = job.parameters["ROLE"] #"dom" or "sub" for dominant/submissive
  
   
  color=False
  if "color" in job.parameters:
      if job.parameters["color"]=="True":
          color=True
      else:
          color=False
  
  makemap=False
  if "map" in job.parameters:
      if job.parameters["map"]=="True":
          makemap=True
      else:
          makemap=False
  
  
  tag = ''
  if color:
      tag+="color "
  if makemap:
      tag+="map "
  
  if role=="sub":
      jobscript =(bts.BATCHSCRIPT(job,notify)+
                  "module load gcc/4.9.1                                          \n"+
                  "module load python/2.7.9                                       \n"+
                  "for jl in {%02d..%02d};                                  \n"%(lats[0],lats[1]-1)+
                  "do \n"+
                  "     for il in {%02d..%02d};                        \n"%(lons[0],lons[1]-1)+
                  "     do \n"+
                  "          ILAT=`printf '%02d' $(( 10#$jl ))`           \n"+
                  "          ILON=`printf '%02d' $(( 10#$il ))`           \n"+
                  "          echo $ILAT $ILON              \n"+
                  "          TAG=${ILAT}_${ILON}                    \n"+
                  "          cd "+workdir+"/sbdart-$TAG                          \n"+
                  "          ./sbdart > "+workdir+"/sbout.$TAG                \n"+
                  "          cd "+workdir+"                                  \n"+
                  "     done                                    \n"+
                  "done \n"+
                  "cp "+dest+"/running/"+token_name+" "+dest+"/finished/ \n"+
                  './release.sh "'+dest+'"                                \n')
      
  else:
      jobscript =(bts.BATCHSCRIPT(job,notify)+
                  "module load gcc/4.9.1                                          \n"+
                  "module load python/2.7.9                                       \n"+
                  "for jl in {%02d..%02d};                                  \n"%(lats[0],lats[1]-1)+
                  "do \n"+
                  "     for il in {%02d..%02d};                        \n"%(lons[0],lons[1]-1)+
                  "     do \n"+
                  "          ILAT=`printf '%02d' $(( 10#$jl ))`           \n"+
                  "          ILON=`printf '%02d' $(( 10#$il ))`           \n"+
                  "          echo $ILAT $ILON              \n"+
                  "          TAG=${ILAT}_${ILON}                    \n"+
                  "          cd "+workdir+"/sbdart-$TAG                          \n"+
                  "          ./sbdart > "+workdir+"/sbout.$TAG                \n"+
                  "          cd "+workdir+"                                  \n"+
                  "     done                                    \n"+
                  "done \n"+
                  "cp "+dest+"/running/"+token_name+" "+dest+"/finished/ \n"+
                  './release.sh "'+dest+'"                                \n'+
                  "python -B checkprogress.py "+dest+" "+lat1+" "+lat2+" "+lon1+" "+lon2+" 0 "+job.top+
                  " "+job.parameters["type"]+" "+job.parameters["gcm"]+" "+tag+"         \n")
  
  rs = open(workdir+"/runsbdart","w")
  rs.write(jobscript)
  rs.close()
  
  os.system("cp -r "+job.top+"/sbdart_locked/release.sh "+workdir+"/")
  os.system("cp -r "+job.top+"/sbdart_locked/checkprogress.py "+workdir+"/")
  os.system("cp "+job.top+"/crawldefs.py "+workdir+"/")
  os.system("cp "+job.top+"/identity.py "+workdir+"/")
  
  for jlon in range(lons[0],lons[1]):
    for jlat in range(lats[0],lats[1]):
      csz,azm,surf,sic,tsurf,altz = analyzecell_plasim(data,jlat,jlon,
                                                       workdir+"/sbdart-%02d_%02d"%(jlat,jlon),
                                                       grav=grav,smooth=smooth,clouds=clouds)
      if uniform:
          surf='uniform'
      latitude = data.variables['lat'][jlat]
      write_input(workdir+"/sbdart-%02d_%02d"%(jlat,jlon),csz,azm,latitude,surf,pCO2,p0,tsurf,altz,
                  flux,wmin=wmin,albedo=unialb,flat=flat,sic=sic,spec=star,smooth=smooth)
      print "Prepped lat %02d lon %02d"%(jlat,jlon)

def submit(job):
  workdir = job.top+"/sbdart_locked/job"+str(job.home)
  
  os.system("cd "+workdir+" && "+bts.SUB+" runsbdart && cd "+job.top)
  
def run(job):
  workdir = ide.SCRATCH+"/sbdart_locked_job"+str(job.home)
  os.system("cd "+workdir+" && bash runsbdart")
   
def _prep(job): #data,lats,lons,pCO2,p0,flux,grav=9.80665
  print job
  
  with open("../.home","r") as homef:
      top = homef.read().split('\n')[0]
  
  #os.system("rm "+top+"/sbdart_locked/*.pyc")
  os.system("cp "+top+"/identity.py "+top+"/sbdart_locked/")
  os.system("cp "+top+"/crawldefs.py "+top+"/sbdart_locked/")
  os.system("cp "+top+"/torque.py "+top+"/sbdart_locked/")
  os.system("cp "+top+"/slurm.py "+top+"/sbdart_locked/")
  os.system("cp "+top+"/batch_system.py "+top+"/sbdart_locked/")
  
  jobname = job[0]
  jobncores = int(job[1])
  jobqueue = job[2]
  pco2 = float(job[3])
  flux = float(job[4])
  p0 = float(job[5])
  grav = float(job[6])
  ntimes = '{'+','.join(job[7].split('^'))+'}'
  lviews = '{'+','.join(job[8].split('^'))+"}"
  #ntimes = '{0,1,2,3}'
  #lviews = '{Z,E,W,N,S}'
  if job[9]=="True":
      highcadence=True
  else:
      highcadence=False
  try:
      star = job[10]
  except:
      star = False
  
  workdir = ide.SCRATCH+"/sbdart_locked_"+jobname
  os.system("mkdir "+workdir)

  source = "clean"
    
  
  iout=5
  
  vws = lviews[1:-1].split(',')
  if vws[-1]=='':
      vws = vws[:-1]
  itimes = ntimes[1:-1].split(',')
  if itimes[-1]=='':
      itimes = itimes[:-1]
  itimes = np.array(itimes).astype(int)     
  
  viewdict = {'Z':0,
              'E':1,
              'W':2,
              'N':3,
              'S':4}    
      
  print "Setting up for times ",itimes,"and views",vws
  
  outtag = "snapshot"
  if highcadence:
      outtag="highcadence"
  
    
  os.system("cp "+top+"/plasim/output/"+jobname+"_%s.nc "%outtag
            +workdir+"/"+jobname+"_%s.nc"%outtag)
  data = nc.Dataset(workdir+"/"+jobname+"_%s.nc"%outtag,"r")  
  
  nlats = len(data.variables['lat'][:])
  nlons = len(data.variables['lon'][:])
  
  lats = [0,nlats]
  lons = [0,nlons]
    
    
  tempdest = workdir+"/output/"
  os.system("mkdir "+tempdest)
  
  
  finaldest = top+"/sbdart_locked/"+jobname
  os.system("mkdir "+finaldest)
    
    
  token_name = "token_%02d-%02d-%02d-%02d.crwl"%(lats[0],lats[1],lons[0],lons[1])
  
  os.system("mkdir "+finaldest+"/running")
  os.system("touch "+finaldest+"/running/"+token_name)
     
  lat1=0
  lat2=nlats
  lon1=0
  lon2=nlons
  
  if jobncores>1:
      mode="batch"
      if (lats[1]-lats[0])%jobncores==0:
          latpairs = []
          dlat = (lats[1]-lats[0])/jobncores
          for ll in range(lats[0],lats[1],dlat):
              latpairs.append([ll,ll+dlat])
          lonpairs = [lons,]
      elif (lats[1]-lats[0]==1) and (lons[1]-lons[0])%jobncores==0:
          latpairs = [lats,]
          lonpairs = []
          dlon = (lons[1]-lons[0])/jobncores
          for ll in range(lons[0],lons[1],dlon):
              lonpairs.append([ll,ll+dlon])
      elif ((lats[1]-lats[0])*(lons[1]-lons[0]))%jobncores == 0:
          if (lats[1]-lats[0])<jobncores:
            nlonchunks = jobncores/(lats[1]-lats[0])
            dlons = (lons[1]-lons[0])/nlonchunks
            latpairs = []
            lonpairs = []
            for ll in range(lats[0],lats[1]):
                latpairs.append([ll,ll+1])
            for ll in range(lons[0],lons[1],dlons):
                lonpairs.append([ll,ll+dlons])
          else:
            nlatchunks = jobncores/(lons[1]-lons[0])
            dlats = (lats[1]-lats[0])/nlatchunks
            latpairs = []
            lonpairs = []
            for ll in range(lats[0],lats[1],dlats):
                latpairs.append([ll,ll+dlats])
            for ll in range(lons[0],lons[1]):
                lonpairs.append([ll,ll+1])
      else:
          mode="single"
      if mode!="single":
          for n in range(len(latpairs)):
              for l in range(len(latpairs[n])):
                  latpairs[n][l] = min(nlats-1,latpairs[n][l])
          for n in range(len(lonpairs)):
              for l in range(len(lonpairs[n])):
                  lonpairs[n][l] = min(nlons-1,lonpairs[n][l])
  else:
      mode="single"
      
  fakecloudnz = None
  fakecloudlwp = 0.0
                                                       
  smooth=False                                        
                                                       
  flat=False
    
  clouds=True
    
  notify = 'ae'
  
  if highcadence:
      isteps = range(len(data.variables['time'][:]))
      ntimes = '{'+','.join(np.array(isteps).astype(str))+'}'
  else:
      isteps= [12,]
      ntimes = '{'+str(isteps[0])+",}"
      
  uniform=False     
  unialb = 0.35

  wmin=0.35 #microns
  
  wmax = 80.0 #microns
  
  
  zout = (0.,100.)
      
  nlats = len(data.variables['lat'][:])
  
  os.system("mkdir "+workdir+"/source/")
  #os.system("tar cvzf "+job.top+"/sbdart_earth/"+source+"/source.tar.gz "+
            #job.top+"/sbdart_earth/"+source+"/* ")
  os.system("rsync "+top+"/sbdart_locked/"+source+"/source.tar.gz "+workdir+"/source/")
  #os.system("rm "+job.top+"/sbdart_earth/"+source+"/source.tar.gz ")
  os.system("tar xvzf "+workdir+"/source/source.tar.gz -C "+workdir+"/source/")
  os.system("rm "+workdir+"/source/source.tar.gz")
  if star:
      #os.system("cp -r "+finaldest+"/"+star+" "+workdir+"/source/solar.dat")
      w,f = np.loadtxt(finaldest+"/"+star,unpack=True)
      ffac = flux/np.trapz(f,x=w)
      f*=ffac
      data = np.array([w,f])
      data = data.T
      np.savetxt(workdir+"/source/solar.dat",data,fmt=['%f','%f'])
  
  for jlat in range(lats[0],lats[1]):
    for jlon in range(lons[0],lons[1]):
      for istep in isteps:
        for vw in vws:
            #print "Lat %02d Lon %02d Angle %1d View %s"%(jlat,jlon,nang,vw)
            #print "mkdir "+workdir+"/sbdart-%02d_%02d_%1d_%s"%(jlat,jlon,nang,vw)
            os.system("mkdir "+workdir+"/sbdart-%02d_%02d_%03d_%s"%(jlat,jlon,istep,vw))
            print "made directory "+workdir+"/sbdart-%02d-%02d-%03d_%s"%(jlat,jlon,istep,vw)
            os.system("cp -r "+workdir+"/source/* "+workdir+"/sbdart-%02d_%02d_%03d_%s/"%(jlat,jlon,istep,vw))

  icefile='seaice.dat'
  waterfile='seawater.dat'
  #if "icespec" in job.parameters:
  #    icefile=job.parameters["icespec"]
  #if "waterspec" in job.parameters:
  #    waterfile=job.parameters["waterspec"]
  
  
  tag = 'color map '
  
  runscript = ("#!/bin/bash \n\n"+
              "cd "+workdir+" \n"+
          "for al in "+ntimes+";                         \n"+
          "do \n"+
          "     for vw in "+lviews+";                \n"+
          "     do \n"+     
          "          for ((jl=$1; jl<=$2; jl++));                       \n"+
          "          do \n"+
          "               for ((il=$3; il<=$4; il++));                  \n"+
          "               do \n"+
          "                    ILAT=`printf '%02d' $(( 10#$jl ))`           \n"+
          "                    ILON=`printf '%02d' $(( 10#$il ))`           \n"+
          "                    IANG=`printf '%03d' $al`           \n"+
          "                    echo $ILAT $ILON $IANG $vw              \n"+
          "                    TAG=${ILAT}_${ILON}_${IANG}_$vw                    \n"+
          "                    cd "+workdir+"/sbdart-$TAG                          \n"+
          "                    ./sbdart > "+workdir+"/output/sbout.$TAG            \n"+
          "                    cd "+workdir+"                                  \n"+
          "               done                       \n"+
          "          done                         \n"+
          "          mv "+workdir+"/output/sbout.* "+finaldest+"/                   \n"+
          "     done                                    \n"+
          "done \n")
  with open(workdir+"/sbdart_batch.sh","w") as rs:
      rs.write(runscript)
  os.system("chmod a+x "+workdir+"/sbdart_batch.sh")
  jobscript =("#!/bin/bash -l \n"+
              "cd "+workdir+" \n"+
              "#ls -d * \n"+
              "cat sbdart_batch.sh  \n\n")
  for ltp in latpairs:
      for lnp in lonpairs:
          jobscript += "./sbdart_batch.sh %2d %2d %2d %2d &    \n\n"%(ltp[0],ltp[1],lnp[0],lnp[1])
  jobscript += ("wait     \n"+
          "cp "+finaldest+"/running/"+token_name+" "+finaldest+"/finished/ \n"+
          './release.sh "'+finaldest+'"                                \n'+
          "rm -rf "+workdir+"/                                  \n"+
          "cd "+finaldest+"  \n"+
          bts.SUB+" runparfix  \n")
  
  
  rs = open(workdir+"/runsbdart","w")
  rs.write(jobscript)
  rs.close()
  
  jobtag = ' '.join(job)
  dummyjob = crd.Job("# PID MODEL JOBNAME STATE NCORES QUEUE","%d pipeline pfix_%s 0 %d %s"%(-9999,jobname,jobncores,jobqueue),-1)
  fixscript = (bts.BATCHSCRIPT(dummyjob,'ae')+
               "cd "+top+"/sbdart_locked/     \n"+
               "python -B buildsbdart.py PARFIX %s  \n"%jobtag)
  with open(finaldest+"/runparfix","w") as fixf:
      fixf.write(fixscript)
  
  jobtag = ' '.join(job)
  dummyjob = crd.Job("# PID MODEL JOBNAME STATE NCORES QUEUE","%d pipeline fix_%s 0 1 %s"%(-9999,jobname,jobqueue),-1)
  fixscript = (bts.BATCHSCRIPT(dummyjob,'ae')+
               "cd "+top+"/sbdart_locked/     \n"+
               "python -B buildsbdart.py FIX %s  \n"%jobtag+
               "cd "+top+"                    \n"+
               "python -B setpostprocess_locked.py "+top+"/sbdart_locked/"+jobname+" "+
                jobname+" "+'^'.join(ntimes[1:-1].split(','))+" "+'^'.join(lviews[1:-1].split(','))+" "+jobqueue+" \n")
  with open(finaldest+"/runfix","w") as fixf:
      fixf.write(fixscript)
  
  os.system("chmod a+x "+workdir+"/runsbdart")
  os.system("cp "+top+"/crawldefs.py "+workdir+"/")
  os.system("cp "+top+"/identity.py "+workdir+"/")
  
  #views = [28.125,118.125,196.875,275.625]
  views = [0.,90.,180.,270.]
  nang = itimes[0]
  
  for jlon in range(lons[0],lons[1]):
    for jlat in range(lats[0],lats[1]):
      for istep in isteps:
        csz,azm,surf,sic,tsurf,altz = analyzecell_plasim_locked(data,vws,jlat,jlon,
                                                         workdir+"/sbdart-%02d_%02d_%03d"%(jlat,jlon,istep),
                                                         sol_lon = views[nang],
                                                         grav=grav,smooth=smooth,clouds=clouds,
                                                         istep=istep)
        if uniform:
            surf='uniform'
        latitude = data.variables['lat'][jlat]
        longitude = data.variables['lon'][jlon]
        #vws = ['Z','N','E','S','W']
        for vv in range(0,len(vws)):
            nv = viewdict[vws[vv]]
            write_input_locked(workdir+"/sbdart-%02d_%02d_%03d_%s"%(jlat,jlon,istep,vws[vv]),
                              nv,nang,csz,azm,latitude,longitude,surf,pco2,p0,tsurf,altz,flux,
                              wmin=wmin,wmax=wmax,albedo=unialb,flat=flat,sic=sic,spec=star,
                              smooth=smooth,iout=iout,zout=zout,waterfile=waterfile,
                              icefile=icefile)
        print "Prepped lat %02d lon %02d Time %03d View %s"%(jlat,jlon,istep,vws[vv])
        
def _run(job):
  jobname = job[0]
  workdir = ide.SCRATCH+"/sbdart_locked_"+jobname
  os.system("cd "+workdir+" && bash runsbdart")

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
  
def getbroken(job):
    name = job[0]
    views = job[8].split('^')
    
    if views[-1]=='':
        views = views[:-1]
    
    if job[9]=="True":
        highcadence=True
    else:
        highcadence=False
        
    with open("../.home","r") as homef:
      top = homef.read().split('\n')[0]
  
    namedir = top+"/sbdart_locked/%s"%name 
    
    data = nc.Dataset(top+"/plasim/output/"+name+"_snapshot.nc","r")
    
    if highcadence:
        isteps = range(len(data.variables['time'][:]))
    else:
        isteps= [12,]
    
    nlats = len(data.variables['lat'][:])
    nlons = len(data.variables['lon'][:])
    
    ofiles = sorted(glob.glob(namedir+"/sbout*"))
    broken = []
    reasons = []
    rejigger = []
    
    for istep in isteps:
        for nlt in range(0,nlats):
            for nln in range(0,nlons):
                for v in views:
                    if namedir+"/sbout.%02d_%02d_%03d_%s"%(nlt,nln,istep,v) not in ofiles:
                        broken.append("sbout.%02d_%02d_%03d_%s"%(nlt,nln,istep,v))
                        reasons.append("missing file")
                        rejigger.append(False)
    for o in ofiles:
        size = os.path.getsize(o)
        if size<2.45e4:
            broken.append(o.split('/')[-1])
            reasons.append("incomplete file")
            rejigger.append(False)
        else:
            try:
                bundle = readradiance(o)
                if np.nanmax(bundle["radiance"])<1.0 and bundle["zens"]<90.:
                    broken.append(o.split('/')[-1])
                    reasons.append("unphysically low output")
                    rejigger.append(False)
            except: #We encounter some sort of error when parsing the file
                broken.append(o.split('/')[-1])
                reasons.append("parsing error")
                with open(o,"r") as tmpo:
                    ddt = tmpo.read()
                if "gwk" in ddt:
                    rejigger.append(True)
                else:
                    rejigger.append(False)
    lats = []
    lons = []
    steps = []
    vws = []
    kb = 0
    for b in broken:
        tag = b.split('sbout')[1][1:]
        parts = tag.split('_')
        if parts[0]!=str(nlats) and parts[1]!=str(nlons):
            lats.append(int(parts[0]))
            lons.append(int(parts[1]))
            steps.append(int(parts[2]))
            vws.append(parts[3])
            print "Found broken output file: sbout_%02d_%02d_%03d_%s --- reason: %s"%(lats[-1],lons[-1],steps[-1],vws[-1],reasons[kb])
        else:
            os.system("rm "+top+"/sbdart_locked/%s/%s"%(name,b))
        kb += 1
    #if len(steps)>0:
        #return lons,lats,steps[-1],vws,rejigger  
    #else:
        #return lons,lats,steps,vws,rejigger
    return lons,lats,steps,vws,rejigger


def par_do_plasim_locked(job,vws,isteps,lons,lats,rejigger):
  name = job[0]
  with open("../.home","r") as homef:
    top = homef.read().split('\n')[0]
  
  jobname = job[0]
  jobncores = int(job[1])
  jobqueue = job[2]
  pco2 = float(job[3])
  flux = float(job[4])
  p0 = float(job[5])
  grav = float(job[6])
  ntimes = job[7].split('^')
  views = job[8].split('^')
  if job[9]=="True":
      highcadence=True
  else:
      highcadence=False
  try:
      star = job[10]
  except:
      star = False
  
  viewdict = {'Z':0,
              'E':1,
              'W':2,
              'N':3,
              'S':4}    
  
  nang = 0
    
  print "Setting up for times ",nang,"and views",vws
  
  outtag = "snapshot"
  if highcadence:
      outtag="highcadence"
  
  
  
  if len(lons)>=jobncores:
      pjobscript = ("#!/bin/bash -l           \n\n"+
                    "cd "+top+"/sbdart_locked/%s/   \n"%name)
      
      ncells = int(len(lons)/jobncores+0.5)
      subjobs = []
      jobtag = ' '.join(job)
      
      
      #We're going to spawn NCORES independent 1-core jobs, each of which will have a
      #specialized list of things to do.
      
      for j in range(jobncores-1):
        
        dummyjob = crd.Job("# PID MODEL JOBNAME STATE NCORES QUEUE","%d pipeline %s 0 1 %s"%(-9999,"fixsbdart_%d_%s"%(j,name),jobqueue),-1)
        jobscript = bts.BATCHSCRIPT(dummyjob,'a')
        jobtask = ''
        
        #Build the list of cells this job will have to do
        for n in range(ncells):
            jlat = lats[j*ncells+n]
            jlon = lons[j*ncells+n]
            vw = vws[j*ncells+n]
            istep = isteps[j*ncells+n]
            jobtask += "%d %d %d %s %d\n"%(jlat,jlon,istep,vw,rejigger[n]*1.0)
            
        jobscript += "cd "+top+"/sbdart_locked/  \n"
        
        jobscript += "python -B buildsbdart.py BATCHFIX %d %s  \n"%(j,jobtag)
        
        #Write the control script that will run these cells
        with open(top+"/sbdart_locked/%s/fixsbdart_%d.sh"%(name,j),"w") as rf:
            rf.write(jobscript)
        
        #Write the list of cells for this job
        with open(top+"/sbdart_locked/%s/laundry%d"%(name,j),"w") as jf:
            jf.write(jobtask)
        
        #This will submit and spawn the job, and store the job id in '$JOB%d'%j
        pjobscript += "JOB%d=$("%j+bts.SUB+" fixsbdart_%d.sh)     \n"%j
        subjobs.append("$JOB%d"%j)
        
      #We couldn't do this one in the loop because it may have fewer elements  
      j=jobncores-1
      dummyjob = crd.Job("# PID MODEL JOBNAME STATE NCORES QUEUE","%d pipeline %s 0 1 %s"%(-9999,"fixsbdart_%d_%s"%(j,name),jobqueue),-1)
      jobscript = bts.BATCHSCRIPT(dummyjob,'a')
      jobtask = ''
      
      for n in range(j*ncells,len(lons)): #may be less than ncells
          jlat = lats[n]
          jlon = lons[n]
          vw = vws[n]
          istep = isteps[n]
          jobtask += "%d %d %d %s %d\n"%(jlat,jlon,istep,vw,rejigger[n]*1.0)
          
      jobscript += "cd "+top+"/sbdart_locked/  \n"
      jobscript += "python -B buildsbdart.py BATCHFIX %d %s  \n"%(j,jobtag)
      
      with open(top+"/sbdart_locked/%s/fixsbdart_%d.sh"%(name,j),"w") as rf:
          rf.write(jobscript)
      with open(top+"/sbdart_locked/%s/laundry%d"%(name,j),"w") as jf:
          jf.write(jobtask)
          
      pjobscript += "JOB%d=$("%j+bts.SUB+" fixsbdart_%d.sh)     \n\n"%j
      subjobs.append("$JOB%d"%j)
      
      #Now that we've spawned all child jobs, submit the follow-up job, but specifying that
      #it should be held until the previous jobs finish.
      pjobscript += bts.HOLD(subjobs)+" runfix  \n"
      with open(top+"/sbdart_locked/%s/fixsbdart.sh"%name,"w") as prf:
          prf.write(pjobscript)
  else:    
      #We have fewer than NCORES cells to fix, so we'll do the same thing, but have fewer
      #than NCORES jobs submitted, with one cell per job.
      pjobscript = ("#!/bin/bash -l           \n\n"+
                    "cd "+top+"/sbdart_locked/%s/   \n"%name)
      subjobs = []
      jobtag = ' '.join(job)
      for j in range(len(lats)):
        
        dummyjob = crd.Job("# PID MODEL JOBNAME STATE NCORES QUEUE","%d pipeline %s 0 1 %s"%(-9999,"fixsbdart_%d_%s"%(j,name),jobqueue),-1)
        jobscript = bts.BATCHSCRIPT(dummyjob,'a')
        jobtask = ''
        jlat = lats[j]
        jlon = lons[j]
        vw = vws[j]
        istep = isteps[j]
        
        jobtask += "%d %d %d %s %d\n"%(jlat,jlon,istep,vw,rejigger[j]*1.0)
    
        jobscript += "cd "+top+"/sbdart_locked/  \n"
        jobscript += "python -B buildsbdart.py BATCHFIX %d %s  \n"%(j,jobtag)
        
        with open(top+"/sbdart_locked/%s/fixsbdart_%d.sh"%(name,j),"w") as rf:
            rf.write(jobscript)
            
        with open(top+"/sbdart_locked/%s/laundry%d"%(name,j),"w") as jf:
            jf.write(jobtask)
        
        pjobscript += "JOB%d=$("%j+bts.SUB+" fixsbdart_%d.sh)     \n"%j
        subjobs.append("$JOB%d"%j)
      
      #When those previous jobs are finished, submit the single-core fixer
      pjobscript += bts.HOLD(subjobs)+" runfix  \n"
      with open(top+"/sbdart_locked/%s/fixsbdart.sh"%name,"w") as prf:
          prf.write(pjobscript)
  
  #Run the script that will submit these child jobs and continue the work.
  os.system("chmod a+x "+top+"/sbdart_locked/%s/fixsbdart.sh"%name)
  os.system("bash "+top+"/sbdart_locked/%s/fixsbdart.sh"%jobname) 
  
def batch_plasim_locked(jid,job):
  #SBDART-process the cells given in laundry<jid>.
  name = job[0]
  workdir = ide.SCRATCH+"/fix_%s_%d"%(name,jid)
  os.system("mkdir "+workdir)
  with open("../.home","r") as homef:
    top = homef.read().split('\n')[0]
  
  jobname = job[0]
  jobncores = int(job[1])
  jobqueue = job[2]
  pco2 = float(job[3])
  flux = float(job[4])
  p0 = float(job[5])
  grav = float(job[6])
  views = job[8].split('^')
  ntimes = job[7].split('^')
  if job[9]=="True":
      highcadence=True
  else:
      highcadence=False
  try:
      star = job[10]
  except:
      star = False
    
  flat=False  
    
  lats = []
  lons = []
  vws = []
  steps = []
  rejigger = []
  #Get the list of cells we need to do
  with open(top+"/sbdart_locked/%s/laundry%d"%(name,jid),"r") as laundry:
      jlist = laundry.read().split('\n')
      if jlist[-1]=="":
          jlist = jlist[:-1]
      for j in jlist:
          parts = j.split()
          lats.append(int(parts[0]))
          lons.append(int(parts[1]))
          steps.append(int(parts[2]))
          vws.append(parts[3])
          rejigger.append(int(parts[4])>0.5)

  viewdict = {'Z':0,
              'E':1,
              'W':2,
              'N':3,
              'S':4}    
  
  print "Setting up for times ",steps,"and views",vws
  
  outtag = "snapshot"
  if highcadence:
      outtag="highcadence"
  
    
  os.system("cp "+top+"/plasim/output/"+jobname+"_%s.nc "%outtag+
            workdir+"/"+name+"_%s.nc"%outtag)
  data = nc.Dataset(workdir+"/"+name+"_%s.nc"%outtag,"r")  
  
  nlats = len(data.variables['lat'][:])
  
  if highcadence:
      isteps = range(len(data.variables['time'][:]))
      ntimes = '{'+','.join(np.array(isteps).astype(str))+'}'
  else:
      isteps= [12,]
      ntimes = '{'+str(isteps[0])+",}"
    
  os.system("mkdir "+workdir+"/source/")

  os.system("rsync "+top+"/sbdart_locked/clean/source.tar.gz "+workdir+"/source/")
  os.system("tar xvzf "+workdir+"/source/source.tar.gz -C "+workdir+"/source/")
  os.system("rm "+workdir+"/source/source.tar.gz")

  if star:
      #os.system("cp -r "+top+"/sbdart_locked/%s/"%jobname+star+" "+workdir+"/source/solar.dat")
      w,f = np.loadtxt(top+"/sbdart_locked/%s/"%jobname+star,unpack=True)
      ffac = flux/np.trapz(f,x=w)
      f*=ffac
      data = np.array([w,f])
      data = data.T
      np.savetxt(workdir+"/source/solar.dat",data,fmt=['%f','%f'])
  
      
     
  jobscript = "#!/bin/bash -l                                                  \n\n"
  
  #Set up the script that will iterate over and process the cells
  n=0
  for jlon in lons:
    jlat = lats[n]
    vw = vws[n]
    istep = steps[n]
    n+=1
    #print "Lat %02d Lon %02d Angle %1d View %s"%(jlat,jlon,nang,vw)
    print "mkdir "+workdir+"/sbdart-%02d_%02d_%03d_%s"%(jlat,jlon,istep,vw)
    os.system("mkdir "+workdir+"/sbdart-%02d_%02d_%03d_%s"%(jlat,jlon,istep,vw))
    os.system("cp -r "+workdir+"/source/* "+workdir+"/sbdart-%02d_%02d_%03d_%s/"%(jlat,jlon,istep,vw))
    
    tag = "%02d_%02d"%(jlat,jlon)
    
    jobscript += ("cd "+workdir+"/sbdart-%02d_%02d_%03d_%s    \n"%(jlat,jlon,istep,vw)+
                "./sbdart > "+workdir+"/sbout.%02d_%02d_%03d_%s   \n"%(jlat,jlon,istep,vw))
    
  jobscript += ("cp "+workdir+"/sbout* "+top+"/sbdart_locked/%s/     \n"%jobname+
                "rm -rf "+workdir+"/ \n")
  
  with open(top+"/sbdart_locked/%s/subjob_%d.sh"%(name,jid),"w") as rf:
      rf.write(jobscript)
  os.system("chmod a+x "+top+"/sbdart_locked/%s/subjob_%d.sh"%(name,jid))
  
  views = [0.,90.,180.,270.]
  
  nang = 0
  #Prepare the boundary data for each cell
  n=0
  for jlon in lons:
    jlat = lats[n]
    vw = vws[n]
    istep = steps[n]
    offset=0
    yoffset=0
    if rejigger[n]: #Sometimes we get a failure mode at the terminator at the poles
        os.system("touch rejigger_%02d_%02d_%03d_%s"%(jlat,jlon,istep,vw))
        if jlon>0:
            offset = -1
        else:
            offset =  1
        if jlat>0:
            yoffset = -1
        else:
            yoffset = 1
    n+=1
    csz,azm,surf,sic,tsurf,altz = analyzecell_plasim_locked(data,[vw,],jlat+yoffset,jlon+offset,
                                                     workdir+"/sbdart-%02d_%02d_%03d"%(jlat,jlon,istep),istep=istep,
                                                     sol_lon = views[nang])
    latitude = data.variables['lat'][jlat+yoffset]
    longitude = data.variables['lon'][jlon+offset]
    #vws = ['Z','N','E','S','W']
    nv = viewdict[vw]
    write_input_locked(workdir+"/sbdart-%02d_%02d_%03d_%s"%(jlat,jlon,istep,vw),
                    nv,nang,csz,azm,latitude,longitude,surf,pco2,p0,tsurf,altz,flux,
                    wmin=0.35,wmax=80.0,sic=sic,spec=star,flat=flat)
    print "Prepped lat %02d lon %02d Time %03d View %s"%(jlat,jlon,istep,vw)

  #Launch the script that iterates over these cells.
  os.system("bash "+top+"/sbdart_locked/%s/subjob_%d.sh"%(jobname,jid)) 


def do_plasim_locked(job,vws,isteps,lons,lats,rejigger):
  name = job[0]
  workdir = ide.SCRATCH+"/fix_%s"%name
  os.system("mkdir "+workdir)
  with open("../.home","r") as homef:
    top = homef.read().split('\n')[0]
  
  jobname = job[0]
  jobncores = int(job[1])
  jobqueue = job[2]
  pco2 = float(job[3])
  flux = float(job[4])
  p0 = float(job[5])
  grav = float(job[6])
  views = job[8].split('^')
  ntimes = job[7].split('^')
  if job[9]=="True":
      highcadence=True
  else:
      highcadence=False
  try:
      star = job[10]
  except:
      star = False
  
  viewdict = {'Z':0,
              'E':1,
              'W':2,
              'N':3,
              'S':4}    
  
  flat = False
    
  print "Setting up for times ",isteps,"and views",vws
  
  outtag = "snapshot"
  if highcadence:
      outtag="highcadence"
  
  os.system("cp "+top+"/plasim/output/"+jobname+"_%s.nc "%outtag
            +workdir+"/"+name+"_%s.nc"%outtag)
  data = nc.Dataset(workdir+"/"+name+"_%s.nc"%outtag,"r")  
  
  nlats = len(data.variables['lat'][:])
  
  os.system("mkdir "+workdir+"/source/")
  #os.system("tar cvzf /mnt/scratch-lustre/paradise/crawler2/sbdart_earth/"+source+"/source.tar.gz "+
            #job.top+"/sbdart_earth/"+source+"/* ")
  os.system("rsync "+top+"/sbdart_locked/clean/source.tar.gz "+workdir+"/source/")
  #os.system("rm /mnt/scratch-lustre/paradise/crawler2/sbdart_earth/"+source+"/source.tar.gz ")
  os.system("tar xvzf "+workdir+"/source/source.tar.gz -C "+workdir+"/source/")
  os.system("rm "+workdir+"/source/source.tar.gz")

  if star:
      #os.system("cp -r "+top+"/sbdart_locked/%s/"%jobname+star+" "+workdir+"/source/solar.dat")
      w,f = np.loadtxt(top+"/sbdart_locked/%s/"%jobname+star,unpack=True)
      ffac = flux/np.trapz(f,x=w)
      f*=ffac
      data = np.array([w,f])
      data = data.T
      np.savetxt(workdir+"/source/solar.dat",data,fmt=['%f','%f'])
  
      
     
  jobscript = "#!/bin/bash -l                                                  \n\n"
  
  #os.system("rsync /mnt/scratch-lustre/paradise/crawler2/sbdart_locked/clean/source.tar.gz "+workdir+"/source/")
  #os.system("tar xvzf "+workdir+"/source/source.tar.gz -C "+workdir+"/source/")
  #os.system("rm "+workdir+"/source/source.tar.gz")\n")
  
  n=0
  for jlon in lons:
    jlat = lats[n]
    vw = vws[n]
    istep = isteps[n]
    n+=1
    #print "Lat %02d Lon %02d Angle %1d View %s"%(jlat,jlon,nang,vw)
    print "mkdir "+workdir+"/sbdart-%02d_%02d_%03d_%s"%(jlat,jlon,istep,vw)
    os.system("mkdir "+workdir+"/sbdart-%02d_%02d_%03d_%s"%(jlat,jlon,istep,vw))
    os.system("cp -r "+workdir+"/source/* "+workdir+"/sbdart-%02d_%02d_%03d_%s/"%(jlat,jlon,istep,vw))
    
    tag = "%02d_%02d"%(jlat,jlon)
    
    jobscript += ("cd "+workdir+"/sbdart-%02d_%02d_%03d_%s    \n"%(jlat,jlon,istep,vw)+
                "./sbdart > "+workdir+"/sbout.%02d_%02d_%03d_%s   \n"%(jlat,jlon,istep,vw))
    
  jobscript += ("cp "+workdir+"/sbout* "+top+"/sbdart_locked/%s/     \n"%jobname+
                "rm -rf "+workdir+"/ \n")
  
  with open(top+"/sbdart_locked/%s/fixsbdart.sh"%name,"w") as rf:
      rf.write(jobscript)
  os.system("chmod a+x "+top+"/sbdart_locked/%s/fixsbdart.sh"%name)
  
  views = [0.,90.,180.,270.]
  
  nang = 0
  
  n=0
  for jlon in lons:
    jlat = lats[n]
    vw = vws[n]
    istep = isteps[n]
    offset=0
    yoffset=0
    if rejigger[n]: #Sometimes we get a failure mode at the terminator at the poles
        os.system("touch rejigger_%02d_%02d_%03d_%s"%(jlat,jlon,istep,vw))
        if jlon>0:
            offset = -1
        else:
            offset =  1
        if jlat>0:
            yoffset = -1
        else:
            yoffset =  1
    n+=1
    csz,azm,surf,sic,tsurf,altz = analyzecell_plasim_locked(data,[vw,],jlat+yoffset,jlon+offset,
                                                     workdir+"/sbdart-%02d_%02d_%03d"%(jlat,jlon,istep),istep=istep,
                                                     sol_lon = views[nang])
    latitude = data.variables['lat'][jlat+yoffset]
    longitude = data.variables['lon'][jlon+offset]
    #vws = ['Z','N','E','S','W']
    nv = viewdict[vw]
    write_input_locked(workdir+"/sbdart-%02d_%02d_%03d_%s"%(jlat,jlon,istep,vw),
                    nv,nang,csz,azm,latitude,longitude,surf,pco2,p0,tsurf,altz,flux,
                    wmin=0.35,wmax=80.0,sic=sic,flat=flat,spec=star)
    print "Prepped lat %02d lon %02d Time %03d View %s"%(jlat,jlon,istep,vw)

  os.system("bash "+top+"/sbdart_locked/%s/fixsbdart.sh"%jobname) 
  lons,lats,isteps,views,rejigger = getbroken(job)
  if len(lons)>0:
      do_plasim_locked(job,views,isteps,lons,lats,rejigger) #Recursive--risky, but this should run
                                                   #we're all good.
  
if __name__=="__main__":
    if len(sys.argv[:])==1:
        job = np.load("jobdat.npy").item()
        mode = "crawler2"
        try:
            top = job.top
        except:
            os.system("echo 'CRITICAL ERROR OPENING jobdat.npy'>sbdartfail.log")
    else:
        if sys.argv[1]=="FIX":
            job = sys.argv[2:]
            mode = "fixup"
        elif sys.argv[1]=="PARFIX":
            job = sys.argv[2:]
            mode = "parfixup"
        elif sys.argv[1]=="BATCHFIX":
            jid = int(sys.argv[2])
            job = sys.argv[3:]
            mode = "batchfix"
        else:
            job = sys.argv[1:]
            mode = "command"
            
    if mode=="crawler2":
        prep(job)
        run(job)
    elif mode=="command":
        _prep(job)
        _run(job)
    elif mode=="fixup":
        name = job[0]
        lons,lats,isteps,views,rejigger = getbroken(job)
        do_plasim_locked(job,views,isteps,lons,lats,rejigger)
    elif mode=="parfixup":
        name = job[0]
        lons,lats,isteps,views,rejigger = getbroken(job)
        par_do_plasim_locked(job,views,isteps,lons,lats,rejigger)
    elif mode=="batchfix":
        name = job[0]
        #We don't need to get the broken jobs for this mode, because they
        #were gotten and stored in the parfixup mode that ran previously
        batch_plasim_locked(jid,job)

    #submit(job)