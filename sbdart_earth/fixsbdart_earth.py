import numpy as np
import netCDF4 as nc
import os, sys, glob
import time

def write_input_earth(workdir,nview,ntime,cszenith,azimuth,latitude,longitude,surface,pCO2,p0,
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
  
  faces = [28.125,118.125,196.875,275.625]
  
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

def analyzecell_plasim_earth(data,views,lat,lon,workdir,grav=9.80665,sol_dec=0.0,
                     sol_lon=0.0,smooth=False,clouds=True,istep=-1,fakecloudnz=None,
                     fakecloudlwp=0.0,cirrus=False):
  #cszenith,azimuth,surface,pCO2,p0,tsurf,altz
  
  surf = "uniform"
  
  lsm = data.variables['lsm'][-1,lat,lon]
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
    if ((data.variables['snd'][istep,lat,lon]+data.variables['glac'][istep,lat,lon])>=0 and 
        data.variables['as'][istep,lat,lon]>0.2):
        surf = "snow"
    else:
        surf = "sand"
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
  
  print("Writing atms.dat")
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
    
  print("Writing usrcld.dat")
  #Need to implement a cloud-free option
  #Write usrcld.dat, which has level data on cloud water content and coverage fraction
  
  if type(fakecloudnz)!=type(None):
    dql[fakecloudnz] = fakecloudlwp
    cld[fakecloudnz] = 1.0
  
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
  
  print("Writing atms.dat")
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
    
  print("Writing usrcld.dat")
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
    if 'pctsrf_sic' in list(data.keys()):
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
  
  print("Writing atms.dat")
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
  
  print("Writing usrcld.dat")
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
    

def do_plasim_earth(name,vws,nang,pCO2,p0,flux,lons,lats): #data,lats,lons,pCO2,p0,flux,grav=9.80665
  workdir = "/mnt/node_scratch/paradise/fix_%s"%name
  os.system("mkdir "+workdir)
  
  viewdict = {'Z':0,
              'E':1,
              'W':2,
              'N':3,
              'S':4}    
      
  print("Setting up for times ",nang,"and views",vws)
  
  os.system("cp /mnt/scratch-lustre/paradise/crawler2/hopper/"+name+"_snapshot.nc "+workdir+"/"+name+"_snapshot.nc")
  data = nc.Dataset(workdir+"/"+name+"_snapshot.nc","r")  
  
  nlats = len(data.variables['lat'][:])
  
  os.system("mkdir "+workdir+"/source/")
  #os.system("tar cvzf /mnt/scratch-lustre/paradise/crawler2/sbdart_earth/"+source+"/source.tar.gz "+
            #job.top+"/sbdart_earth/"+source+"/* ")
  os.system("rsync /mnt/scratch-lustre/paradise/crawler2/sbdart_earth/clean/source.tar.gz "+workdir+"/source/")
  #os.system("rm /mnt/scratch-lustre/paradise/crawler2/sbdart_earth/"+source+"/source.tar.gz ")
  os.system("tar xvzf "+workdir+"/source/source.tar.gz -C "+workdir+"/source/")
  os.system("rm "+workdir+"/source/source.tar.gz")
  #if star:
      #os.system("cp -r /mnt/scratch-lustre/paradise/crawler2/sbdart_earth/"+star+" "+workdir+"/source/solar.dat")
  
      
  jobscript = "#!/bin/bash -l                                                  \n\n"
  
  #os.system("rsync /mnt/scratch-lustre/paradise/crawler2/sbdart_earth/clean/source.tar.gz "+workdir+"/source/")
  #os.system("tar xvzf "+workdir+"/source/source.tar.gz -C "+workdir+"/source/")
  #os.system("rm "+workdir+"/source/source.tar.gz")
  
  n=0
  for jlat in lats:
    jlon = lons[n]
    n+=1
    #print "Lat %02d Lon %02d Angle %1d View %s"%(jlat,jlon,nang,vw)
    #print "mkdir "+workdir+"/sbdart-%02d_%02d_%1d_%s"%(jlat,jlon,nang,vw)
    os.system("mkdir "+workdir+"/sbdart-%02d_%02d_%1d_%s"%(jlat,jlon,nang,vws[0]))
    os.system("cp -r "+workdir+"/source/* "+workdir+"/sbdart-%02d_%02d_%1d_%s/"%(jlat,jlon,nang,vws[0]))
    
    tag = "%02d_%02d"%(jlat,jlon)
    
    jobscript += ("cd "+workdir+"/sbdart-%02d_%02d_%1d_%s    \n"%(jlat,jlon,nang,vws[0])+
                  "./sbdart > "+workdir+"/sbout.%02d_%02d_%1d_%s   \n"%(jlat,jlon,nang,vws[0]))
      
  jobscript += ("cp "+workdir+"/sbout* /mnt/scratch-lustre/paradise/crawler2/sbdart_earth/%s/     \n"%name+
                "rm -rf "+workdir+"/ \n")
  
  with open("/mnt/scratch-lustre/paradise/crawler2/sbdart_earth/%s/fixsbdart.sh"%name,"w") as rf:
      rf.write(jobscript)
  os.system("chmod a+x /mnt/scratch-lustre/paradise/crawler2/sbdart_earth/%s/fixsbdart.sh"%name)
  
  views = [28.125,118.125,196.875,275.625]
  
  n=0
  for jlon in lons:
    jlat = lats[n]
    n+=1
    csz,azm,surf,sic,tsurf,altz = analyzecell_plasim_earth(data,vws,jlat,jlon,
                                                     workdir+"/sbdart-%02d_%02d_%1d"%(jlat,jlon,nang),
                                                     sol_lon = views[nang])
    latitude = data.variables['lat'][jlat]
    longitude = data.variables['lon'][jlon]
    #vws = ['Z','N','E','S','W']
    nv = viewdict[vws[0]]
    write_input_earth(workdir+"/sbdart-%02d_%02d_%1d_%s"%(jlat,jlon,nang,vws[0]),
                      nv,nang,csz,azm,latitude,longitude,surf,pCO2,p0,tsurf,altz,flux,
                      wmin=0.35,wmax=80.0)
    print("Prepped lat %02d lon %02d Angle %1d View %s"%(jlat,jlon,nang,vws[0]))

  os.system("bash /mnt/scratch-lustre/paradise/crawler2/sbdart_earth/%s/fixsbdart.sh"%name) 

  
def getbroken(name):
    namedir = "/mnt/scratch-lustre/paradise/crawler2/sbdart_earth/%s"%name 
    ofiles = glob.glob(namedir+"/sbout*")
    broken = []
    for o in ofiles:
        size = os.path.getsize(o)
        if size==0:
            broken.append(o.split('/')[-1])
    lats = []
    lons = []
    angs = []
    vws = []
    for b in broken:
        tag = b.split('sbout')[1][1:]
        parts = tag.split('_')
        if parts[0]!='32' and parts[1]!='64':
            lats.append(int(parts[0]))
            lons.append(int(parts[1]))
            angs.append(int(parts[2]))
            vws.append(parts[3])
            print("Found broken output file: sbout_%02d_%02d_%d_%s"%(lats[-1],lons[-1],angs[-1],vws[-1]))
        else:
            os.system("rm /mnt/scratch-lustre/paradise/crawler2/sbdart_earth/%s/%s"%(name,b))
    return lons,lats,angs[-1],vws[-1]
  
if __name__=="__main__":
    mode = sys.argv[1]
    if mode=="PREP":
        finput = sys.argv[2]
        with open(finput,"r") as rf:
            llist = rf.read().split('\n')[1:]
        for l in llist:
            if l!="":
                parts = l.split()
                name = parts[0]
                pco2 = parts[1]
                flux = parts[2]
                p0 = parts[3]
                
                runscript = ("#!/bin/bash -l                                                 \n"+
                             "#PBS -l nodes=1:ppn=1                                          \n"+
                             "#PBS -q workq                                                  \n"+
                             "#PBS -m a                                                      \n"+
                             "#PBS -r n                                                      \n"+
                             "#PBS -l walltime=48:00:00                                      \n"+
                             "#PBS -N "+name+"_fix                                           \n"+
                             "# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE"+
                             " nodes,ppn,walltime and my_job_name VALUES                     \n"+
                             "cd $PBS_O_WORKDIR                                              \n"+
                             "module load gcc/4.9.1      \n"+
                             "module load python/2.7.9   \n"+
                             "python fixsbdart_earth.py RUN %s %s %s %s \n"%(name,pco2,flux,p0))
                with open("runfix_%s"%name,"w") as rs:
                    rs.write(runscript)
                os.system("qsub runfix_%s"%name)
                
    elif mode=="RUN":
        name = sys.argv[2]
        pco2 = float(sys.argv[3])
        flux = float(sys.argv[4])
        p0 = float(sys.argv[5])
        lons,lats,nang,view = getbroken(name)
        do_plasim_earth(name,(view,),nang,pco2,p0,flux,lons,lats)
    #submit(job)