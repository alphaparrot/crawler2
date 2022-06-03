import matplotlib
matplotlib.use('Agg')
import numpy as np
import glob
import sys
from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc
import netCDF4 as netc
#import gcmtools as gt
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, LogLocator, NullFormatter)

# Usage: python plasimtransit.py output.nc xH2 xHe xCO2 xN2 xO2 gascon
#           NOTE: gascon is optional, and will override the mmw computed from abundances, but
#                 modification from humidity will still be applied.



if __name__=="__main__":

    smws = {'mH2': 2.01588,
          'mHe': 4.002602,
          'mN2': 28.0134,
          'mO2': 31.9988,
          'mCO2':44.01,
          'mAr': 39.948,
          'mNe': 20.1797,
          'mKr': 83.798,
          'mH2O':18.01528}
	
    
    atmosphere = Radtrans(line_species = ['H2O', 'CO2'], \
                        cloud_species = ['H2O(c)_cm'], \
                        rayleigh_species = ['H2','He','N2','O2'], \
                        continuum_opacities = ['H2-H2', 'H2-He'], \
                        #wlen_bords_micron = [0.3, 15])
                        wlen_bords_micron = [0.4, 6])
    ugascon=None
    output = sys.argv[1]
    plarad = float(sys.argv[2])
    grav = float(sys.argv[3])
    starrad = float(sys.argv[4])
    xH2 = float(sys.argv[5])
    xHe = float(sys.argv[6])
    xCO2 = float(sys.argv[7])
    xN2 = float(sys.argv[8])
    xO2 = float(sys.argv[9])
    try:
        ugascon = float(sys.argv[10])
    except:
        ugascon = None
    #surfp = gt.spatialmath("ps",file=output)
    ncfile = netc.Dataset(output,"r")
    ln  = ncfile.variables['lon']
    lt  = ncfile.variables['lat']
    lev = ncfile.variables['lev']
    ts  = ncfile.variables['ts']
    prw = ncfile.variables['prw']
    ps  = ncfile.variables['ps']
    clt = ncfile.variables['clt']
    clf = ncfile.variables['cl']
    sic = ncfile.variables['sic']
    hus = ncfile.variables['hus']
    ta  = ncfile.variables['ta']
    dlat = abs(np.gradient(lt))
    latwt = 0.0
    for lat in range(len(lt)):
        latwt += 2*dlat[lat]
    transit = np.zeros(2708)
    prefix = output[:-3]
    if "notime" in sys.argv[:]:
        ntimes=[3,]
    else:
        ntimes=list(range(ps.shape[0]))
    alltransits = np.zeros((len(ntimes)*len(lt)*2,2708))
    snaptransits = np.zeros((len(ntimes),2708))
    nx=0
    lon=16
    for lat in range(len(lt)):
        print("Latitude %d"%(lat+1))
        for nt in ntimes:
            print("Time %d"%(nt+1))
            p_sim=ps[nt,lat,lon]*lev
            t_sim=ta[nt,:,lat,lon]
            h2o_spechum_sim=hus[nt,:,lat,lon]
            cloud = clf[nt,:,lat,lon]   
            if len(lev) == 10:
               m1 = np.log10((h2o_spechum_sim[1]+1.0e-18)/(h2o_spechum_sim[0]+1.0e-18))\
                        /np.log10((p_sim[1]+1.0e-18)/(p_sim[0]+1.0e-18))  
               new_p = np.concatenate((np.logspace(-1,np.log10(p_sim[0]),num=50),
                           np.logspace(np.log10(p_sim[0]),np.log10(p_sim[1]),num=50)[1:],
                           p_sim[2:]))
               new_q = np.zeros_like(new_p)
               new_q[:50] = np.minimum(h2o_spechum_sim[0]*(new_p[:50]/p_sim[0])**m1,h2o_spechum_sim[0])
               new_q[50:98] = np.minimum(h2o_spechum_sim[1]*(new_p[50:98]/p_sim[1])**m1,h2o_spechum_sim[1])
                   #new_q[50:] = np.interp(new_p[50:],p_sim,h2o_spechum_sim)
               new_q[98:] = h2o_spechum_sim[1:]
               
               new_t = np.zeros_like(new_p)
               new_t[:50] = t_sim[0]
               new_t[50:] = np.interp(new_p[50:],p_sim,t_sim)
               
               new_cloud = np.zeros_like(new_p)
               new_cloud[50:] = np.interp(new_p[50:],p_sim,cloud)
                   #pressures = p_sim / 1e3 # mbar to bar
            else:
               new_p = p_sim[:]
               new_q = h2o_spechum_sim[:]
               new_t = t_sim[:]
               new_cloud = cloud[:]
            pressures = new_p / 1e3
                
            atmosphere.setup_opa_structure(pressures)
                
            R_pl = plarad* nc.r_earth 
            gravity = grav * 1e2 # m to cm
            P0 = ps[nt,lat,lon]*1.0e-3 #0.1

                #temperature = t_sim
            temperature = new_t
                
            gasesx={}
            gasesx['H2'] = xH2*(1-new_q)
            gasesx['He'] = xHe*(1-new_q)
            gasesx['CO2'] = xCO2*(1-new_q)
            gasesx['N2'] = xN2*(1-new_q)
            gasesx['O2'] = xO2*(1-new_q)
            gasesx['H2O'] = new_q
            totalfrac = np.zeros(new_q.shape)
            for z in list(gasesx.keys()):
                totalfrac += gasesx[z]
            if totalfrac.min()==0. or np.isnan(totalfrac.min()):
                print(totalfrac)
                print(gasesx)
            #totalfrac = np.nansum([gasesx[z] for z in gasesx.keys()],axis=0)
            for z in list(gasesx.keys()):
                if not (gasesx[z].max()==gasesx[z].min()==0.):
                    gasesx[z] = gasesx[z]/totalfrac
            abundances = {}
            abundances['H2'] = gasesx['H2']
            abundances['He'] = gasesx['He']
            abundances['N2'] = gasesx['N2']
            abundances['O2'] = gasesx['O2']
                #abundances['H2O'] = h2o_vmr_sim * 20. / 2.33 #h2o_spechum_sim #*1e10
                #abundances['H2O'] = h2o_spechum_sim
            abundances['H2O'] = gasesx['H2O']
            #abundances['CO_all_iso'] = 0.0 * np.ones_like(temperature)
            abundances['CO2'] = gasesx['CO2']
            #abundances['CH4'] = 0.0 * np.ones_like(temperature)
            #abundances['Na'] = 0.0 * np.ones_like(temperature)
            #abundances['K'] = 0.0 * np.ones_like(temperature)
                #abundances['H2O(c)'] = 0.0 * np.ones_like(temperature)
                #abundances['H2O(c)'] = cloud*h2o_spechum_sim
            abundances['H2O(c)'] = new_cloud*new_q
                #abundances['H2O(c)'][1] = 1e6*abundances['H2O'][1]
            if ugascon:
                mgascon = 8314.46261815324/ugascon
                mmwd = (1-gasesx['H2O'])/mgascon + gasesx['H2O']/smws['mH2O'] #same shape as new_q
                MMW = 1.0/mmwd
            else:
                MMW = np.ones_like(temperature)
                for k in range(len(MMW)):
                    mmwd = 0
                    for x in list(gasesx.keys()):
                        mmwd += gasesx[x][k]/smws['m'+x]
                    mmwt = 1.0/mmwd
                    MMW[k] = mmwt

            Kzz = np.ones_like(temperature)* 0.0
            fsed = 2.
            radius = {}
            radius['H2O(c)'] = 0.0008*np.ones_like(temperature) # I.e. an 8-micron particle size (0.00005 cm)

            sigma_lnorm = 1.05
                
            atmosphere.calc_transm(temperature, abundances, gravity, MMW, R_pl=R_pl, P0_bar=P0, \
                            radius = radius, sigma_lnorm = sigma_lnorm)
            transit = np.nansum([transit,atmosphere.transm_rad*dlat[lat]],axis=0)
            alltransits[nx,:] = atmosphere.transm_rad[:]
            snaptransits[nt,:] = np.nansum([snaptransits[nt,:],atmosphere.transm_rad*dlat[lat]],axis=0)
            nx+=1

    lon=48
    for lat in range(len(lt))[::-1]:
        print("Latitude %d"%(lat+1))
        for nt in ntimes:
            print("Time %d"%(nt+1))
            p_sim=ps[nt,lat,lon]*lev
            t_sim=ta[nt,:,lat,lon]
            h2o_spechum_sim=hus[nt,:,lat,lon]
            cloud = clf[nt,:,lat,lon]   
            
            if len(lev) == 10:
                m1 = np.log10((h2o_spechum_sim[1]+1.0e-18)/(h2o_spechum_sim[0]+1.0e-18))\
                         /np.log10((p_sim[1]+1.0e-18)/(p_sim[0]+1.0e-18))  
                new_p = np.concatenate((np.logspace(-1,np.log10(p_sim[0]),num=50),
                            np.logspace(np.log10(p_sim[0]),np.log10(p_sim[1]),num=50)[1:],
                            p_sim[2:]))
                new_q = np.zeros_like(new_p)
                new_q[:50] = np.minimum(h2o_spechum_sim[0]*(new_p[:50]/p_sim[0])**m1,h2o_spechum_sim[0])
                new_q[50:98] = np.minimum(h2o_spechum_sim[1]*(new_p[50:98]/p_sim[1])**m1,h2o_spechum_sim[1])
                    #new_q[50:] = np.interp(new_p[50:],p_sim,h2o_spechum_sim)
                new_q[98:] = h2o_spechum_sim[1:]
                
                new_t = np.zeros_like(new_p)
                new_t[:50] = t_sim[0]
                new_t[50:] = np.interp(new_p[50:],p_sim,t_sim)
                
                new_cloud = np.zeros_like(new_p)
                new_cloud[50:] = np.interp(new_p[50:],p_sim,cloud)
                    #pressures = p_sim / 1e3 # mbar to bar
            else:
                new_p = p_sim[:]
                new_q = h2o_spechum_sim[:]
                new_t = t_sim[:]
                new_cloud = cloud[:]
            pressures = new_p / 1e3
                
            atmosphere.setup_opa_structure(pressures)
                
            R_pl = plarad* nc.r_earth 
            gravity = grav * 1e2 # m to cm
            P0 = ps[nt,lat,lon]*1.0e-3 #0.1

                #temperature = t_sim
            temperature = new_t
                
            gasesx={}
            gasesx['H2'] = xH2*(1-new_q)
            gasesx['He'] = xHe*(1-new_q)
            gasesx['CO2'] = xCO2*(1-new_q)
            gasesx['N2'] = xN2*(1-new_q)
            gasesx['O2'] = xO2*(1-new_q)
            gasesx['H2O'] = new_q
            totalfrac = np.zeros(new_q.shape)
            for z in list(gasesx.keys()):
                totalfrac += gasesx[z]
            if totalfrac.min()==0. or np.isnan(totalfrac.min()):
                print(totalfrac)
                print(gasesx)
            #totalfrac = np.nansum([gasesx[z] for z in gasesx.keys()],axis=0)
            for z in list(gasesx.keys()):
                if not (gasesx[z].max()==gasesx[z].min()==0.):
                    gasesx[z] = gasesx[z]/totalfrac
            abundances = {}
            abundances['H2'] = gasesx['H2']
            abundances['He'] = gasesx['He']
            abundances['N2'] = gasesx['N2']
            abundances['O2'] = gasesx['O2']
                #abundances['H2O'] = h2o_vmr_sim * 20. / 2.33 #h2o_spechum_sim #*1e10
                #abundances['H2O'] = h2o_spechum_sim
            abundances['H2O'] = gasesx['H2O']
            #abundances['CO_all_iso'] = 0.0 * np.ones_like(temperature)
            abundances['CO2'] = gasesx['CO2']
            #abundances['CH4'] = 0.0 * np.ones_like(temperature)
            #abundances['Na'] = 0.0 * np.ones_like(temperature)
            #abundances['K'] = 0.0 * np.ones_like(temperature)
                #abundances['H2O(c)'] = 0.0 * np.ones_like(temperature)
                #abundances['H2O(c)'] = cloud*h2o_spechum_sim
            abundances['H2O(c)'] = new_cloud*new_q
                #abundances['H2O(c)'][1] = 1e6*abundances['H2O'][1]
            if ugascon:
                mgascon = 8314.46261815324/ugascon
                mmwd = (1-gasesx['H2O'])/mgascon + gasesx['H2O']/smws['mH2O'] #same shape as new_q
                MMW = 1.0/mmwd
            else:
                MMW = np.ones_like(temperature)
                for k in range(len(MMW)):
                    mmwd = 0
                    for x in list(gasesx.keys()):
                        mmwd += gasesx[x][k]/smws['m'+x]
                    mmwt = 1.0/mmwd
                    MMW[k] = mmwt

            Kzz = np.ones_like(temperature)* 0.0
            fsed = 2.
            radius = {}
            radius['H2O(c)'] = 0.0008*np.ones_like(temperature) # I.e. an 8-micron particle size (0.00005 cm)

            sigma_lnorm = 1.05
                
            atmosphere.calc_transm(temperature, abundances, gravity, MMW, R_pl=R_pl, P0_bar=P0, \
                            radius = radius, sigma_lnorm = sigma_lnorm)
            transit = np.nansum([transit,atmosphere.transm_rad*dlat[lat]],axis=0)
            alltransits[nx,:] = atmosphere.transm_rad[:]
            snaptransits[nt,:] = np.nansum([snaptransits[nt,:],atmosphere.transm_rad*dlat[lat]],axis=0)
            nx+=1
    transit /= (latwt*len(ntimes))
    snaptransits /= latwt
    np.save("wavelengths.npy",nc.c/atmosphere.freq/1e-4)
    np.save(prefix+"_transit.npy",transit)
    np.save(prefix+"_transits.npy",alltransits)
    np.save(prefix+"_transit_snapshots.npy",snaptransits)
    fig,ax=plt.subplots(figsize=(10,6))
    plt.plot(nc.c/atmosphere.freq/1e-4, (transit/(starrad*6.9551e10))**2*1e6, label = 'cloudy', zorder = 1)
    #plt.xlim(0.6,5.0)
    plt.xlim(0.4,6.0)
    #plt.ylim(2700,3000)
        # For the minor ticks, use no labels; default NullFormatter.
        #plt.xlim(1.0e-3,1.0e2)
    plt.xscale('log')
    ax.xaxis.set_major_locator(LogLocator(subs=np.linspace(1.0,9.5,num=18)))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%1.1f'))
        
    ax.xaxis.set_minor_locator(LogLocator(subs = np.linspace(1.0,9.9,num=90)))
    ax.xaxis.set_minor_formatter(NullFormatter())
    plt.xlabel('Wavelength [microns]')
        #plt.ylabel(r'Transit radius ($\rm R_{Earth}$)')
    plt.ylabel("Transit Depth [ppm]")
    plt.title(prefix+"; T_s = %d K"%t_sim[-1])
    plt.savefig(prefix+"_transit.png",bbox_inches='tight')
    plt.savefig(prefix+"_transit.pdf",bbox_inches='tight')
    plt.close('all')
    
    fig,ax=plt.subplots(figsize=(10,6))
    for nnx in range(nx):
        plt.plot(nc.c/atmosphere.freq/1e-4, (alltransits[nnx,:]/1e5), alpha=0.2,color='k', zorder = 1)
    #plt.xlim(0.6,5.0)
    plt.xlim(0.4,6.0)
    # For the minor ticks, use no labels; default NullFormatter.
    #plt.xlim(1.0e-3,1.0e2)
    plt.xscale('log')
    ax.xaxis.set_major_locator(LogLocator(subs=np.linspace(1.0,9.5,num=18)))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%1.1f'))
    
    ax.xaxis.set_minor_locator(LogLocator(subs = np.linspace(1.0,9.9,num=90)))
    ax.xaxis.set_minor_formatter(NullFormatter())
    plt.xlabel('Wavelength [microns]')
    #plt.ylabel(r'Transit radius ($\rm R_{Earth}$)')
    plt.ylabel("Transit Radius [km]")
    plt.title(prefix+"; T_s = %d K"%t_sim[-1])
    plt.savefig(prefix+"_transits.png",bbox_inches='tight')
    plt.savefig(prefix+"_transits.pdf",bbox_inches='tight')
    plt.close('all')
    
    fig,ax=plt.subplots(figsize=(10,6))
    for nnt in range(len(ntimes)):
        plt.plot(nc.c/atmosphere.freq/1e-4, (snaptransits[nnt,:]/(starrad*6.9551e10))**2*1e6, alpha=0.3,color='k', zorder = 1)
    #plt.xlim(0.6,5.0)
    plt.xlim(0.4,6.0)
    # For the minor ticks, use no labels; default NullFormatter.
    #plt.xlim(1.0e-3,1.0e2)
    plt.xscale('log')
    ax.xaxis.set_major_locator(LogLocator(subs=np.linspace(1.0,9.5,num=18)))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%1.1f'))
    
    ax.xaxis.set_minor_locator(LogLocator(subs = np.linspace(1.0,9.9,num=90)))
    ax.xaxis.set_minor_formatter(NullFormatter())
    plt.xlabel('Wavelength [microns]')
    #plt.ylabel(r'Transit radius ($\rm R_{Earth}$)')
    plt.ylabel("Transit Depth [ppm]")
    plt.title("Snapshots: "+prefix+"; T_s = %d K"%t_sim[-1])
    plt.savefig(prefix+"_transit_snapshots.png",bbox_inches='tight')
    plt.savefig(prefix+"_transit_snapshots.pdf",bbox_inches='tight')
    plt.close('all')
    
