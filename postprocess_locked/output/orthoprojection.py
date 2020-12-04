import matplotlib
matplotlib.use('Agg')
import colormatch
import numpy as np
import matplotlib.pyplot as plt
import sys
import netCDF4 as nc
import os

def orthographic(lon,lat,imap,l0=0,p0=0,ny=36,nx=36,interp='bilinear'):
    xymap = np.zeros((ny,nx,3))
    xymap[:] = 0.0
    coords = np.zeros((ny,nx,2))
    coords2 = np.zeros((ny,nx,2,3))
    coords3 = np.zeros((ny,nx,2))
    rad=0.5*(8*nx/18.0+8*ny/18.0)
    
    zlat = np.zeros(np.array(lat.shape)+2)
    zlon = np.zeros(np.array(lon.shape)+2)
    dlat = np.diff(lat)
    dlon = np.diff(lon)
    zlon[0] = lon[0]-dlon[0]
    zlon[1:-1] = lon[:]
    zlon[-1] = lon[-1]+dlon[-1]
    zlat[0] = lat[0]-dlat[0]
    zlat[-1] = lat[-1]+dlat[-1]
    zlat[1:-1] = lat[:]
    zmap = np.zeros(np.array(imap.shape)+np.array((2,2,0)))
    zmap[1:-1,1:-1,:] = imap[:,:,:]
    zmap[0,:,:] = zmap[1,:,:]
    zmap[-1,:,:] = zmap[-2,:,:]
    zmap[:,0,:] = zmap[:,-2,:]
    zmap[:,-1,:] = zmap[:,1,:]
    
    if l0>180.0:
        l0 -= 360.0

    if lon.max()>180.0:
        zlon[zlon>180]-=360.0    

    p0 *= np.pi/180.0
    l0 *= np.pi/180.0
    xx = np.arange(nx)-nx/2
    yy = np.arange(ny)-ny/2
    for j in range(0,ny):
        jy = yy[j]
        for i in range(0,nx):
            ix = xx[i]
            if (ix**2+jy**2)<=rad**2:
                rho = np.sqrt(ix**2+jy**2)
                cc = np.arcsin(rho/rad)
                if rho==0:
                    phi = p0
                    lamb = l0
                else:
                    phi = np.arcsin(np.cos(cc)*np.sin(p0) + jy*np.sin(cc)*np.cos(p0)/rho)
                    lamb = l0 + np.arctan2(ix*np.sin(cc),(rho*np.cos(cc)*np.cos(p0)-jy*np.sin(cc)*np.sin(p0)))
                phi *= 180.0/np.pi
                lamb *= 180.0/np.pi
                if lamb<-180:
                    lamb += 360
                if lamb>180:
                    lamb -= 360
                jlat = np.argmin(abs(phi-zlat))
                jlon = np.argmin(abs(lamb-zlon))
                if interp=="bilinear":
                    xlat1 = np.where(zlat<phi)[0]
                    xlat2 = np.where(zlat>=phi)[0]
                    xlon1 = np.where(zlon<lamb)[0]
                    xlon2 = np.where(zlon>=lamb)[0]
                    if len(xlat1)>0 and len(xlat2)>0 and len(xlon1)>0 and len(xlon2)>0:
                        jlat1 = np.where(zlat==zlat[xlat1].max())[0][0]
                        jlat2 = np.where(zlat==zlat[xlat2].min())[0][0]
                        jlon1 = np.where(zlon==zlon[xlon1].max())[0][0]
                        jlon2 = np.where(zlon==zlon[xlon2].min())[0][0]
                        p1 = zlat[jlat1]
                        p2 = zlat[jlat2]
                        l1 = zlon[jlon1]
                        l2 = zlon[jlon2]
                        #print "interpolating;",jlat1,jlat2,jlon1,jlon2,lamb#,zlon[np.where(zlon<lamb)]
                        dl = l2-l1
                        dp = p2-p1
                        try:
                            if dl==0 and dp>0:
                                fxy1 = zmap[jlat1,jlon,:]
                                fxy2 = zmap[jlat2,jlon,:]
                                xymap[j,i,:] = (p2-phi)/dp*fxy1 + (phi-p1)/dp*fxy2
                            elif dp==0 and dl>0:
                                fx1y = zmap[jlat,jlon1,:]
                                fx2y = zmap[jlat,jlon2,:]
                                xymap[j,i,:] = (l2-lamb)/dl*fx1y + (lamb-l1)/dl*fx2y
                            elif dp==0 and dl==0:
                                xymap[j,i,:] = zmap[jlat,jlon,:]
                            else:
                                fxy1 = (l2-lamb)/dl*zmap[jlat1,jlon1,:] + (lamb-l1)/dl*zmap[jlat1,jlon2,:]
                                fxy2 = (l2-lamb)/dl*zmap[jlat2,jlon1,:] + (lamb-l1)/dl*zmap[jlat2,jlon2,:]
                                xymap[j,i,:] = (p2-phi)/dp*fxy1 + (phi-p1)/dp*fxy2
                        except:
                            print p1,p2,l1,l2,jlat1,jlat2,jlon1,jlon2
                            raise
                    else:
                        #print jlat1,jlat2,jlon1,jlon2
                        jlat1=-200
                        jlat2=-200
                        xymap[j,i,:] = zmap[jlat,jlon,:]
                else:
                    #print "not bilinear"
                    xymap[j,i,:] = zmap[jlat,jlon,:]
                    
    return xymap


def orthographic_single(lon,lat,imap,l0=0,p0=0,ny=36,nx=36,interp='bilinear'):
    xymap = np.zeros((ny,nx))
    xymap[:] = 0.0
    coords = np.zeros((ny,nx,2))
    coords2 = np.zeros((ny,nx,2,3))
    coords3 = np.zeros((ny,nx,2))
    rad=0.5*(8*nx/18.0+8*ny/18.0)

    zlat = np.zeros(np.array(lat.shape)+2)
    zlon = np.zeros(np.array(lon.shape)+2)
    dlat = np.diff(lat)
    dlon = np.diff(lon)
    zlon[0] = lon[0]-dlon[0]
    zlon[1:-1] = lon[:]
    zlon[-1] = lon[-1]+dlon[-1]
    zlat[0] = lat[0]-dlat[0]
    zlat[-1] = lat[-1]+dlat[-1]
    zlat[1:-1] = lat[:]
    zmap = np.zeros(np.array(imap.shape)+np.array((2,2)))
    zmap[1:-1,1:-1] = imap[:,:]
    zmap[0,:] = zmap[1,:]
    zmap[-1,:] = zmap[-2,:]
    zmap[:,0] = zmap[:,-2]
    zmap[:,-1] = zmap[:,1]
    
    if l0>180.0:
        l0 -= 360.0

    if zlon.max()>180.0:
        zlon[zlon>180]-=360.0    

    p0 *= np.pi/180.0
    l0 *= np.pi/180.0
    xx = np.arange(nx)-nx/2
    yy = np.arange(ny)-ny/2
    for j in range(0,ny):
        jy = yy[j]
        for i in range(0,nx):
            ix = xx[i]
            if (ix**2+jy**2)<=rad**2:
                rho = np.sqrt(ix**2+jy**2)
                cc = np.arcsin(rho/rad)
                if rho==0:
                    phi = p0
                    lamb = l0
                else:
                    phi = np.arcsin(np.cos(cc)*np.sin(p0) + jy*np.sin(cc)*np.cos(p0)/rho)
                    lamb = l0 + np.arctan2(ix*np.sin(cc),(rho*np.cos(cc)*np.cos(p0)-jy*np.sin(cc)*np.sin(p0)))
                phi *= 180.0/np.pi
                lamb *= 180.0/np.pi
                if lamb<-180:
                    lamb += 360
                if lamb>180:
                    lamb -= 360
                jlat = np.argmin(abs(phi-zlat))
                jlon = np.argmin(abs(lamb-zlon))
                if interp=="bilinear":
                    xlat1 = np.where(zlat<phi)[0]
                    xlat2 = np.where(zlat>=phi)[0]
                    xlon1 = np.where(zlon<lamb)[0]
                    xlon2 = np.where(zlon>=lamb)[0]
                    
                    if len(xlat1)>0 and len(xlat2)>0 and len(xlon1)>0 and len(xlon2)>0:
                        #jlat1 = jlat1[0]
                        #jlat2 = jlat2[-1]
                        #jlon1 = jlon1[-1]
                        #jlon2 = jlon2[0]
                        try:
                            jlat1 = np.where(zlat==zlat[xlat1].max())[0][0]
                            jlat2 = np.where(zlat==zlat[xlat2].min())[0][0]
                            jlon1 = np.where(zlon==zlon[xlon1].max())[0][0]
                            jlon2 = np.where(zlon==zlon[xlon2].min())[0][0]
                        except:
                            print xlat1
                            print xlat2
                            print xlon1
                            print xlon2
                            raise
                        p1 = zlat[jlat1]
                        p2 = zlat[jlat2]
                        l1 = zlon[jlon1]
                        l2 = zlon[jlon2]
                        dl = l2-l1
                        dp = p2-p1
                        #print "interpolating!"
                        try:
                            if dl==0 and dp>0:
                                fxy1 = zmap[jlat1,jlon]
                                fxy2 = zmap[jlat2,jlon]
                                xymap[j,i] = (p2-phi)/dp*fxy1 + (phi-p1)/dp*fxy2
                            elif dp==0 and dl>0:
                                fx1y = zmap[jlat,jlon1]
                                fx2y = zmap[jlat,jlon2]
                                xymap[j,i] = (l2-lamb)/dl*fx1y + (lamb-l1)/dl*fx2y
                            elif dp==0 and dl==0:
                                xymap[j,i] = zmap[jlat,jlon]
                            else:
                                fxy1 = (l2-lamb)/dl*zmap[jlat1,jlon1] + (lamb-l1)/dl*zmap[jlat1,jlon2]
                                fxy2 = (l2-lamb)/dl*zmap[jlat2,jlon1] + (lamb-l1)/dl*zmap[jlat2,jlon2]
                                xymap[j,i] = (p2-phi)/dp*fxy1 + (phi-p1)/dp*fxy2
                        except:
                            print p1,p2,l1,l2,jlat1,jlat2,jlon1,jlon2
                            raise
                    else:
                        #print jlat1,jlat2,jlon1,jlon2
                        jlat1=-200
                        jlat2=-200
                        xymap[j,i] = zmap[jlat,jlon]
                else:
                    #print "not using bilinear"
                    xymap[j,i] = zmap[jlat,jlon]
                    
    return xymap

def getphase(phasecurve,nphase,phase,ntime):
    ln = phasecurve.variables['lon'][:]
    lt = phasecurve.variables['lat'][:]
    color = phasecurve.variables['colors'][ntime,nphase,:,:,:]
    color /= 5.0*np.mean(color)
    atimes = np.array([0.,90.,180.,270.])
    if phase=='N':
        p0 = 90.0
    elif phase=='S':
        p0 = -90.0
    elif phase=="E":
        atimes += 90.0
        p0 = 0.0
    elif phase=="W":
        atimes += -90.0
        p0 = 0.0
    else:
        p0 = 0.0
    atimes[atimes>180] -= 360.0
    proj= orthographic(ln,lt,np.minimum(color,1.0),l0=atimes[0],p0=p0,nx=200,ny=200)
    return proj

if __name__=="__main__":
    filename = sys.argv[1]
    times = sys.argv[2]
    times = times.split('^')
    if times[-1]=='':
        times = times[:-1]
    times = np.array(times).astype(int)
    pc = nc.Dataset(filename,"r")
    tag = filename.split("_phases.nc")[0]
    os.system("mkdir "+tag)
    phases = pc.variables['phase'][:]
    for ntime in times:
        for p in range(0,len(phases)):
            proj = getphase(pc,p,phases[p],ntime)
            f,a=plt.subplots(figsize=(14,12))
            plt.imshow(proj,interpolation='gaussian',origin='lower')
            plt.xticks([])
            plt.yticks([])
            plt.savefig(tag+"/"+tag+"_%03d_%s.png"%(ntime,phases[p]),bbox_inches='tight')
        plt.close('all')
