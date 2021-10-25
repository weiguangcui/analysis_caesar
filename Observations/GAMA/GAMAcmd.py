
import pylab as plt
import numpy as np
from astropy.io import fits
import sys

GAMAFILE = 'GAMAStellarMasses.dat'
logMmin = 9
color = 'uminusr'

colors = ('r', 'g', 'b', 'c', 'm', 'y', 'k')

def cmdData():  # Schawinski+14 Fig 2
    ax = plt.axis()
    msrange = np.arange(ax[0],ax[1],0.1)
    gvhi = -0.25+0.25*msrange
    gvlo = -0.5+0.25*msrange
    plt.fill_between(msrange,gvlo,gvhi,color='g',alpha=0.3,label='Schawinki+14 Green Valley')
    #plt.plot(msrange,gvhi,'--',lw=2,c='k',label='Schawinki+14 Green Valley')
    #plt.plot(msrange,gvlo,'--',lw=2,c='k')

def GAMA_printfields():
    PARFILE = 'GAMAStellarMasses.par'
    name = np.loadtxt(PARFILE,usecols=(0,),dtype=(np.str),unpack=True)
    print name

def GAMA_loaddata(field):
    PARFILE = 'GAMAStellarMasses.par'    
    GAMAFILE = 'GAMAStellarMasses.dat'
    name,col = np.loadtxt(PARFILE,usecols=(0,1),dtype=(np.str),unpack=True)
    colnum = col[name==field][0]
    colnum = np.int(colnum)-1
    data = np.loadtxt(GAMAFILE,usecols=(colnum,))
    print field,colnum,len(data)
    return data

def GAMA_contour():
    zlim = 0.15
    mstar = GAMA_loaddata('logmstar')
    fluxscale = GAMA_loaddata('fluxscale')
    mstar = mstar + np.log10(fluxscale) - 2*np.log10(0.7/0.7)  # correct for aperture & cosmology
    colorvec = GAMA_loaddata(color)[mstar>logMmin]
    z_GAMA = GAMA_loaddata('Z')[mstar>logMmin]
    fluxscale = fluxscale[mstar>logMmin]
    mstar = mstar[mstar>logMmin]
    ax = plt.axis()
    mstar = mstar[z_GAMA<zlim]
    colorvec = colorvec[z_GAMA<zlim]
    print 'GAMA: Ndata=',len(mstar),' z= ',np.average(z_GAMA[z_GAMA<zlim])
    #plt.plot(mstar,colorvec,'o',ms=1)
    plt.hexbin(mstar,colorvec,cmap=plt.cm.BuPu,label='GAMA, z<%g' % (zlim))
    counts,xedges,yedges = np.histogram2d(mstar,colorvec,bins=(50,50),range=[[ax[0],ax[1]],[ax[2],ax[3]]])
    plt.contour(counts.transpose(),extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],linewidths=(1,1.5,2,2.5,3),colors=colors)


logMmin = 9
plt.xlim(logMmin-0.1,12)
plt.ylim(0.3,3.3)

GAMA_printfields()
GAMA_contour()

cmdData()

plt.ylabel(r'$u-r$',fontsize=20)
plt.xlabel(r'$\log M_*$' ,fontsize=20)
plt.legend(loc='lower right',fontsize=18)

plt.show()

