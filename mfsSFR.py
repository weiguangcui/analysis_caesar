import pylab as plt
import numpy as np
import sys
import os
import caesar
import function as fu


SNAPNUM = int(sys.argv[1])
WIND = sys.argv[2:]

path = '/home/cossim/HELUCID/'
# mlim = 8.7

colors = ('b', 'crimson', 'g', 'c', 'm', 'y', 'k')
# modeldict = {'fh_qr':'Mufasa','s48':'Simba','s43nojet':'No-jet'}

def cal_mass_function(mass, vol_Mpc, nbin=30, minmass=1):
    # everything in log
    #not very useful, used in cosmic_variance()
    maxmass     = np.nanmax(mass)
    lbin        = np.linspace(minmass, maxmass, nbin+1)
    bin_        = 10** lbin
    step        = ( maxmass - minmass ) / nbin
    x           = 10 ** ( (lbin + step/2.)[:nbin])
    hist        = np.histogram( 10**mass, bins=bin_, range=(bin_.min(),bin_.max()) )[0]
    y           =  hist / (vol_Mpc * step)
    return x, y, bin_, step


def OSFRF(fname,ax):
    infile = 'Observations/%s' % (fname)
    if os.path.isfile(infile):
        m,phi,err = np.loadtxt(infile,unpack=True)
        ax.errorbar(m,phi,yerr=err,ls='--',label='Zhao+20',color='k',lw=2,marker='P',ms=8,zorder=-99,alpha=0.6)
    else: print('File not found %s',infile)

def massFunc(objs,labels,ax,jwind):
    galpos = np.array([g.pos.d for g in objs[jwind].galaxies])
    if len(galpos) == 0: return 0
    gmasses=np.array([g.masses['stellar'].d for g in objs[jwind].galaxies]) 
    gsfr = np.array([g.sfr.d for g in objs[jwind].galaxies])
    idm=(gmasses>10**9.)&(gsfr>10**-1.5)
    volume = objs[jwind].simulation.boxsize.to('Mpccm').d**3 * objs[jwind].simulation.ndm/5400**3
    for j in range(2):
        if j==0: 
            curType = 'SFR' 
            mass = np.log10(gsfr[idm])
            ofn='SFRF-Zhao2020.txt'
        elif j==1: 
            curType = 'sSFR' #: mass = np.array([g.sfr.d/g.masses['stellar'].d for g in objs[j].galaxies]) 
            mass = np.log10(gsfr[idm]/gmasses[idm])
            ofn='sSFRF-Zhao2020.txt'
        
#         mlim = np.log10(32*objs[j].simulation.critical_density*objs[j].simulation.boxsize**3*objs[j].simulation.omega_baryon/objs[j].simulation.effective_resolution**3*objs[j].simulation.scale_factor**3/objs[j].simulation.hubble_constant**3) # galaxy mass resolution limit: 32 gas particle masses

#         x,y,sig = fu.cosmic_variance(mass, galpos, objs[j].simulation.boxsize, volume, nbin=20, minmass=mass.min())
        x,y,bins,step=cal_mass_function(mass,volume,nbin=10, minmass=mass.min())
#         print(x,y,mass.min(),mass.max())
        ncol = 1
        icol = 0
        irow = j
#         print('z=',objs[0].simulation.redshift,labels[jwind],'mlim=',mass.min(),'j=',j,'irow=',irow,'icol=',icol,'ncol=',ncol)
        if ncol==1: ax0 = ax[irow]
        else: ax0 = ax[irow][icol]
        if jwind==0: ltype = '-'
        else: ltype = '--'
        ax0.plot(np.log10(x),y,ltype,color=colors[jwind],label=labels[jwind])
#         elo = sig
#         ehi = sig
#         if jwind >= 0:
#             ax0.errorbar(np.log10(x)+j*0.001,np.log10(y),yerr=[elo,ehi],color=colors[jwind])

        if jwind==len(WIND)-1: 
#             smf = obs.SMF(objs[j].simulation.redshift)
#             if smf.cond:
#                 ax0.errorbar(smf.x, smf.y, yerr=smf.yerr, fmt='x', label=smf.name, zorder=100)
#                 ax0.legend(loc='upper right', fontsize=8).draw_frame(False)
#             gsmf_eagle(objs[j].simulation.redshift,ax0)
            #gsmf_Illustris(objs[j].simulation.redshift,ax[j])
            #gsmf_MBII(z)
            OSFRF(ofn,ax0)
            ax0.set_yscale('log')           
            ax0.legend(loc='lower left',fontsize=8)
            ax0.annotate('z=%g'%np.round(objs[jwind].simulation.redshift,1), xy=(0.8, 0.75), xycoords='axes fraction',size=12,bbox=dict(boxstyle="round", fc="w"))
#             ax0.set_ylabel(r'$\log \Phi [Mpc^{-3}$]',fontsize=16)
#             if irow == nrowmax-1: ax0.set_xlabel(r'$\log M_* [M_\odot]$',fontsize=16)
    return mass.min()

###########################################################################

if __name__ == '__main__':
    sims = []
    labels = []
    for j in range(0,len(WIND)):
#         for k in range(0,len(SNAPNUM)):
        caesarfile = path+'z0.025/%s/Caesar_snap_%03d.hdf5' % (WIND[j],SNAPNUM)

        if os.path.isfile(caesarfile):
            sims.append(caesar.load(caesarfile))
            if sims[-1].simulation.redshift > 0.15:
                print('Warrning the simulation redshift %d may different very much from observation at z=0.1')
            #labels.append('Simba '+WIND[j])
            try:
                labels.append(modeldict[WIND[j]])
            except:
                labels.append(WIND[j])
        else:
            print('Could not find caesar file %s'%caesarfile)
            continue
        ncol = 2
        nrow = 1
        if j==0: 
            fig,ax = plt.subplots(ncols=ncol,nrows=nrow,sharey=True,figsize=(3+3*ncol,2+2*nrow))
            print(len(sims),'ncol=',ncol,'nrow=',nrow)
        mlim = massFunc(sims,labels,ax,j)
        if len(sims) == 0: sys.exit('No caesar files found!')

    plt.minorticks_on()
    ax[0].set_xlabel(r'$\log\ SFR\ [M_\odot\ yr^{-1}]$',fontsize=16)
    ax[1].set_xlabel(r'$\log\ sSFR\ [yr^{-1}]$',fontsize=16)
    ax[0].set_ylabel(r'$\log\ \Phi\ [Mpc^{-3}$]',fontsize=16)
#     plt.xlim(mlim,)
#     plt.ylim(-6.1,-1.1)
    plt.subplots_adjust(hspace=.0,wspace=0.)
    ax[0].grid()
    plt.grid()
    
    figname = 'test_plot/mf_sSFRF.pdf' #% (MODEL)
    plt.savefig(figname,bbox_inches='tight')
    plt.show()

