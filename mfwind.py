import pylab as plt
import numpy as np
import sys
import os
import caesar
import function as fu
import OBSSMF as obs

MODEL = sys.argv[1]
WIND = sys.argv[2:]
TYPES = ['GSMF']
path = '/home/cossim/HELUCID/'
nrowmax = 3
mlim = 8.7

if "m12.5n128" in MODEL:
    SNAPNUM = [19,27,32,41,46,54,64,77]
elif "nifty" in MODEL:
    SNAPNUM = [128]
else:
    #SNAPNUM = [46,62,70,85,105,135]
    SNAPNUM = [11,14,15,17,22, 45, 76, 153]
    #SNAPNUM = [36,51,62,78,105,151]
#     SNAPNUM = np.arange(1,9)
    #SNAPNUM = [78,105,151]

colors = ('b', 'crimson', 'g', 'c', 'm', 'y', 'k')
# modeldict = {'fh_qr':'Mufasa','s48':'Simba','s43nojet':'No-jet'}

def haloMF(obj,mbin,ax):
    theory = MFtheory(obj.simulation.redshift,obj.simulation.hubble, obj.snapAttribs.O0, obj.snapAttribs.Ol)
    THEORY=limitTheory(mbin,theory)
    ax.plot(np.log10(THEORY['M']),np.log10(THEORY['ST']),'k--',label='ST')

def gsmf_Illustris(zbin,ax):
    zbin = np.round(zbin,0)
    infile = 'Observations/Illustris/data/stellar_mass_function_z%g.txt' % (zbin)
    if os.path.isfile(infile):
        m,phi = np.loadtxt(infile,unpack=True)
        ax.plot(np.log10(m),np.log10(phi),'--',label='Illustris',color='r')
    else: print('illustris: file not found %s',infile)

def gsmf_eagle(zbin,ax):
    zbin = np.round(zbin,0)
    infile = 'Observations/eagle/gsmf/Ref_gsmf_z%gp0.txt' % (zbin)
    if zbin == 0:
        infile = 'Observations/eagle/gsmf/Ref_gsmf_z%gp1.txt' % (zbin)
    if os.path.isfile(infile):
        m,phi = np.loadtxt(infile,usecols=(0,2),unpack=True)
        ax.plot(np.log10(m),np.log10(phi),':',label='EAGLE',color='c')
    else: print('eagle: file not found %s',infile)

def gsmf_MBII(zbin,ax):
    infile = 'Observations/MBII/massfunction.txt'
    if os.path.isfile(infile):
        m,phi_z4,phi_z3,phi_z2,phi_z1,phi_z0,phi_Ill = np.loadtxt(infile,usecols=(0,1,2,3,4,5,7),unpack=True)
        if zbin >= 4: phi = phi_z4
        elif zbin == 3: phi = phi_z3
        elif zbin == 2: phi = phi_z2
        elif zbin == 1: phi = phi_z1
        elif zbin == 0: phi = phi_z0
        elif zbin == -1: phi = phi_Ill
        ax.plot(np.log10(m),np.log10(phi),'--',label='MB-II',color='c')


def massFunc(objs,labels,ax,jwind):
    for j in range(0,len(objs)):
        for curType in TYPES:
            galpos = np.array([g.pos.d for g in objs[j].galaxies])
            if len(galpos) == 0: continue
            if curType == 'GSMF': mass = np.array([g.masses['stellar'].d for g in objs[j].galaxies])
            elif curType == 'HI': mass = np.array([g.masses['HI'].d for g in objs[j].galaxies]) 
            elif curType == 'H2': mass = np.array([g.masses['H2'].d for g in objs[j].galaxies]) 
            elif curType == 'SFR': mass = np.array([g.sfr.d for g in objs[j].galaxies]) 
            elif curType == 'Halo': 
                mass = np.array([h.masses['virial'].d for h in objs[j].halos]) 
                galpos = np.array([h.pos.d for h in objs[j].halos]) 
            volume = objs[j].simulation.boxsize.to('Mpccm').d**3 * objs[j].simulation.ndm/5400**3 #need this fraction to correct the zoomin simulations!!
            mlim = np.log10(32*objs[j].simulation.critical_density*objs[j].simulation.boxsize**3*objs[j].simulation.omega_baryon/objs[j].simulation.effective_resolution**3*objs[j].simulation.scale_factor**3/objs[j].simulation.hubble_constant**3) # galaxy mass resolution limit: 32 gas particle masses

            x,y,sig = fu.cosmic_variance(mass, galpos, objs[j].simulation.boxsize, volume, nbin=20, minmass=mlim-1)
            xmf,ymf, bins, nsp= fu.mass_function(mass, volume, nbin=20, minmass=mlim-1)
            ncol = int((len(objs)-1)/nrowmax)+1
            icol = int(j/nrowmax)
            irow = j%nrowmax
            print('z=',objs[j].simulation.redshift,labels[jwind],'mlim=',mlim,'j=',j,'irow=',irow,'icol=',icol,'ncol=',ncol)
            if ncol==1: ax0 = ax[irow]
            else: ax0 = ax[irow][icol]
            if jwind==0: ltype = '-'
            else: ltype = '--'
            ax0.plot(np.log10(x)+jwind*0.001,np.log10(y),ltype,color=colors[jwind],label=labels[jwind])
            ax0.plot(np.log10(xmf),np.log10(ymf),'r--',label=labels[jwind])
            print(xmf,ymf)
            elo = sig
            ehi = sig
            if jwind >= 0:
                ax0.errorbar(np.log10(x)+j*0.001,np.log10(y),yerr=[elo,ehi],color=colors[jwind])

            if jwind==len(WIND)-1: 
                smf = obs.SMF(objs[j].simulation.redshift)
                if smf.cond:
                    ax0.errorbar(smf.x, smf.y, yerr=smf.yerr, fmt='x', label=smf.name, zorder=100)
                    ax0.legend(loc='upper right', fontsize=8).draw_frame(False)
                gsmf_eagle(objs[j].simulation.redshift,ax0)
                #gsmf_Illustris(objs[j].simulation.redshift,ax[j])
                #gsmf_MBII(z)
                ax0.legend(loc='lower left',fontsize=8)
                ax0.annotate('z=%g'%np.round(objs[j].simulation.redshift,1), xy=(0.8, 0.75), xycoords='axes fraction',size=12,bbox=dict(boxstyle="round", fc="w"))
                ax0.set_ylabel(r'$\log \Phi [Mpc^{-3}$]',fontsize=16)
                if irow == nrowmax-1: ax0.set_xlabel(r'$\log M_* [M_\odot]$',fontsize=16)
    return mlim

###########################################################################

if __name__ == '__main__':

    for j in range(0,len(WIND)):
        sims = []
        labels = []
        for k in range(0,len(SNAPNUM)):
            caesarfile = path+'%s/%s/Caesar_snap_%03d.hdf5' % (MODEL,WIND[j],SNAPNUM[k])
            if os.path.isfile(caesarfile):
                sims.append(caesar.load(caesarfile))
                #labels.append('Simba '+WIND[j])
                try:
                    labels.append(modeldict[WIND[j]])
                except:
                    labels.append(WIND[j])
            else:
                print('Could not find caesar file %s'%caesarfile)
                continue
        ncol = int((len(sims)-1)/nrowmax)+1
        nrow = min(len(sims),nrowmax)
        if j==0: 
            fig,ax = plt.subplots(ncols=ncol,nrows=nrow,sharey=True,figsize=(3+3*ncol,2+2*nrow))
        print(len(sims),'ncol=',ncol,'nrow=',nrow)
        if ncol==nrow==1:
            mlim = massFunc(sims,labels,[ax],j)
        else:
            mlim = massFunc(sims,labels,ax,j)
        if len(sims) == 0: sys.exit('No caesar files found!')

    plt.minorticks_on()
    plt.xlabel(r'$\log M_* [M_\odot]$',fontsize=16)
    plt.xlim(mlim,)
    plt.ylim(-6.1,-1.1)
    plt.subplots_adjust(hspace=.0)

    figname = 'test_plot/mfwind_%s.pdf' % (MODEL)
    plt.savefig(figname,bbox_inches='tight')
    plt.show()

