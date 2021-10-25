
import caesar
import sys
import pylab as plt
import numpy as np
import plotmedian as pm

# define input file
MODEL = sys.argv[1]
SNAP = int(sys.argv[2])
WIND = sys.argv[3:]
mmin = 5.8e8
SOLAR = 0.0134

colors = ('b', 'crimson', 'g', 'c', 'm', 'y', 'k')
# modeldict = {'fh_qr':'Mufasa','s49':'Simba','s43nojet':'No-jet','s48':'Old Simba','s49nox':'No X-ray','s50':'Simba','s50fedd':'Simba','s50j7k':'Simba'}
modeldict = {}
for i in WIND:
    modeldict[i]=i
    
def plot_data(z,ltype):
    zfact = z/(1+z)
    M1 = 11.59+1.195*zfact
    N = 0.0351-0.0247*zfact
    beta = 1.376-0.826*zfact
    gamma = 0.608+0.329*zfact
    mhalo = np.arange(10.5,14,0.1)
    mratio = 10**mhalo/10**M1
    msmh_moster = 2*N / (mratio**-beta + mratio**gamma)
    #print mh,msmh_moster
    plt.plot(mhalo,np.log10(msmh_moster),ltype,lw=3,color='r',label='Moster+13 z=%g'%(0.01*int(z*100)))

def plot_peeples14(ax):	# Peeples+14 ApJ 786 54: f_coldgas vs M*
    logMs = np.asarray([6.7,7.1,7.6,8.2,8.6,9.1,9.6,10.1,10.6,11,11.4])
    fgas = np.asarray([7.4,3.7,5.6,5.8,3.6,1.3,0.6,0.5,0.23,0.17,0.067])
    fghi = np.asarray([12,19.5,15.7,8.4,4.3,2.2,1.7,1.2,0.52,0.34,0.12])
    fglo = np.asarray([4.6,2.1,2.2,4.8,2.6,1.0,0.32,0.083,0.056,0.082,0.057])
    logfg = np.log10(fgas)
    elfghi = np.log10(fghi)-logfg
    elfglo = logfg-np.log10(fglo)
    ax.errorbar(logMs,logfg,yerr=[elfglo,elfghi],fmt='x',ms=6,color='k',lw=2,label='P14')

def plot_saintonge17(ax):	# Saintonge+17: f_H2 vs M* from Table 5
    # Main sequence only
    logMs = np.asarray([9.406,9.643,9.854,10.07,10.25,10.46,10.68,10.82,11.10,11.31])
    logfH2 = np.asarray([-0.93,-0.89,-0.96,-0.89,-1.05,-0.99,-1.10,-1.26,-1.41,-1.68])
    elfH2 = 3*np.asarray([0.049,0.033,0.023,0.017,0.023,0.018,0.014,0.013,0.010,0.039])
    #ax.errorbar(logMs,logfH2,yerr=[elfH2,elfH2],fmt='+',ms=6,color='cyan',lw=2,label='S17-SF')
    # All
    logMs = np.asarray([9.407,9.648,9.848,10.04,10.24,10.46,10.67,10.87,11.07,11.27])
    logfH2 = np.asarray([-1.01,-1.08,-0.97,-0.90,-1.08,-1.05,-1.34,-1.41,-1.66,-2.02])
    elfH2 = 3*np.asarray([0.041,0.035,0.020,0.011,0.013,0.009,0.011,0.009,0.012,0.041])
    ax.errorbar(logMs,logfH2,yerr=[elfH2,elfH2],fmt='s',ms=6,fillstyle='none',color='k',lw=3,label='Saintonge+17')

def plot_catinella12(ax):	# Catinella+12: <f_HI> vs M* from Table 1, assuming lower limits=0
    logMs = np.asarray([10.16,10.46,10.77,11.06,11.30])
    fHI = np.asarray([0.262,0.136,0.078,0.038,0.019])
    efHI = 3*np.asarray([0.035,0.019,0.009,0.007,0.006])
    logfHI = np.log10(fHI)
    ehi = np.log10(fHI+efHI)-logfHI
    elo = logfHI-np.log10(fHI-efHI)
    ax.errorbar(logMs,logfHI,yerr=[elo,ehi],fmt='s',ms=6,fillstyle='none',color='k',lw=3,label='Catinella+12')

#fig,(ax1,ax2,ax3) = plt.subplots(3, sharex=True, sharey=True, figsize=(6,6))
fig,(ax2,ax3) = plt.subplots(2, sharex=True, sharey=True, figsize=(6,6))
for iwind in range(0,len(WIND)):
# load in input file
    infile = '/home/cossim/HELUCID/%s/%s/Caesar_snap_%03d.hdf5' % (MODEL,WIND[iwind],SNAP)
    sim = caesar.load(infile)
    redshift = sim.simulation.redshift

    icent = np.asarray([i.central for i in sim.galaxies])
    pos = np.asarray([i.pos for i in sim.galaxies])
    ids = np.asarray([i.GroupID for i in sim.galaxies if i.central == 1])
    ms = np.asarray([i.masses['stellar'] for i in sim.galaxies if i.central == 1])
    mHI = np.asarray([i.masses['HI'] for i in sim.galaxies if i.central == 1])
    mH2 = np.asarray([i.masses['H2'] for i in sim.galaxies if i.central == 1])
    sfr = np.asarray([i.sfr for i in sim.galaxies if i.central == 1])
    met = np.asarray([i.metallicities['sfr_weighted'] for i in sim.galaxies if i.central == 1])
    #mh = np.asarray([i.halo.masses['total'] for i in sim.galaxies if i.central == 1])
    ms_sat = np.asarray([i.masses['stellar'] for i in sim.galaxies if i.central != 1])
    mHI_sat = np.asarray([i.masses['HI'] for i in sim.galaxies if i.central != 1])
    mH2_sat = np.asarray([i.masses['H2'] for i in sim.galaxies if i.central != 1])
    sfr_sat = np.asarray([i.sfr for i in sim.galaxies if i.central != 1])
    met_sat = np.asarray([i.metallicities['sfr_weighted'] for i in sim.galaxies if i.central != 1])
    #central_galaxy_halo_masses = [i.halo.masses['total'] for i in sim.galaxies if i.masses['total'] > 1.0e11]

    ssfr = np.log10(1.e9*sfr/ms+10**(-2.5+0.3*redshift))
    ssfr_sat = np.log10(1.e9*sfr_sat/ms_sat+10**(-2.5+0.3*redshift))
    fgas = np.log10(mHI/ms+mH2/ms+1.e-3)
    fgas_sat = np.log10(mHI_sat/ms_sat+mH2_sat/ms_sat+1.e-3)
    fHI = np.log10(mHI/ms+1.e-3)
    fHI_sat = np.log10(mHI_sat/ms_sat+1.e-3)
    fH2 = np.log10(mH2/ms+1.e-3)
    fH2_sat = np.log10(mH2_sat/ms_sat+1.e-3)
    Rmol = np.log10(mH2/(mHI+1.e-3) + 1.e-3)
    Rmol_sat = np.log10(mH2_sat/(mHI_sat+1.e-3) + 1.e-3)
    yflg = ((np.append(ssfr,ssfr_sat)>-2))
    #fig,(ax1,ax2,ax3,ax4) = plt.subplots(4, sharex=True, sharey=True)

    colorcode = r'$\Delta log Z/Z_\odot $'
    colorcode = r'$\Delta$sSFR'
    if colorcode == r'$\Delta$sSFR':
        cvec = ssfr
        cvec_sat = ssfr_sat
        massbin,cvecbin,ebinlo,ebinhi = pm.runningmedian(np.log10(ms[sfr>0]),ssfr[sfr>0])
    elif colorcode == r'$\Delta log Z/Z_\odot $':
        cvec = np.log10(met/SOLAR)
        cvec_sat = np.log10(met_sat/SOLAR)
        massbin,cvecbin,ebinlo,ebinhi = pm.runningmedian(np.log10(ms[met>0]),np.log10(met[met>0]/SOLAR),bins=20)
    cvec = cvec - np.interp(np.log10(ms),massbin,cvecbin)
    cvec_sat = cvec_sat - np.interp(np.log10(ms_sat),massbin,cvecbin)
    colormap = plt.cm.jet_r
    
    '''
    yflg = ((np.append(mH2,mH2_sat)>1.e6))
    if iwind==0:
        im1 = ax1.scatter(np.log10(ms),fgas, c=cvec, s=8, lw=0, cmap=colormap)
        im1 = ax1.scatter(np.log10(ms_sat),fgas_sat, c=cvec_sat, s=2, lw=0, cmap=colormap)
        fig.colorbar(im1,ax=ax1,label=colorcode)
        plot_peeples14(ax1)
        ax1.set_ylabel(r'$\log\ f_{gas}$',fontsize=16)
    pm.plotmedian(np.append(np.log10(ms),np.log10(ms_sat)),np.append(fgas,fgas_sat),yflag=yflg,c=colors[iwind],lw=2-min(iwind,1),ax=ax1,bins=20,label=modeldict[WIND[iwind]],stat='mean')
    ax1.legend(loc='lower left')
    '''
    
    yflg = ((np.append(mH2,mH2_sat)>1.e6))
    if iwind==0:
        im2 = ax2.scatter(np.log10(ms),fH2, c=cvec, s=8, lw=0, cmap=colormap)
        fig.colorbar(im2,ax=ax2,label=colorcode)
        im2 = ax2.scatter(np.log10(ms_sat),fH2_sat, c=cvec_sat, s=2, lw=0, cmap=colormap)
        ax2.scatter(np.log10(ms[sfr>100]),fH2[sfr>100], marker='x', c='k', s=20, label='Simba SFR>100')
        if redshift<0.5: plot_saintonge17(ax2)
        ax2.set_ylabel(r'$\log\ f_{H2}$',fontsize=16)
    print(np.append(np.log10(ms),np.log10(ms_sat)),np.append(fH2,fH2_sat))
    pm.plotmedian(np.append(np.log10(ms),np.log10(ms_sat)),np.append(fH2,fH2_sat),yflag=yflg,c=colors[iwind],lw=2-min(iwind,1),ax=ax2,bins=16,label=modeldict[WIND[iwind]],stat='mean') #,pos=pos,boxsize=sim.simulation.boxsize* sim.simulation.ndm/5400**3)
    #yflg = ((np.append(mH2,mH2_sat)>1.e6)&(np.append(ssfr,ssfr_sat)>-2))
    #pm.plotmedian(np.append(np.log10(ms),np.log10(ms_sat)),np.append(fH2,fH2_sat),yflag=yflg,c='b',lw=2-min(iwind,1),ax=ax2,bins=20,label=modeldict[WIND[iwind]],stat='mean')
    ax2.legend(loc='lower left')
    
    yflg = ((np.append(fHI,fHI_sat)>-2)&(np.append(mHI,mHI_sat)>5.e8))
    if iwind==0:
        im3 = ax3.scatter(np.log10(ms),fHI, c=cvec, s=8, lw=0, cmap=colormap)
        fig.colorbar(im3,ax=ax3,label=colorcode)
        im3 = ax3.scatter(np.log10(ms_sat),fHI_sat, c=cvec_sat, s=2, lw=0, cmap=colormap)
        ax3.scatter(np.log10(ms[sfr>100]),fHI[sfr>100], marker='x', c='k', s=20, label='Simba SFR>100')
        if redshift<0.5: plot_catinella12(ax3)
        ax3.set_ylabel(r'$\log\ f_{HI}$',fontsize=16)
    pm.plotmedian(np.append(np.log10(ms),np.log10(ms_sat)),np.append(fHI,fHI_sat),yflag=yflg,c=colors[iwind],lw=2-min(iwind,1),ax=ax3,bins=16,label=modeldict[WIND[iwind]],stat='mean') #,pos=pos,boxsize=sim.simulation.boxsize* sim.simulation.ndm/5400**3)
    ax3.legend(loc='lower left')
    
    #im4 = ax4.scatter(np.log10(ms),Rmol, c=cvec, s=5, lw=0, cmap=colormap)
    #im4 = ax4.scatter(np.log10(ms_sat),Rmol_sat, c=cvec_sat, s=2, lw=0, cmap=colormap)
    #pm.plotmedian(np.append(np.log10(ms),np.log10(ms_sat)),np.append(Rmol,Rmol_sat),yflag=yflg,c='g',ax=ax4)
    #ax4.set_ylabel(r'$R_{\rm mol}$',fontsize=16)
    #ax4.set_ylim(-3.1,2.0)
    #fig.colorbar(im4,ax=ax4,label='sSFR')
    
plt.subplots_adjust(hspace=.0)
#plt.rc('text', usetex=True)
plt.minorticks_on()
plt.ylim(-3.4+0.5*redshift,1.0)
plt.xlim(np.log10(mmin),12.3)
plt.xlabel(r'$\log\ M_{*}$' ,fontsize=16)

plt.annotate('z=%g'%(np.round(redshift,0)), xy=(0.8, 0.9), xycoords='axes fraction',size=14,bbox=dict(boxstyle="round", fc="w"))

plt.savefig('test_plot/fgas_%s.pdf'%MODEL, bbox_inches='tight', format='pdf')

plt.show()

