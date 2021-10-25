import caesar
from readgadget import *
import sys
import pylab as plt
import numpy as np

if len(sys.argv) < 3:
    print('usage: MODEL SNAP WIND1 WIND2 ...')
    exit()

# define input file
MODEL = sys.argv[1]
SNAP = int(sys.argv[2])
WIND = sys.argv[3:]
#plt.rc('text', usetex=True)
mlim = 8.7
mmax = 12.5

def plot_data(redshift):
    if redshift < 0.5:
        infile = 'Observations/KormendyHo2013/KH13.dat'
        hubtype,mbh,mbhlo,mbhhi,sig,esig = np.loadtxt(infile,usecols=(2,11,12,13,14,15),unpack=True)
        #print mbh,mbhlo,mbhhi,sig,esig
        embh = [np.log10(mbh)-np.log10(mbh-mbhlo),np.log10(mbhhi+mbh)-np.log10(mbh)]
        elogsig = [np.log10(sig)-np.log10(sig-esig),np.log10(sig+esig)-np.log10(sig)]
        plt.errorbar(np.log10(sig),np.log10(mbh*1.e6),fmt='.',yerr=embh,xerr=elogsig,lw=1,color='grey',label='Kormendy+Ho')

# load in input file
fig,ax = plt.subplots()
for iwind in range(0,len(WIND)):
#     infile = '/home/wgcui/%s/%s/CaesarGroups_%03d.hdf5' % (MODEL,WIND[iwind],SNAP)
    infile = '/home/cossim/HELUCID/%s/%s/Caesar_snap_%03d.hdf5' % (MODEL,WIND[iwind],SNAP)
    sim = caesar.load(infile)
    redshift = sim.simulation.redshift
    
    ids = np.asarray([i.GroupID for i in sim.galaxies if i.central == 1])
    ms = np.asarray([i.masses['stellar'] for i in sim.galaxies if i.central == 1])
    mbh = np.asarray([i.masses['bh'] for i in sim.galaxies if i.central == 1])
    sfr = np.asarray([i.sfr for i in sim.galaxies if i.central == 1])
    met = np.asarray([i.metallicities['sfr_weighted'] for i in sim.galaxies if i.central == 1])
    sigmastar = np.asarray([i.velocity_dispersions['stellar'] for i in sim.galaxies if i.central == 1])
    ms_sat = np.asarray([i.masses['stellar'] for i in sim.galaxies if i.central != 1])
    mbh_sat = np.asarray([i.masses['bh'] for i in sim.galaxies if i.central != 1])
    sigmastar_sat = np.asarray([i.velocity_dispersions['stellar'] for i in sim.galaxies if i.central != 1])

    snapfile = '/home/cossim/HELUCID/%s/snapdir_%03d/snap_%03d' % (MODEL,SNAP, SNAP)
    sv = readsnap(snapfile,'vel','star',units=1,suppress=1) # physical km/s
    cents = np.asarray([i for i in sim.galaxies if i.central == 1])
    sigv1d = []
    for g in cents:
        svgal = np.array([sv[k] for k in g.slist])
        svcent = np.mean(svgal,axis=0)
        svgal = (svgal-svcent)*(svgal-svcent)
        sigv = np.sqrt(sum(svgal)/len(svgal))
        sigv1d.append(np.sqrt((sigv[0]*sigv[0]+sigv[1]*sigv[1]+sigv[2]*sigv[2])/3))
        #print np.log10(g.masses['stellar']),sigv,sigv1d[:1],g.velocity_dispersions['stellar']
    sats = np.asarray([i for i in sim.galaxies if i.central != 1])
    sigv1d_sat = []
    for g in sats:
        svgal = np.array([sv[k] for k in g.slist])
        svcent = np.mean(svgal,axis=0)
        svgal = (svgal-svcent)*(svgal-svcent)
        sigv = np.sqrt(sum(svgal)/len(svgal))
        sigv1d_sat.append(np.sqrt((sigv[0]*sigv[0]+sigv[1]*sigv[1]+sigv[2]*sigv[2])/3))

    logms = np.log10(ms)
    logmbh = np.log10(mbh)
    logsig = np.log10(sigv1d)
    logms_sat = np.log10(ms_sat)
    logmbh_sat = np.log10(mbh_sat)
    logsig_sat = np.log10(sigv1d_sat)
    ssfr = 1.e9*sfr/ms
    ssfr = np.log10(ssfr+10**(-2.5+0.3*redshift))
    pixcolor = ssfr
    pixsize = 8*(np.log10(ms/min(ms))+1)
    if iwind==0: 
        ax.plot(logsig_sat,logmbh_sat, 'x', c='grey', ms=1, label='Sats')
        im = ax.scatter(logsig,logmbh, c=pixcolor, s=pixsize, lw=0, cmap=plt.cm.jet_r, label='Centrals')
        fig.colorbar(im,ax=ax,label=r'sSFR')
        im.set_clim(-2.6+0.3*redshift,0.0+0.3*redshift)
    else: im = ax.scatter(logsig,logmbh, c='k', s=pixsize, lw=0)

plot_data(0)
#plot_data(0,'.')

plt.annotate('z=%g,%s'%(np.round(redshift,1),WIND[0]), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))

plt.minorticks_on()
plt.xlim(1.8,2.7)
plt.ylim(6,11)
plt.ylabel(r'$\log\ M_{BH}$',fontsize=16)
plt.xlabel(r'$\log\ \sigma_{*,1D}$' ,fontsize=16)
plt.legend(loc='lower right')

plt.savefig('./test_plot/mbhsigma_%s.pdf'%WIND[0], bbox_inches='tight', format='pdf')

plt.show()

