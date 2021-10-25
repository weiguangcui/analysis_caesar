
import caesar
import sys
import pylab as plt
import numpy as np

# define input file
MODEL = sys.argv[1]
WIND = sys.argv[2]
SNAP = int(sys.argv[3])
# infile = '/home/wgcui/%s/%s/CaesarGroups_%03d.hdf5' % (MODEL,WIND,SNAP)
infile = '/home/cossim/HELUCID/%s/%s/Caesar_snap_%03d.hdf5' % (MODEL,WIND,SNAP)

def GSWLC_data(mbin,ssfrlim=-2.5):
    ms,sigms,sfr,sigsfr,flag_sed,flag_mgs = np.loadtxt("Observations/GSWLC/GSWLC-D2.dat",usecols=(9,10,11,12,19,23),unpack=True)
    ssfr = sfr-ms+9
    ssfr = np.maximum(ssfr,ssfrlim)
    select = (flag_sed==0)&(flag_mgs==1)
    ssfrbin = np.arange(-3,1,0.2)
    for i in range(0,len(mbin)):
        ssfrp = ssfr[(ms>mbin[i][0]) & (ms<mbin[i][1])]
        hist,bin_edges = np.histogram(ssfrp,ssfrbin,normed=True)
        hist /= np.sum(hist)
        bin_cent = 0.5*(bin_edges[:-1]+bin_edges[1:])
        plt.plot(bin_cent,hist,':',lw=2,c=colors[i])

colors = ('navy', 'g', 'crimson', 'm', 'r', 'y', 'k')

# load in input file
sim = caesar.load(infile)
redshift = sim.simulation.redshift

ids = np.asarray([i.GroupID for i in sim.galaxies if i.central <= 1])
ms = np.asarray([i.masses['stellar'] for i in sim.galaxies if i.central <= 1])
sfr = np.asarray([i.sfr for i in sim.galaxies if i.central <= 1])
#mh = np.asarray([i.halo.masses['total'] for i in sim.galaxies if i.central == 1])
ms_sat = np.asarray([i.masses['stellar'] for i in sim.galaxies if i.central != 1])

#logmh = np.log10(mh)
ssfrlim = -2.5+0.3*redshift
logms = np.log10(ms)
#msmh = np.log10(ms/mh)
ssfr = 1.e9*sfr/ms
ssfr = np.log10(ssfr+10**ssfrlim)
pixcolor = ssfr
pixsize = 5*(np.log10(ms/min(ms))+1)
fig,ax = plt.subplots()
mbin = [[9,10],[10,11],[11,15]]
ssfrbin = np.arange(-3,0.5,0.2)
for i in range(0,len(mbin)):
    ssfrp = ssfr[(logms>mbin[i][0]) & (logms<mbin[i][1])]
    sfrp = np.log10(sfr[(logms>mbin[i][0]) & (logms<mbin[i][1])]+1.e-2)
    hist,bin_edges = np.histogram(ssfrp,ssfrbin,normed=True)
    hist /= np.sum(hist)
    bin_cent = 0.5*(bin_edges[:-1]+bin_edges[1:])
    if mbin[i][1]<13: label = r'$10^{%.4g}<M_*<10^{%.4g}$' % (mbin[i][0],mbin[i][1])
    else: label = r'$M_*>10^{%.4g}$' % (mbin[i][0])
    plt.plot(bin_cent,hist,'-',lw=2,c=colors[i],label=label)

GSWLC_data(mbin,ssfrlim)

plt.minorticks_on()
plt.xlabel(r'$\log$ sSFR (Gyr$^{-1}$)',fontsize=16)
plt.ylabel(r'fraction' ,fontsize=16)
plt.annotate('Solid: Simba, z=%g'%(np.round(redshift,1)), xy=(0.95, 0.7), xycoords='axes fraction',size=14,horizontalalignment='right') #,bbox=dict(boxstyle="round", fc="w"))
plt.annotate('Dotted: GSWLC-D2', xy=(0.95, 0.6), xycoords='axes fraction',size=14,horizontalalignment='right') #,bbox=dict(boxstyle="round", fc="w"))
plt.legend(loc='upper right')

plt.savefig('test_plot/mainseqhist_%s.pdf'%WIND)
plt.show()

