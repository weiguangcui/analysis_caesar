
import caesar
import sys
import pylab as plt
import numpy as np

if len(sys.argv) < 3:
    print('usage: MODEL SNAP WIND1 WIND2 ...')
    exit()

colorvar = 'ssfr'
#colorvar = 'mbhdot'

# define input file
MODEL = sys.argv[1]
SNAP = int(sys.argv[2])
WIND = sys.argv[3:]
#plt.rc('text', usetex=True)
mlim = 8.7
mmax = 12.9

modeldict = {'fh_qr':'Mufasa','s49':'Simba','s50nojet':'No-jet','s48':'Old Simba','s50nox':'No-Xray','s50':'Simba','s50fedd':'Simba','s50j7k':'Simba'}

def plot_data(redshift):
    if redshift < 0.5:
        lms = np.arange(9.5,mmax,0.3)
    #    lmbh = 7.45 + 1.05*(lms-11) # Reines & Volonteri 2015
    #    plt.plot(lms,lmbh,'--',c='c',label='RV15')
        lmbh = 8.69 + 1.16*(lms-11)	# Kormendy & Ho 2013
        plt.plot(lms,lmbh,'--',c='m',label='KH13')
        lmbh = 8.20 + 1.12*(lms-11) # Haring & Rix 2004
        #plt.plot(lms,lmbh,'--',c='k',label='HR04')
        lmbh = 8.18 + 1.00*(lms-11) # MBH=0.0015*M* (Schramm+Silverman 2013)
        #plt.plot(lms,lmbh,'--',c='k',lw=2,label=r'$M_{BH}/M_*=0.0015$')
        lms = np.arange(9,11,0.3)
        lmbh = 8.40 + 1.84*(lms-11) # Bentz+18
        plt.plot(lms,lmbh,':',c='k',lw=3,label=r'Bentz+18')
    if redshift > 1.5 and redshift < 3.5:
        lms = np.arange(mlim,12,0.3)
        lmbh = 7.45 + 1.05*(lms-11)
        plt.plot(lms,lmbh,'--',c='m',label='RV15')


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
    mdot = np.asarray([i.bhmdot for i in sim.galaxies if i.central == 1])
    sfr = np.asarray([i.sfr for i in sim.galaxies if i.central == 1])
    met = np.asarray([i.metallicities['sfr_weighted'] for i in sim.galaxies if i.central == 1])
    ms_sat = np.asarray([i.masses['stellar'] for i in sim.galaxies if i.central != 1])
    mbh_sat = np.asarray([i.masses['bh'] for i in sim.galaxies if i.central != 1])
    mdot_sat = np.asarray([i.bhmdot for i in sim.galaxies if i.central != 1])
    #print ms,mbh
    
    logms = np.log10(ms)
    logmbh = np.log10(mbh)
    logms_sat = np.log10(ms_sat)
    logmbh_sat = np.log10(mbh_sat)
    Zmet = np.log10(met/0.0189+1.e-6)
    ssfr = 1.e9*sfr/ms
    ssfr = np.log10(ssfr+10**(-2.5+0.3*redshift))
    if colorvar == 'ssfr': pixcolor = ssfr
    elif colorvar == 'mbhdot': pixcolor = np.log10(mdot+1.e-10)
    pixsize = 2*(np.log10(ms/min(ms))+1)
    if iwind==0: 
        plt.plot(logms[sfr>400],logmbh[sfr>400], 'x', c='k', ms=8, label='Simba SFR>400')
        im = ax.scatter(logms,logmbh, c=pixcolor, s=pixsize, lw=0, cmap=plt.cm.jet_r, label=WIND[iwind])
        #fig.colorbar(im,ax=ax,label=r'$\log\ M_h$')
        if colorvar == 'ssfr': 
            fig.colorbar(im,ax=ax,label=r'sSFR')
            im.set_clim(-2.9+0.3*redshift,0.2+0.3*redshift)
        elif colorvar == 'mbhdot': 
            fig.colorbar(im,ax=ax,label=r'$\dot{M}_{BH}$')
            im.set_clim(-3,0.5)
        #ax.plot(logms_sat,logmbh_sat, 'o', c='grey', ms=0.5, label='Satellites')
    else: im = ax.scatter(logms,logmbh, c='k', s=pixsize, lw=0)
    #plt.plot(logms,Zmet,ms=1,lw=0,c='r')
    #print logms,Zmet

plot_data(0)
#plot_data(0,'.')

plt.annotate('z=%g,%s'%(np.round(redshift,1),WIND[0]), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))

plt.minorticks_on()
plt.xlim(mlim,mmax)
plt.ylim(4.5,)
plt.ylabel(r'$\log\ M_{BH}$',fontsize=16)
plt.xlabel(r'$\log\ M_{*}$' ,fontsize=16)
plt.legend(loc='lower right')

plt.savefig('./test_plot/mbhms_%s.pdf'%(WIND[0]), bbox_inches='tight', format='pdf')

plt.show()

