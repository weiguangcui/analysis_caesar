import caesar
import sys
import os
import pylab as plt
import matplotlib.colors as colors
import numpy as np
import plotmedian as pm

# define input file
MODEL = sys.argv[1]
SNAP = int(sys.argv[2])
WIND = sys.argv[3:]

clrs = ('navy', 'g', 'crimson', 'm', 'r', 'y', 'k')

mmin = 5.8e8
mmin = 1.e10

if WIND == 'fh_qr': windlabel = 'Mufasa'
else: windlabel = 'Simba-SF'

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def ssfr_Illustris(zbin):
    infile = 'Observations/Illustris/data/sfr_vs_stellar_mass_z%g.txt' % (zbin)
    if os.path.isfile(infile):
        m,sfr = np.loadtxt(infile,unpack=True)
        plt.plot(np.log10(m),np.log10(sfr/m)+9,'--',label='Illustris',color='orange',lw=2)

def ssfr_eagle(zbin):
    zbin = np.round(zbin,0)
    infile = 'Observations/eagle/sfr/Ref_sfr_z%gp0.txt' % (zbin)
    if zbin == 0:
        infile = 'Observations/eagle/sfr/Ref_sfr_z%gp1.txt' % (zbin)
    if os.path.isfile(infile):
        m,sfr = np.loadtxt(infile,usecols=(0,1),unpack=True)
        plt.plot(m,sfr-m+9,':',label='EAGLE',color='m',lw=3)

def ssfr_data(redshift):
    mslo = 8.7
    ms_data = np.arange(mslo,11.5,0.1)
    if redshift<=0.5:
        ms_data = np.arange(8.5,11.3,0.1)
        ssfr_data = 5.96e-2 * 10**(-0.35*(ms_data-11.03)) * np.exp(-10**(ms_data-11.03)) #Salim+07
        plt.plot(ms_data,np.log10(ssfr_data),'-',lw=4,label='Salim+07',color='k')
#        ssfr_data = (0.84-0.026*tage)*ms_data-(6.51-0.11*tage)-ms_data+9  # Speagle+14
#        plt.plot(ms_data,ssfr_data,'--',lw=4,label='Speagle+14',color='k')
        ms_data = np.arange(9,11,0.1)
        ssfr_data = 0.8*(ms_data-10)-0.23 -ms_data+9    #Chang+15
        #plt.plot(ms_data,ssfr_data,'--',lw=4,label='Chang+15',color='k')
        ms_data = [9,9.5,10,10.5,11]
        ssfr_data = np.asarray([-0.653,-0.814,-0.822,-1.056,-1.457])    # Bauer+13 SFG
        #ssfr_data = np.asarray([-0.7,-0.8,-0.8,-1.0,-1.4])
        #ssfr_data = np.asarray([-0.4,-0.8,-1.0,-1.4,-1.6])
        #ssfr_data = ssfr_data - 0.15   # correct from z~0.15 to z=0 based on (1+z)^2 evolution
        plt.plot(ms_data,ssfr_data,'o',ms=8,mfc='gray',label='Bauer+13',color='k')
    elif redshift>0.5 and redshift<1.5:
        ssfr_lo = -27.40+5.02*ms_data-0.22*ms_data*ms_data  # Whitaker+14
        ssfr_hi = -26.03+4.62*ms_data-0.19*ms_data*ms_data
        ssfr_data = np.log10(0.5*(10**ssfr_lo+10**ssfr_hi))-ms_data+9
        plt.plot(ms_data,ssfr_data,'-',lw=4,label='Whitaker+14',color='k')
    elif redshift>=1.5 and redshift <2.5:
        ssfr_lo = -24.04+4.17*ms_data-0.16*ms_data*ms_data  # Whitaker+14
        ssfr_hi = -19.99+3.44*ms_data-0.13*ms_data*ms_data
        ssfr_lo = ssfr_hi  # use 2-2.5 since something seems wrong with 1.5-2 fit; doesnt match
        ssfr_data = np.log10(0.5*(10**ssfr_lo+10**ssfr_hi))-ms_data+9
        #print ms_data,ssfr_lo-ms_data+9,ssfr_hi-ms_data+9,ssfr_data
        plt.plot(ms_data,ssfr_data,'--',lw=3,label='Whitaker+14',color='k')
    elif redshift>=2.5 and redshift <4.5:
        ms_data = np.asarray([9.00,9.25,9.50,9.75,10.00,10.25,10.5])
        ms_data= ms_data-0.25 # to Chabrier IMF
        ssfr_data = np.asarray([0.71,0.90,1.01,1.04,1.35,1.51,1.87])
        ssfr_sig = [0.36,0.41,0.35,0.24,0.27,0.25,0.24]
        ssfr_data = ssfr_data-ms_data+9
        plt.plot(ms_data,ssfr_data,'o',ms=8,mfc='gray',label='Salmon+14',color='k')
        plt.errorbar(ms_data,ssfr_data,yerr=[ssfr_sig,ssfr_sig],color='k')
    elif redshift>=4.5 and redshift <5.5:
        ms_data = np.asarray([9.00,9.25,9.50,9.75,10.00,10.25,10.5])
        ms_data= ms_data-0.25 # to Chabrier IMF
        ssfr_data = np.asarray([0.88,1.04,1.12,1.23,1.46,1.62,1.85])
        ssfr_sig = [0.42,0.38,0.41,0.43,0.31,0.37,0.33]
        ssfr_data = ssfr_data-ms_data+9
        plt.plot(ms_data,ssfr_data,'o',ms=8,mfc='gray',label='Salmon+14',color='k')
        plt.errorbar(ms_data,ssfr_data,yerr=[ssfr_sig,ssfr_sig],color='k')
    elif redshift>=5.5 and redshift <6.5:
        ms_data = np.asarray([9.00,9.25,9.50,9.75,10.00,10.25])
        ms_data= ms_data-0.25 # to Chabrier IMF
        ssfr_data = np.asarray([0.92,1.07,1.27,1.40,1.47,1.79])
        ssfr_sig = np.asarray([0.19,0.21,0.35,0.26,0.07,0.35])
        ssfr_data = ssfr_data-ms_data+9
        plt.plot(ms_data,ssfr_data,'o',ms=8,mfc='gray',label='Salmon+14',color='k')
        plt.errorbar(ms_data,ssfr_data,yerr=[ssfr_sig,ssfr_sig],color='k')
    else:
        ssfr_data = (0.84-0.026*tage)*ms_data-(6.51-0.11*tage)-ms_data+9  # Speagle+14
        plt.plot(ms_data,ssfr_data,'-',lw=4,label='Speagle+14',color='k')

def GSWLC_data(ssfrlim=-2):
    ms,sigms,sfr,sigsfr,flag_sed,flag_mgs = np.loadtxt("Observations/GSWLC/GSWLC-X2.dat",usecols=(9,10,11,12,19,23),unpack=True)
    ssfr = sfr-ms+9
    select = (flag_sed==0)&(flag_mgs==1)
    plt.hexbin(ms[select],ssfr[select],gridsize=50,cmap=plt.cm.Greys)
    select = (flag_sed==0)&(flag_mgs==1)&(ssfr>ssfrlim)
    pm.plotmedian(ms[select],ssfr[select],ax=plt,bins=10,c='k',stat='median',label='GSWLC')

# load in input file
fig,ax = plt.subplots()
for iwind in range(len(WIND)):
    infile = '/home/cossim/HELUCID/%s/%s/Caesar_snap_%03d.hdf5' % (MODEL,WIND[iwind],SNAP)
    sim = caesar.load(infile)
    redshift = sim.simulation.redshift
    volume = sim.simulation.boxsize.to('Mpccm').d**3

    pos = np.asarray([i.pos for i in sim.galaxies])
    ids = np.asarray([i.GroupID for i in sim.galaxies if i.central == 1])
    ms = np.asarray([i.masses['stellar'] for i in sim.galaxies if i.central == 1])
    mbh = np.asarray([i.masses['bh'] for i in sim.galaxies if i.central == 1])
    mdot = np.asarray([i.bhmdot*1.99e33/3.1416e7 for i in sim.galaxies if i.central == 1])
    mHI = np.asarray([i.masses['HI'] for i in sim.galaxies if i.central == 1])
    sfr = np.asarray([i.sfr for i in sim.galaxies if i.central == 1])
    met = np.asarray([i.metallicities['sfr_weighted'] for i in sim.galaxies if i.central == 1])

    #mh = np.asarray([i.halo.masses['total'] for i in sim.galaxies if i.central == 1])
    ms_sat = np.asarray([i.masses['stellar'] for i in sim.galaxies if i.central != 1])
    mbh_sat = np.asarray([i.masses['bh'] for i in sim.galaxies if i.central != 1])
    mdot_sat = np.asarray([i.bhmdot*1.99e33/3.1416e7 for i in sim.galaxies if i.central != 1])
    #mh_sat = np.asarray([i.halo.masses['total'] for i in sim.galaxies if i.central != 1])
    mHI_sat = np.asarray([i.masses['HI'] for i in sim.galaxies if i.central != 1])
    sfr_sat = np.asarray([i.sfr for i in sim.galaxies if i.central != 1])
    met_sat = np.asarray([i.metallicities['sfr_weighted'] for i in sim.galaxies if i.central != 1])

    ssfr = np.log10(1.e9*sfr/ms+10**(-2.7+0.3*redshift))
    ssfr_sat = np.log10(1.e9*sfr_sat/ms_sat+10**(-2.7+0.3*redshift))

    LX = 0.1*mdot*9.e20+1.e40
    LX_sat = 0.1*mdot_sat*9.e20+1.e40
    Sey = 6.e43

    sel = (ms>1.e10)&(ssfr>-1.8)
#     print 1.*len(ms[sel&(LX>Sey)])/len(ms[sel])
    sel = (ms>1.e11)&(ssfr>-1.8)
#     print 1.*len(ms[sel&(LX>Sey)])/len(ms[sel])

    colorvec = np.log10(mHI/ms+1.e-3)
    colorvec = np.log10(met[met>0])
    colorvec = np.log10(mbh/ms+1.e-5)
    colorvec = np.log10(LX)
    colorvec_sat = np.log10(LX_sat)
    print(colorvec)
    logms = np.log10(ms)
    massbin,cvecbin,ebinlo,ebinhi = pm.runningmedian(logms[ms>1.e10],colorvec[ms>1.e10],bins=12)
    cvec = colorvec - np.interp(logms,massbin,cvecbin)
    cvec_sat = colorvec_sat - np.interp(np.log10(ms_sat),massbin,cvecbin)

    print(redshift,len(ssfr[ssfr<-0.82])/volume)

    sfrlim = -1.8+0.3*redshift

    if iwind ==0:
        cmap = plt.get_cmap('jet')
        new_cmap = truncate_colormap(cmap, 0.2, 0.8)

        pixcolor = colorvec
        pixsize = 1*np.log10(ms/min(ms))+1
        im = ax.scatter(np.log10(ms),ssfr, c=pixcolor, s=pixsize, lw=0, alpha=1.0, cmap=new_cmap)
        #fig.colorbar(im,ax=ax,label=r'$\log\ M_{BH}/M_*$')
        fig.colorbar(im,ax=ax,label=r'$\log\ L_{X}$')
        #cv_med = -2.7
        #im.set_clim(cv_med-0.3,cv_med+0.3)
        im.set_clim(43,46)

        pixcolor = colorvec_sat
        im = ax.scatter(np.log10(ms_sat),ssfr_sat, c=pixcolor, s=0.5, lw=0, alpha=1.0, cmap=new_cmap)

    ms = np.append(ms,ms_sat)
    ssfr = np.append(ssfr,ssfr_sat)+0.01*np.random.rand(len(ms))
    #sfrlim = -1.5-0.5*(np.log10(ms)-9)+0.3*redshift
    plt.plot(np.log10(ms),sfrlim*np.ones(len(ms)),':',c=clrs[-1-iwind])
    pm.plotmedian(np.log10(ms[ssfr>sfrlim]),ssfr[ssfr>sfrlim],ax=ax,bins=12,c=clrs[iwind],ltype='-',lw=3,label=WIND[iwind],stat='median',pos=pos[ssfr>sfrlim]) #,boxsize=sim.simulation.boxsize)

if redshift<0.3: GSWLC_data(sfrlim)
else: ssfr_data(redshift)
    
ssfr_eagle(redshift)

plt.minorticks_on()
plt.xlim(np.log10(mmin),max(np.log10(ms))+0.01)
plt.ylim(min(ssfr)-0.1,0.4+0.3*redshift)
plt.ylabel(r'$\log\ sSFR$',fontsize=16)
plt.xlabel(r'$\log\ M_{*}$' ,fontsize=16)
plt.legend(loc='upper right')

plt.annotate('z=%g'%(np.round(redshift,1)), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))

plt.savefig('test_plot/mainseq_z%g.pdf'%np.round(redshift,0), bbox_inches='tight',dvi=300)
plt.show()

