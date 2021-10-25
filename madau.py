from readgadget import *
from pylab import *
import sys
import os.path
from astropy.cosmology import FlatLambdaCDM

colors = ('g', 'c', 'crimson', 'navy', 'm', 'y', 'k')

xvar = str(sys.argv[1])	# 'z' plots vs. redshift, 'a' vs. expansion factor, 't' vs. time


pref = sys.argv[2:]
binsize = 1000
header=readheader('/home/cossim/HELUCID/'+pref[0]+'/snapdir_000/snap_000.0.hdf5','header')
boxstr = [header['boxsize']/1000.*header['npartTotal'][0]**(1/3.)/5400]
zlim = [0]
boxsize = np.asarray(boxstr)

cosmo = FlatLambdaCDM(H0=header['h']*100, Om0=header['O0'])
#rc('text', usetex=True)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

for k in range(0,len(boxsize)):
    for i in range(0,len(pref)):
#         thislabel = 'm'+str(boxstr[k])+'n'+nside+'/'+pref[i]
#         modlabel = 'm'+str(boxstr[k])+'n'+nside
#         modlabel = 'Simba'
        s = '/home/cossim/HELUCID/'+pref[i]+'/sfr.txt'
        if os.path.isfile(s):
            print(s,boxsize[k])
        else:
            raise ValueError('Can not find file: ', s)
        aex,mformed,sfrsys,sfr,mstar = loadtxt(s,unpack=True)
        a = []
        sfrate = []
        for j in range(0,len(aex)-len(aex)%binsize,binsize):
            tmp = sum(aex[j:j+binsize]/binsize)
            a.append(tmp)
            tmp = sum(sfr[j:j+binsize]/binsize)
            sfrate.append(tmp)
        alim = 1./(1+zlim[k])
        a = np.asarray(a)
        sfrate = np.asarray(sfrate)
        vol = 0.7*0.7*0.7/(boxsize[k]*boxsize[k]*boxsize[k])
        sfrate = vol*sfrate[(a>0.1) & (a<alim)]
        z = [1./x-1 for x in a if x>0.1 and x<alim]
        aexp = [x for x in a if x>0.1 and x<alim]
        z = np.asarray(z)
        t = cosmo.age(z).value
        lz = np.log10(1+z)
        ax1.set_ylabel(r'$\log\ SFRD (M_\odot yr^{-1} Mpc^{-3})$',fontsize=20)
        if xvar == 'a':
            ax1.set_xlabel('Expansion factor',fontsize=20)
            ax1.plot(aexp,log10(sfrate),label=pref[i],color=colors[i])
        elif xvar == 'z':
            ax1.set_xlabel('Redshift',fontsize=20)
            ax1.plot(z,log10(sfrate),label=pref[i],color=colors[i])
        elif xvar == 't':
            ax1.plot(t,log10(sfrate),label=pref[i],color=colors[i])
            ax1.set_xlabel('Cosmic Time (Gyr)',fontsize=20)
            topticks1 = np.array([0,1,2,4,6])  # desired redshift labels
            topticks2 = cosmo.age(topticks1).value  # tick locations in time
            ax2.set_xticks(topticks2)
            ax2.set_xticklabels(topticks1)
            ax2.set_xlabel('Redshift',fontsize=20)
        elif xvar == 'lz':
            ax1.plot(lz,log10(sfrate),label=pref[i],color=colors[i])
            ax1.set_xlabel(r'$\log(1+z)$',fontsize=20)
            topticks1 = np.array([0,1,2,4])  # desired redshift labels
            topticks2 = np.log10(1+topticks1)  # tick locations in log(1+z)
            print(max(lz),topticks2)
            ax2.set_xticks(topticks2)
            ax2.set_xticklabels(topticks1)
            ax2.set_xlabel('Redshift',fontsize=20)
        else:
            sys.exit('file %s not found'%s)

# Plot fit from Madau & Dickinson 2014
z = np.arange(0,8,0.01)
a = 1./(1+z)
lz = np.log10(1+z)
t = cosmo.age(z).value
imf_factor = 1.7 # Salp->Chabrier
sfrd_md14 = 0.015*(1+z)**2.7/(1+((1+z)/2.9)**5.6) / imf_factor 
ax1.plot(eval(xvar),log10(sfrd_md14),color='k')
z,lsfrd,ehi,elo=np.loadtxt("Observations/md14_data.txt",unpack=True)
a = 1./(1+z)
lz = np.log10(1+z)
t = cosmo.age(z).value
lsfrd = lsfrd-np.log10(imf_factor)
ax1.plot(eval(xvar),lsfrd,'o',color='k',ms=3,label='MD14')
ax1.errorbar(eval(xvar),lsfrd,fmt='o',color='k',yerr=[-elo,ehi])

#savefig("madau.eps")
if xvar == 'a':
    ax1.set_xlim(0.1,1)
elif xvar == 'z':
    ax1.set_xlim(min(zlim),8.5)
elif xvar == 't':
    ax1.set_xlim(0,max(cosmo.age(zlim).value))
elif xvar == 'lz':
    maxlz = 0.8
    ax1.set_xlim(0,maxlz)
    ax2.set_xlim(0,maxlz)
# ax1.set_ylim(-2.6,-0.8)
ax1.grid()

if xvar == 'z':
    ax1.legend(loc='upper right',fontsize=12)
else:
    ax1.legend(loc='lower center',fontsize=12)
ax1.tick_params(labelsize=16)
ax2.tick_params(labelsize=16)

figname = 'test_plot/madau_'+xvar+'.pdf'
savefig(figname,bbox_inches='tight')
show()
