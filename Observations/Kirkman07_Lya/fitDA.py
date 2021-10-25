
import numpy as np


z_K07,DA_K07,sig_K07 = np.loadtxt('Kirkman_07.dat',unpack=True)
z_D16,DA_D16,sig_D16 = np.loadtxt('Danforth_16.dat',unpack=True)
select_K07 = (z_K07>=0) & (z_K07<5.01)
D16_step = 1  # use every N'th value for D16
z = np.concatenate((z_D16[::D16_step],z_K07[select_K07]))
DA = np.concatenate((DA_D16[::D16_step],DA_K07[select_K07]))
sig = np.concatenate((sig_D16[::D16_step],sig_K07[select_K07]))
logz1 = np.log10(1.+z)
DA = np.log10(DA)

select = (z>-1) & (z<100)
#m,b = np.polyfit(logz1[select],DA[select],1,w=1./sig[select])
#m,b = np.polyfit(logz1[select],DA[select],1,w=1./np.sqrt(sig[select]))
m,b = np.polyfit(logz1[select],DA[select],1)
DA_fit = m*logz1[select]+b
print('slope=',m)
print('intercept=',b)
print('DA(z=0,z=1,z=2)=',np.interp(0,logz1,DA_fit),np.interp(np.log10(1+1),logz1,DA_fit),np.interp(np.log10(1+2),logz1,DA_fit))
diff = m*logz1[select]+b-DA[select]
print('sig=',np.sqrt(np.sum(diff*diff)/len(diff)))
print('diff=',diff)
print('fit=',m*logz1+b)
print('Danforth=',np.log10(0.014)+2.2*logz1)
print('vs Danforth=',m*logz1+b-np.log10(0.014)-2.2*logz1)

