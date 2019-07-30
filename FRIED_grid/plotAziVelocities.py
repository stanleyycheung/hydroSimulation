import numpy as np


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pylab

#dataFin100_1p3=np.loadtxt("1p3Msol_100AU_1000G0_20pc_finish.dat")
#dataAbund100_1p3=np.loadtxt("1p3Msol_100AU_1000G0_20pc_abundances.dat")
data=np.loadtxt("radial0019.dat")

au=1.5e3
mu=1.3
mH=1.67e-24
kerg=1.380626e-16

plt.axis(fontsize=24)
fig, axes= plt.subplots(figsize=(11,7))
axes.set_xlabel('R, AU', fontsize=18)
axes.set_ylabel('kms$^{-1}$', fontsize=18)

axes.tick_params(labelsize=17)
#axes2.tick_params(labelsize=16)

bigG=6.67259e-8
Mstar=1.989e33
axes.set_ylim(0,5)
#axes2=axes.twinx()
#axes2.set_ylabel('Temperature', fontsize=18)

#axes.plot(dataFin50[:,0]/au, np.log10(dataFin50[:,2]/mu/mH))
#axes.plot(dataFin100[:,0]/au, np.log10(dataFin100[:,2]/mu/mH))
axes.plot(data[:,0]/au, (data[:,5]), linewidth=4, label="Azimuthal velocity")
axes.plot(data[:,0]/au, (bigG*Mstar/(data[:,0]*1.e10))**0.5/1.e5, linewidth=4, label="Keplerian velocity", color="red")

#axes.plot(data[:,0]/au, ((kerg*data[:,2]/mu/mH)**0.5)/1.e5, linewidth=4, color="red", label="sound speed")

axes.legend(loc=0)


#axes2.plot(dataFin100[:,0]/au,(dataFin100[:,3]), color="red")
#axes.plot(dataFin200[:,0]/au, np.log10(dataFin200[:,2]/mu/mH))

plt.savefig("velocities_azi.pdf")
plt.show()
